import os
import gzip
import io
import regex
import numpy as np
import faiss
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict, Counter
from datasketch import MinHash
import math
import multiprocessing
from joblib import Parallel, delayed
from operator import itemgetter
import random

__all__ = [
    'zipped',
    'ensure_open',
    'open_file', 
    'natural_key',
    'check_genome_compression',
    'fetch_sequence',
    'set_faiss_threads',
    'apply_minhash',
    'get_nn',
    'group_similar_hashes',
    'extract_grouped_entries',
    'TRANS_TABLE',
    'revcomp',
    'collapse_group',
    'final_pass_collapse',
    'sample_reads',
    'check_consensus',
    'generate_consensus'
]

TRANS_TABLE = bytes.maketrans(b"acgtumrwsykvhdbnACGTUMRWSYKVHDBN", 
                                b"TGCAAKYWSRMBDHVNTGCAAKYWSRMBDHVN")

def revcomp(seq):
    return seq.translate(TRANS_TABLE)[::-1]

def zipped(file):
    if os.path.exists(file):
        with open(file, 'rb') as zip_test:
            return zip_test.read(2) == b'\x1f\x8b'
    else:
        return file.endswith('.gz')
    
def ensure_open(file):
    if isinstance(file, (str, bytes)):
        if zipped(file):
            return gzip.open(file, 'rt')
        else:
            return open(file, 'rt')
    elif isinstance(file, io.IOBase) and not file.closed:
        return file
    else:
        raise ValueError("Invalid file input or file is closed")
    
def open_file(filename, mode):
    if zipped(filename):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
    
def natural_key(chrom):
    pattern = regex.compile(r'(\d+)')
    return [int(text) if text.isdigit() else text for text in pattern.split(chrom)]

def check_genome_compression(genome_fasta):
    with open(genome_fasta, 'rb') as f:
        magic_number = f.read(4)
        if magic_number == b'\x1f\x8b\x08\x00':
            return 'gzip'
        elif magic_number == b'\x1f\x8b\x08\x04':
            return 'bgzf'
    return 'uncompressed'

def fetch_sequence(coordinates, genome, U3, shift):
    chrom = coordinates['chrom']
    if coordinates['strand'] == '-':
        strand = '-' if not U3 else '+'
    else:
        strand = '+' if not U3 else '-'
        
    chr_len = genome.get_reference_length(chrom)
    
    if strand == '-':
        start = max(((coordinates['end'] - shift) - 40), 0)
        end = min((coordinates['end'] - shift) + 15, chr_len)
    else:
        start = max(((coordinates['start'] + shift) - 15), 0)
        end = min(((coordinates['start'] + shift) + 40), chr_len)
        
    sequence = genome.fetch(chrom, start, end)
    
    if strand == '-':
        sequence = revcomp(sequence)
    return sequence.upper()

def apply_minhash(key, value, string, num_perm, token_size):
    minhash = MinHash(num_perm = num_perm)
    tokens = [string[i:i + token_size] for i in range(len(string) - token_size + 1)]
    
    for token in tokens:
        minhash.update(token.encode('utf8'))
    
    minhash_values = np.array(minhash.hashvalues, dtype = np.uint32)
    
    packed_bits = minhash_values.view(np.uint8)
    
    value['hash'] = packed_bits
    
    return key, value

def set_faiss_threads(nthr):
    if hasattr(faiss, 'omp_set_num_threads'):
        faiss.omp_set_num_threads(nthr)
        return True
    return False

def get_nn(input_dict, num_perm, nthr, len_diff):
    set_faiss_threads(nthr)
    batch_size = min(1000, len(input_dict))
    nn_k = min(len(input_dict), 1000)
    dim = num_perm * 4
    
    db = np.zeros((len(input_dict), dim), dtype='uint8')
    for i, (_, value) in enumerate(input_dict.items()):
        db[i] = value['hash']
    
    hash_index = faiss.IndexBinaryFlat(dim * 8)
    hash_index.add(db)

    all_distances = []
    all_indices = []
    
    for start_idx in range(0, len(db), batch_size):
        end_idx = min(start_idx + batch_size, len(db))
        query_batch = db[start_idx:end_idx]
        
        distances, indices = hash_index.search(query_batch, k=nn_k)
        all_distances.append(distances)
        all_indices.append(indices)
    
    return np.vstack(all_distances), np.vstack(all_indices)

def group_similar_hashes(hash_distances, hash_indices, nthr, similarity, num_perm, input_dict):
    hamming_threshold = (num_perm * 32) * (1 - similarity)
    
    n_sequences = len(input_dict)
    
    rows, cols = [], []
    valid_queries = hash_indices[:, 0] < n_sequences
    
    for i in np.where(valid_queries)[0]:
        mask = hash_distances[i] <= hamming_threshold
        valid_idx = hash_indices[i][mask]
        valid_idx = valid_idx[valid_idx < n_sequences]
        if len(valid_idx) > 0:
            rows.extend([hash_indices[i, 0]] * len(valid_idx))
            cols.extend(valid_idx)
    
    graph = csr_matrix(
        ([1] * len(rows), (rows, cols)),
        shape=(n_sequences, n_sequences),
        dtype=np.int8
    )
    graph = graph.maximum(graph.T)
    
    n_components, labels = connected_components(
        csgraph=graph,
        directed=False,
        return_labels=True
    )
    
    # More efficient component extraction
    unique_entries = []
    grouped_entries = []
    
    # Use bincount for single-item components
    label_counts = np.bincount(labels)
    
    for component_id in range(n_components):
        if label_counts[component_id] == 1:
            # This is a single-item component
            unique_entries.append(int(np.where(labels == component_id)[0][0]))
        elif label_counts[component_id] > 1:
            # This is a multi-item component
            grouped_entries.append(sorted(np.where(labels == component_id)[0].astype(int).tolist()))
            
    return unique_entries, grouped_entries

def extract_grouped_entries(groups, input_dict):
    dict_values = list(input_dict.values())
    grouped_entries = []
    for group in groups:
        entries = [dict_values[index] for index in group]
        grouped_entries.append(entries)

    return grouped_entries

def collapse_group(chrom_strand_tuple, pos_counts, len_diff, min_count, 
                    count_fc, read_mapping):
    positions = np.array([p for p, _ in pos_counts])
    counts = np.array([c for _, c in pos_counts])
    
    abundant_mask = counts >= min_count
    abundant_positions = positions[abundant_mask]
    abundant_counts = counts[abundant_mask]
    
    sort_idx = np.argsort(abundant_counts)[::-1]
    abundant_positions = abundant_positions[sort_idx]
    abundant_counts = abundant_counts[sort_idx]
    
    keep_abundant = np.ones(len(abundant_positions), dtype=bool)
    for i in range(1, len(abundant_positions)):
        nearby_mask = np.abs(abundant_positions[:i] - abundant_positions[i]) <= len_diff
        if any(nearby_mask):
            max_nearby_count = np.max(abundant_counts[:i][nearby_mask])
            if max_nearby_count > abundant_counts[i] * count_fc:
                keep_abundant[i] = False
    
    final_abundant_positions = abundant_positions[keep_abundant]
    
    read_updates = {}
    for abundant_pos in final_abundant_positions:
        collapse_mask = np.abs(positions - abundant_pos) <= len_diff
        for pos in positions[collapse_mask]:
            for read_name in read_mapping[(chrom_strand_tuple[0], pos, chrom_strand_tuple[1])]:
                read_updates[read_name] = abundant_pos
    
    return read_updates

def final_pass_collapse(kept_frags, len_diff, nthr, min_count, count_fc):
    positions = np.array([(frag['chrom'], 
                            frag['end'] if frag['strand'] == '-' else frag['start'],
                            frag['strand'])
                            for frag in kept_frags.values()],
                            dtype=[('chrom', 'U25'), ('pos', 'i8'), ('strand', 'U1')])
    unique_pos, counts = np.unique(positions, return_counts=True)
    
    read_mapping = defaultdict(list)
    for read_name, frag in kept_frags.items():
        pos = frag['end'] if frag['strand'] == '-' else frag['start']
        read_mapping[(frag['chrom'], pos, frag['strand'])].append(read_name)
    
    grouped_data = defaultdict(list)
    for pos, count in zip(unique_pos, counts):
        grouped_data[(pos['chrom'], pos['strand'])].append((pos['pos'], count))
    
    results = Parallel(n_jobs=nthr)(
        delayed(collapse_group)(
            key, 
            pos, 
            min_count, 
            count_fc, 
            len_diff, 
            read_mapping)
        for key, pos in grouped_data.items()
    )
    
    for read_updates in results:
        for read_name, pos in read_updates.items():
            if kept_frags[read_name]['strand'] == '-':
                kept_frags[read_name]['end'] = pos
            else:
                kept_frags[read_name]['start'] = pos
    
    return kept_frags

def sample_reads(filename, n_reads=5000):
    sequences = []
    read_count = 0
    
    with open_file(filename, 'rt') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # Sequence lines
                if read_count < n_reads:
                    # Fill the reservoir until we have n_reads
                    sequences.append(line.strip())
                else:
                    # Use reservoir sampling for the rest
                    j = random.randint(0, read_count)
                    if j < n_reads:
                        sequences[j] = line.strip()
                
                read_count += 1
            
            line_num += 1
    
    return sequences

def check_consensus(r1_file, r2_file, ltr_seq, linker_seq, sample_size=5000, threshold=0.02, error_rate=0.3):
    print(f"Sampling {sample_size} reads from input files...", flush=True)
    r1_sequences = sample_reads(r1_file, sample_size)
    r2_sequences = sample_reads(r2_file, sample_size)
    
    # Compile patterns using the hierarchical approach from crop.py
    ltr_errors = math.floor(len(ltr_seq) * error_rate)
    linker_errors = math.floor(len(linker_seq) * error_rate)
    
    patterns = {
        'ltr': {
            'perfect': regex.compile(ltr_seq),
            'mismatch': regex.compile(f'({ltr_seq}){{s<={ltr_errors}}}'),
            'indel': regex.compile(f'({ltr_seq}){{e<={ltr_errors}}}')
        },
        'linker': {
            'perfect': regex.compile(linker_seq),
            'mismatch': regex.compile(f'({linker_seq}){{s<={linker_errors}}}'),
            'indel': regex.compile(f'({linker_seq}){{e<={linker_errors}}}')
        }
    }
    
    def find_pattern_match(seq, pattern_type):
        match = patterns[pattern_type]['perfect'].search(seq)
        if match:
            return match, 'perfect'
        
        match = patterns[pattern_type]['mismatch'].search(seq)
        if match:
            return match, 'mismatch'
        
        match = patterns[pattern_type]['indel'].search(seq)
        if match:
            return match, 'indel'
        
        return None, None
    
    # Extract matching sequences
    r1_matches = []
    r2_matches = []
    
    print("Searching for LTR/linker sequences in sampled reads...", flush=True)
    
    # Process R1 sequences (expected to contain LTR)
    for seq in r1_sequences:
        match, match_type = find_pattern_match(seq, 'ltr')
        if match:
            r1_matches.append({
                'sequence': seq,
                'match_start': match.start(),
                'match_end': match.end(),
                'matched_text': match.group(),
                'match_type': match_type
            })
    
    # Process R2 sequences (expected to contain linker)
    for seq in r2_sequences:
        match, match_type = find_pattern_match(seq, 'linker')
        if match:
            r2_matches.append({
                'sequence': seq,
                'match_start': match.start(),
                'match_end': match.end(),
                'matched_text': match.group(),
                'match_type': match_type
            })
    
    print(f'Found {len(r1_matches)} R1 sequences and {len(r2_matches)} R2 sequences', flush=True)
    
    if not r1_matches or not r2_matches:
        if not r1_matches and not r2_matches:
            print(f'No matching LTR or linker sequences found in the sampled reads. Check the given inputs.')
        if not r1_matches:
            print(f'No matching LTR sequences found in the sampled reads. Check the given input.')
        if not r2_matches:
            print(f'No matching linker sequences found in the sampled reads. Check the given input.')
        sys.exit()
    
    # Count sequence occurrences
    r1_counts = {}
    for match in r1_matches:
        seq = match['matched_text']
        r1_counts[seq] = r1_counts.get(seq, 0) + 1
    
    r2_counts = {}
    for match in r2_matches:
        seq = match['matched_text']
        r2_counts[seq] = r2_counts.get(seq, 0) + 1
    
    # Calculate frequencies
    r1_total = len(r1_matches)
    r2_total = len(r2_matches)
    
    r1_freqs = {seq: count/r1_total for seq, count in r1_counts.items()}
    r2_freqs = {seq: count/r2_total for seq, count in r2_counts.items()}
    
    # Format report
    report = []
    report.append("")
    
    report.append(f"Target LTR sequence: {ltr_seq}")
    report.append(f"Target Linker sequence: {linker_seq}")
    report.append(f"Sample size: {sample_size} reads")
    report.append(f"Consensus threshold: {threshold}")
    report.append(f"Sequences with LTR match: {len(r1_matches)}")
    report.append(f"Sequences with linker match: {len(r2_matches)}")
    report.append("")
    
    # Report match type statistics
    r1_match_types = {}
    for match in r1_matches:
        match_type = match['match_type']
        r1_match_types[match_type] = r1_match_types.get(match_type, 0) + 1
    
    r2_match_types = {}
    for match in r2_matches:
        match_type = match['match_type']
        r2_match_types[match_type] = r2_match_types.get(match_type, 0) + 1
    
    report.append("R1 Match Types:")
    for match_type in ['perfect', 'mismatch', 'indel']:
        if match_type in r1_match_types:
            count = r1_match_types[match_type]
            report.append(f"  {match_type}: {count} ({count/r1_total:.2%})")
        else:
            report.append(f"  {match_type}: 0 (0.00%)")

    report.append("")

    report.append("R2 Match Types:")
    for match_type in ['perfect', 'mismatch', 'indel']:
        if match_type in r2_match_types:
            count = r2_match_types[match_type]
            report.append(f"  {match_type}: {count} ({count/r2_total:.2%})")
        else:
            report.append(f"  {match_type}: 0 (0.00%)")
    
    report.append("")
    
    # Report R1 sequence frequencies
    report.append("R1 Sequence Frequencies:")
    r1_significant = {seq: freq for seq, freq in r1_freqs.items() if freq >= threshold}
    r1_sorted = sorted(r1_significant.items(), key=lambda x: x[1], reverse=True)
    
    r1_significant_total = sum(r1_significant.values())
    r1_other = 1 - r1_significant_total
    
    for seq, freq in r1_sorted:
        report.append(f"{freq:.4f}\t{seq}")
    
    if r1_other > 0:
        report.append(f"{r1_other:.4f}\tOther")
    
    report.append("")
    
    # Report R2 sequence frequencies
    report.append("R2 Sequence Frequencies:")
    r2_significant = {seq: freq for seq, freq in r2_freqs.items() if freq >= threshold}
    r2_sorted = sorted(r2_significant.items(), key=lambda x: x[1], reverse=True)
    
    r2_significant_total = sum(r2_significant.values())
    r2_other = 1 - r2_significant_total
    
    for seq, freq in r2_sorted:
        report.append(f"{freq:.4f}\t{seq}")
    
    if r2_other > 0:
        report.append(f"{r2_other:.4f}\tOther")
    
    return "\n".join(report)

def generate_consensus(r1_file, r2_file, consensus_length=50, sample_size=5000, threshold=0.8):
    result = {}
    
    print(f'Sampling {sample_size} reads from input files...')
    r1_sequences = sample_reads(r1_file, sample_size)
    r2_sequences = sample_reads(r2_file, sample_size)
    
    # Process R1 file
    print('Calculating R1 consensus sequence...')
    trimmed_r1_sequences = []
    for seq in r1_sequences:
        trimmed_seq = seq[:consensus_length]
        if len(trimmed_seq) == consensus_length:
            trimmed_r1_sequences.append(trimmed_seq)
    
    if trimmed_r1_sequences:
        r1_consensus = ''
        
        for i in range(consensus_length):
            bases = [seq[i] for seq in trimmed_r1_sequences if i < len(seq)]
            base_counts = Counter(bases)
            total = len(bases)
            
            if total == 0:
                r1_consensus += 'N'
                continue
            
            frequencies = {
                'A': base_counts.get('A', 0) / total,
                'C': base_counts.get('C', 0) / total,
                'G': base_counts.get('G', 0) / total,
                'T': base_counts.get('T', 0) / total
            }
            
            most_common = max(frequencies.items(), key=lambda x: x[1])
            base, freq = most_common
            
            if freq >= threshold:
                r1_consensus += base
            else:
                r1_consensus += 'N'
        
        result['R1'] = r1_consensus
        
        if r1_consensus[-4:] != 'NNNN':
            print(f"Warning: The R1 consensus sequence does not end on 'NNNN' (-{r1_consensus[-4:]}). "
                    f"Consider increasing consensus_length beyond {consensus_length} bp.")
    
    # Process R2 file if provided
    print('Calculating R2 consensus sequence...')
    trimmed_r2_sequences = []
    for seq in r2_sequences:
        trimmed_seq = seq[:consensus_length]
        if len(trimmed_seq) == consensus_length:
            trimmed_r2_sequences.append(trimmed_seq)
    
    if trimmed_r2_sequences:
        r2_consensus = ''
        
        for i in range(consensus_length):
            bases = [seq[i] for seq in trimmed_r2_sequences if i < len(seq)]
            base_counts = Counter(bases)
            total = len(bases)
            
            if total == 0:
                r2_consensus += 'N'
                continue
            
            frequencies = {
                'A': base_counts.get('A', 0) / total,
                'C': base_counts.get('C', 0) / total,
                'G': base_counts.get('G', 0) / total,
                'T': base_counts.get('T', 0) / total
            }
            
            most_common = max(frequencies.items(), key=lambda x: x[1])
            base, freq = most_common
            
            if freq >= threshold:
                r2_consensus += base
            else:
                r2_consensus += 'N'
        
        result['R2'] = r2_consensus
        
        if r2_consensus[-4:] != 'NNNN':
            print(f"Warning: The R2 consensus sequence does not end on 'NNNN' (-{r2_consensus[-4:]}). "
                    f"Consider increasing consensus_length beyond {consensus_length} bp.")
    
    return result