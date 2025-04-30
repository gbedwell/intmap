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
import pandas as pd
from scipy import stats
from scipy.interpolate import BSpline
import statsmodels.api as sm
import patsy
import sys
import warnings

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
    'generate_consensus',
    'sonic_abundance'
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

def final_pass_collapse(kept_frags, len_diff, nthr, min_count, count_fc, ltr_cufp, linker_cufp):
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
                    
    results = Parallel(n_jobs=1)(
        delayed(collapse_group)(
            chrom_strand_tuple = key, 
            pos_counts = pos, 
            min_count = min_count, 
            count_fc = count_fc, 
            len_diff = len_diff, 
            read_mapping = read_mapping)
        for key, pos in grouped_data.items()
    )
    
    for read_updates in results:
        for read_name, pos in read_updates.items():
            if kept_frags[read_name]['strand'] == '-':
                kept_frags[read_name]['end'] = pos
            else:
                kept_frags[read_name]['start'] = pos
    
    reads_to_remove = set()
    
    if ltr_cufp:
        ltr_umi_groups = defaultdict(list)
        all_n_ltr_umis = True
        
        for read_name, frag in kept_frags.items():
            if frag['ltr_umi'] != 'N':
                all_n_ltr_umis = False
                break
        
        if not all_n_ltr_umis:
            for read_name, frag in kept_frags.items():
                ltr_umi_groups[frag['ltr_umi']].append(read_name)
            
            # Process each LTR UMI group
            for ltr_umi, read_names in ltr_umi_groups.items():
                if len(read_names) > 1:
                    sorted_reads = sorted(
                        read_names,
                        key=lambda x: (kept_frags[x].get('count', 0), kept_frags[x].get('mean_qual', 0)),
                        reverse=True
                    )
                    
                    reads_to_remove.update(sorted_reads[1:])
    
    if linker_cufp:
        linker_umi_groups = defaultdict(list)
        all_n_linker_umis = True
        
        for read_name, frag in kept_frags.items():
            if frag['linker_umi'] != 'N':
                all_n_linker_umis = False
                break
        
        if not all_n_linker_umis:
            remaining_reads = {read for read in kept_frags if read not in reads_to_remove}
            
            for read_name in remaining_reads:
                frag = kept_frags[read_name]
                linker_umi_groups[frag['linker_umi']].append(read_name)
            
            for linker_umi, read_names in linker_umi_groups.items():
                if len(read_names) > 1:
                    sorted_reads = sorted(
                        read_names,
                        key=lambda x: (kept_frags[x].get('count', 0), kept_frags[x].get('mean_qual', 0)),
                        reverse=True
                    )
                    
                    reads_to_remove.update(sorted_reads[1:])
    
    removed_frags = {}
    for read_name in reads_to_remove:
        removed_frags[read_name] = kept_frags[read_name]
        del kept_frags[read_name]
    
    return kept_frags, removed_frags

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
    report.append(f"Allowed error rate: {error_rate}")
    report.append(f"Reporting threshold: {threshold}")
    report.append(f"Sequences with LTR match: {len(r1_matches)} ({(len(r1_matches) / sample_size) * 100:.2f}%)")
    report.append(f"Sequences with linker match: {len(r2_matches)} ({(len(r2_matches) / sample_size) * 100:.2f}%)")
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

def pad_tab(tab, start_at=None, end_at=None):
    tab_keys = [int(k) for k in tab.keys()]
    
    if start_at is None:
        start_at = 0
    if end_at is None:
        end_at = max(tab_keys) - 1 + 10
    
    new_tab_names = list(range(start_at, end_at + 1))
    
    if not all(k in new_tab_names for k in tab_keys):
        raise ValueError("Invalid lengths - perhaps some are 0")
    
    new_tab = {str(i): 0.0 for i in new_tab_names}
    
    for k, v in tab.items():
        new_tab[str(k)] = v
    
    result = {
        'y': list(new_tab.values()),
        'x': new_tab_names,
        'orig': [str(x) in tab.keys() for x in new_tab_names]
    }
    
    return result

def mstep(x, theta, phi):
    x = np.asarray(x)
    phi = np.asarray(phi)
    theta = np.asarray(theta)
    
    exptp = np.exp(-np.outer(phi, theta))
    
    denom = 1 - exptp
    
    # Calculate gradient
    grad = np.sum(-(1-x) * phi[:, np.newaxis] + x * phi[:, np.newaxis] * exptp / denom, axis = 0)
    
    # Calculate curvature (second derivative)
    curv = -np.sum(phi[:, np.newaxis]**2 * exptp / denom, axis=0)
    
    # Update theta using Newton-Raphson step
    theta = theta - grad / curv
    
    return {"theta": theta, "grad": grad, "info": curv}

def estep(x, theta, phi):
    x = np.asarray(x)
    phi = np.asarray(phi)
    theta = np.asarray(theta)
    
    # Compute lambda as outer product of phi and theta
    lambda_mat = np.outer(phi, theta)
    
    # Compute exp(-lambda)
    expnl = np.exp(-lambda_mat)
    
    # Compute denominator
    denom = 1 - expnl
    
    # Compute result
    res = np.zeros_like(lambda_mat, dtype=float)
    non_zero_denom = (denom != 0)
    res[non_zero_denom] = x[non_zero_denom] * lambda_mat[non_zero_denom] / denom[non_zero_denom]
    
    zero_denom = ~non_zero_denom
    if np.any(zero_denom):
        res[zero_denom] = x[zero_denom]
    
    # Handle case where denominator is zero
    wd = np.where(denom == 0)
    if len(wd[0]) > 0:
        res[wd] = x[wd]
    
    # Compute log-likelihood
    loglik = np.sum(-lambda_mat[x == 0]) + np.sum(np.log1p(-expnl[x == 1]))

    return res, loglik

def phi_update_default(obj):
    class SuppressSpecificWarnings:
        def __enter__(self):
            # Store the original filters
            self.original_filters = warnings.filters.copy()
            
            warnings.filterwarnings("ignore", category=RuntimeWarning, 
                                    message="overflow encountered in exp")
            warnings.filterwarnings("ignore", category=RuntimeWarning, 
                                    message="divide by zero encountered in divide")
            warnings.filterwarnings("ignore", category=RuntimeWarning, 
                                        message="invalid value encountered in multiply")
            warnings.filterwarnings("ignore", category=RuntimeWarning, 
                                        message="overflow encountered in divide")
            
            return self
            
        def __exit__(self, exc_type, exc_val, exc_tb):
            warnings.filters = self.original_filters
            
    row_sums = np.sum(obj, axis = 1)
    row_dict = {str(i): row_sums[i] for i in range(len(row_sums))}
    tmp_frame = pad_tab(row_dict)
    
    df = pd.DataFrame({
        'y': tmp_frame['y'],
        'x': tmp_frame['x']
    })
    
    df['y'] = df['y'] + 1e-5
    
    significant_indices = df[df['y'] > 1e-5]['x'].values
    significant_values = df[df['y'] > 1e-5]['y'].values
    
    knot_configs = []
    
    if len(significant_indices) < 4:
        knot_configs.append({
        'knots': list(np.percentile(significant_indices, [50]).astype(int)),
        'degree': max(1, len(significant_indices) - 1)
        })

    knot_configs.append({
        'knots': [50, 100],
        'degree': 3
    })
    
    knot_configs.append({
        'knots': list(np.percentile(significant_indices, [25, 75]).astype(int)),
        'degree': 3
    })
    
    x_range = df['x'].max() - df['x'].min()
    knot_configs.append({
        'knots': [x_range // 3, 2 * x_range // 3],
        'degree': 3
    })
    
    with SuppressSpecificWarnings():
        for config in knot_configs:
            try:
                knot_formula = f"bs(x, knots={config['knots']}, degree={config['degree']}, include_intercept=True)"
                X = patsy.dmatrix(knot_formula, df)
                
                # Fit initial Poisson GLM
                family = sm.families.Poisson()
                model = sm.GLM(df['y'], X, family=family)
                initial_fit = model.fit(tol=1e-8)
                
                # Calculate dispersion parameter (scale) using Pearson chi-square
                mu = np.maximum(initial_fit.mu, 1e-5)
                pearson_chi2 = ((df['y'] - mu) ** 2 / mu).sum()
                scale = max(1.0, pearson_chi2 / initial_fit.df_resid)
                
                # Refit the model with the estimated scale parameter
                fit = model.fit(scale=scale, tol=1e-8)
                preds = np.maximum(np.exp(fit.predict()), 1e-5)

                break
                
            except Exception as e:
                if config == knot_configs[-1]:
                    # print(f"All GLM fitting attempts failed. Last error: {str(e)}")
                    preds = df['y'].values  # Just use the input values as a last resort
                else:
                    continue
                    
    orig_indices = np.where(tmp_frame['orig'])[0]
    orig_preds = preds[orig_indices]
    return orig_preds / np.sum(orig_preds)

def maxEM(slmat):
    # Initialize phi and theta
    phi_old = np.sum(slmat + 1/slmat.size, axis = 1)
    phi_old = phi_old / np.sum(phi_old)
    theta_old = np.sum(slmat, axis = 0)
    
    # Do one EM step to start
    # phi_update assures that constraints on phi are met
    Y, _ = estep(slmat, theta_old, phi_old)
    theta_old = np.sum(Y, axis = 0)
    
    row_sums = np.sum(Y, axis = 1)
    phi_old = phi_update_default(Y)
    
    phi_min = np.finfo(float).eps
    phi_old = phi_old + phi_min
    phi_old = phi_old / np.sum(phi_old)
    
    llk_old = float('-inf')
    i = 0
    not_done = True
    
    min_reps = 3
    max_reps = 2000
    max_abs_le = 0.1
    max_rel_le = 1e-5
    
    phi_new = phi_old.copy()
    theta_new = theta_old.copy()
    
    while not_done:
        # M-step
        mres = mstep(slmat, theta_new, phi_new)
        theta_new = mres["theta"]
        
        # E-step
        Y, llk = estep(slmat, theta_new, phi_new)
        
        # Check increase in loglik - use one full EM step if llk doesn't increase
        fail = 0
        if llk > llk_old:
            llk_old = llk
            theta_old = theta_new.copy()
        else:
            fail += 1
            theta_new = theta_old.copy()
            Y, _ = estep(slmat, theta_new, phi_new)
        
        # Update phi
        phi_new = phi_update_default(Y) + phi_min
        phi_new = phi_new / np.sum(phi_new)
        
        # Add small value to avoid numerical issues
        phi_new = phi_new + phi_min
        phi_new = phi_new / np.sum(phi_new)
        
        # Check if phi update improved likelihood
        Y, llk = estep(slmat, theta_new, phi_new)
        if llk > llk_old:
            llk_old = llk
            phi_old = phi_new.copy()
        else:
            phi_new = phi_old.copy()
            fail += 1
        
        # Check convergence criteria
        adjs = mres["grad"]
        max_abs = np.max(np.abs(adjs))
        max_rel = np.max(np.abs(adjs / (theta_old + 1e-10)))
        
        i += 1
        not_done = (fail < 2) and ((i < min_reps) or (max_abs > max_abs_le) or (max_rel > max_rel_le))
        
        if i >= max_reps:
            print("Warning: iteration limit reached")
            break
    
    return {
        "theta": theta_new,
        "phi": phi_new,
        "iter": i
    }
    
def estimate_abundance(location_indices, lengths, min_length=25):
    location_indices = np.array(location_indices)
    lengths = np.array(lengths)
    
    if np.min(lengths) < min_length:
        raise ValueError(f"Minimum length ({np.min(lengths)}) is less than required ({min_length})")
    
    if len(location_indices) != len(lengths):
        raise ValueError("Lengths of location_indices and lengths must be equal")
    
    # Find the actual min and max lengths in the data
    min_length_actual = int(np.min(lengths))
    max_length_actual = int(np.max(lengths))
    
    # Create a cross-tabulation matrix
    n_locations = int(np.max(location_indices)) + 1
    length_range = np.arange(min_length_actual, max_length_actual + 1)
    
    # Initialize the matrix
    slmat = np.zeros((len(length_range), n_locations), dtype=int)
    
    # Fill the matrix with counts (not just binary indicators)
    for loc_idx, length in zip(location_indices, lengths):
        length_idx = int(length) - min_length_actual
        slmat[length_idx, int(loc_idx)] += 1
    
    # Call maxEM with the matrix
    result = maxEM(slmat)
    
    # Add additional fields to match R implementation
    result['obs'] = np.bincount(location_indices.astype(int), minlength=n_locations)
    
    # Store the original data
    result['data'] = {
        'locations': location_indices,
        'lengths': lengths
    }
    
    return result

def sonic_abundance(frag_dict, U3, shift, min_frag_len):
    locations = []
    lengths = []
    
    # Extract locations and lengths from fragment dictionary
    for key, value in frag_dict.items():
        if value['strand'] == '-':
            site_start = (value['end'] - shift) - 1
            site_end = (value['end'] - shift)
            site_strand = '-' if not U3 else '+'
        else:
            site_start = (value['start'] + shift)
            site_end = (value['start'] + shift) + 1
            site_strand = '+' if not U3 else '-'

        frag_length = value['end'] - value['start']
        
        site_location = f"{value['chrom']}:{site_start if site_strand == '+' else site_end}:{site_strand}"
        
        locations.append(site_location)
        lengths.append(frag_length)
    
    # Check if we have enough data
    if len(locations) > 0:
        try:
            # Convert locations to numeric indices
            unique_locations = np.unique(locations)

            # Check if each location has exactly one fragment
            # if len(unique_locations) == len(locations):
            #     print('Each integration site has exactly one fragment. Skipping MLE.', flush = True)
            #     return []
            
            location_indices = np.array([np.where(unique_locations == loc)[0][0] for loc in locations])
            
            # Estimate abundance
            abundance_result = estimate_abundance(location_indices, lengths, min_length = min_frag_len)

            # Format results
            abundance_estimates = []
            for i, loc in enumerate(unique_locations):
                chrom, position, strand = loc.split(':')
                position = int(position)
                
                if strand == '+':
                    start = position
                    end = position + 1
                else:
                    end = position
                    start = position - 1
                
                abundance_estimates.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': '.',
                    'count': abundance_result['theta'][i],
                    'strand': strand
                })
            
                abundance_estimates.sort(key=lambda item: (natural_key(item['chrom']), item['start'], item['end']))

            return abundance_estimates
            
        except Exception as e:
            print(f'Warning: Could not estimate abundance: {str(e)}', flush=True)
            import traceback
            traceback.print_exc()
            return []
    else:
        print('Not enough data for abundance estimation', flush=True)
        return []