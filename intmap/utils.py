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
import hashlib

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
    
    unique_entries = []
    grouped_entries = []
    
    label_counts = np.bincount(labels)
    
    for component_id in range(n_components):
        if label_counts[component_id] == 1:
            unique_entries.append(int(np.where(labels == component_id)[0][0]))
        elif label_counts[component_id] > 1:
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
    positions = np.array([p for p, _, _ in pos_counts])
    counts = np.array([c for _, c, _ in pos_counts])
    dup_counts = np.array([dc for _, _, dc in pos_counts])
    
    abundant_mask = counts >= min_count
    abundant_positions = positions[abundant_mask]
    abundant_counts = counts[abundant_mask]
    abundant_dup_counts = dup_counts[abundant_mask]
    
    sort_idx = np.argsort(abundant_counts)[::-1]
    abundant_positions = abundant_positions[sort_idx]
    abundant_counts = abundant_counts[sort_idx]
    abundant_dup_counts = abundant_dup_counts[sort_idx]

    keep_abundant = np.ones(len(abundant_positions), dtype=bool)
    collapsed = set()
    for i in range(1, len(abundant_positions)):
        if i in collapsed:
            continue
            
        valid_positions_mask = np.array([j not in collapsed for j in range(i)])
        valid_positions = abundant_positions[:i][valid_positions_mask]
        
        if len(valid_positions) > 0:
            nearby_mask = np.abs(valid_positions - abundant_positions[i]) <= len_diff
            if any(nearby_mask):
                valid_counts = abundant_counts[:i][valid_positions_mask]
                valid_dup_counts = abundant_dup_counts[:i][valid_positions_mask]
                
                max_nearby_count = np.max(valid_counts[nearby_mask])
                max_nearby_dup_count = np.max(valid_dup_counts[nearby_mask])
                
                if max_nearby_count >= abundant_counts[i] * count_fc:
                    keep_abundant[i] = False
                    collapsed.add(i)
                elif max_nearby_dup_count >= abundant_dup_counts[i] * count_fc:
                    keep_abundant[i] = False
                    collapsed.add(i)
    
    final_abundant_positions = abundant_positions[keep_abundant]
    
    read_updates = {}
    for abundant_pos in final_abundant_positions:
        collapse_mask = np.abs(positions - abundant_pos) <= len_diff
        for pos in positions[collapse_mask]:
            for read_name in read_mapping[(chrom_strand_tuple[0], pos, chrom_strand_tuple[1])]:
                read_updates[read_name] = abundant_pos
    
    return read_updates

def consolidate_new_dups(kept_dict, include_umis = False):
      fragment_data = []
      read_names = []

      for read_name, frag in kept_dict.items():
          if include_umis:
              fragment_data.append((frag['chrom'], frag['start'], frag['end'], frag['strand'], frag['ltr_umi'], frag['linker_umi']))
          else:
              fragment_data.append((frag['chrom'], frag['start'], frag['end'], frag['strand']))
              
          read_names.append(read_name)
      
      if include_umis:
          fragment_data = np.array(fragment_data, dtype=[('chrom', 'U25'), ('start', 'i8'), ('end', 'i8'), ('strand', 'U1'), ('ltr_umi', 'U50'), ('linker_umi', 'U50')])
      else:
          fragment_data = np.array(fragment_data, dtype=[('chrom', 'U25'), ('start', 'i8'), ('end', 'i8'), ('strand', 'U1')])
      
      read_names = np.array(read_names)

      unique_frags, inverse_indices, counts = np.unique(fragment_data, return_inverse = True, return_counts = True)

      duplicate_mask = counts > 1
      duplicate_indices = np.where(duplicate_mask)[0]

      reads_to_remove = set()

      for group_idx in duplicate_indices:
          group_mask = inverse_indices == group_idx
          group_reads = read_names[group_mask]
          
          best_read = max(group_reads, 
                        key=lambda x: (
                          kept_dict[x].get('count', 0), 
                          kept_dict[x].get('map_qual', 0),
                          kept_dict[x].get('mean_qual', 0)
                          ))
          
          total_count = sum(kept_dict[read].get('count', 0) for read in group_reads)
          kept_dict[best_read]['count'] = total_count
          
          for read in group_reads:
              if read != best_read:
                  reads_to_remove.add(read)
                  
      return reads_to_remove
  
def final_pass_collapse(kept_frags, len_diff, nthr, min_count, count_fc, ltr_cufp, linker_cufp,
                        ltr_umi_len, linker_umi_len):
    positions = np.array([(frag['chrom'], 
                            frag['end'] if frag['strand'] == '-' else frag['start'],
                            frag['strand'])
                            for frag in kept_frags.values()],
                            dtype=[('chrom', 'U25'), ('pos', 'i8'), ('strand', 'U1')])
    
    frag_counts = np.array([frag['count'] for frag in kept_frags.values()])
    unique_pos, inverse_indices, counts = np.unique(positions, return_counts = True, return_inverse = True)

    dup_counts = np.bincount(inverse_indices, weights = frag_counts)
    
    read_mapping = defaultdict(list)
    for read_name, frag in kept_frags.items():
        pos = frag['end'] if frag['strand'] == '-' else frag['start']
        read_mapping[(frag['chrom'], pos, frag['strand'])].append(read_name) 
    
    grouped_data = defaultdict(list)
    for pos, count, dup_count in zip(unique_pos, counts, dup_counts):
        grouped_data[(pos['chrom'], pos['strand'])].append((pos['pos'], count, dup_count))
                    
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
    
    if ltr_umi_len == 0 and linker_umi_len == 0:
        reads_to_remove = consolidate_new_dups(kept_frags, include_umis = False)
    else:
        reads_to_remove = consolidate_new_dups(kept_frags, include_umis = True)
    
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
            
            for ltr_umi, read_names in ltr_umi_groups.items():
                if len(read_names) > 1:
                    sorted_reads = sorted(
                        read_names,
                        key=lambda x: (
                          kept_frags[x].get('count', 0), 
                          kept_frags[x].get('map_qual', 0),
                          kept_frags[x].get('mean_qual', 0)),
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
                        key=lambda x: (
                          kept_frags[x].get('count', 0), 
                          kept_frags[x].get('map_qual', 0),
                          kept_frags[x].get('mean_qual', 0)),
                        reverse=True
                    )
                    
                    reads_to_remove.update(sorted_reads[1:])
    
    if len(reads_to_remove) > 0:
        print(f"Removing {len(reads_to_remove)} additional reads in the final pass.", flush=True)
        
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
            if line_num % 4 == 1:
                if read_count < n_reads:
                    sequences.append(line.strip())
                else:
                    j = random.randint(0, read_count)
                    if j < n_reads:
                        sequences[j] = line.strip()
                
                read_count += 1
            
            line_num += 1
    
    return sequences

def filter_contaminants(seqs, contams, contam_error_rate=0.3, min_seqs_threshold=0.1):
    if not contams:
        return seqs, {'original_count': len(seqs), 'filtered_count': 0, 'remaining_count': len(seqs)}
    
    original_count = len(seqs)
    
    contam_patterns = []
    for contam in contams:
        if len(contam) < 10:
            raise ValueError('Contaminants must be >= 10 nucleotides long.')
        
        contam_check = contam.replace(' ', '')
        contam_char = regex.sub('[ATGC]', '', contam_check.upper())
        if contam_char != '':
            raise ValueError('Contaminants may only contain A,T,G, and C nucleotides.')
        
        max_errors = math.floor(len(contam) * contam_error_rate)
        pattern = regex.compile(f'({contam}){{e<={max_errors}}}', regex.IGNORECASE)
        contam_patterns.append(pattern)
    
    kept_seqs = []
    contam_count = 0
    
    for seq in seqs:
        is_contaminated = False
        for pattern in contam_patterns:
            if pattern.search(seq):
                is_contaminated = True
                break
        
        if not is_contaminated:
            kept_seqs.append(seq)
        else:
            contam_count += 1
    
    remaining_count = len(kept_seqs)
    filtered_count = original_count - remaining_count
    
    remaining_frac = remaining_count / original_count
    if remaining_frac < min_seqs_threshold:
        print(f"WARNING: Contaminant filtering removed {filtered_count} seqs "
              f"({(1-remaining_frac)*100:.1f}%). "
              f"Consider adjusting contaminant seqs or error rates.", flush=True)
    
    filtering_stats = {
        'original_count': original_count,
        'filtered_count': filtered_count,
        'remaining_count': remaining_count,
        'remaining_fraction': remaining_frac
    }
    
    return kept_seqs, filtering_stats


def check_consensus(r1_file, r2_file, ltr_seq, linker_seq, sample_size=5000, threshold=0.02, 
                    error_rate=0.3, contams=None, contam_error_rate=0.3, 
                    min_seqs_threshold=0.1):

    print(f"Sampling {sample_size} reads from input files...", flush=True)
    r1_sequences = sample_reads(r1_file, sample_size)
    r2_sequences = sample_reads(r2_file, sample_size)
    
    if contams:
        r1_filter_stats = {'original_count': len(r1_sequences), 'filtered_count': 0, 'remaining_count': len(r1_sequences)}
        # r2_filter_stats = {'original_count': len(r2_sequences), 'filtered_count': 0, 'remaining_count': len(r2_sequences)}
        
        print(f"Filtering contaminants from sampled R1 reads...", flush=True)
        
        r1_sequences, r1_filter_stats = filter_contaminants(
            r1_sequences, contams, contam_error_rate, min_seqs_threshold
        )
        # r2_sequences, r2_filter_stats = filter_contaminants(
        #     r2_sequences, contams, contam_error_rate, min_seqs_threshold
        # )
        
        # print(f"R1: Filtered {r1_filter_stats['filtered_count']} sequences, "
        #       f"{r1_filter_stats['remaining_count']} remaining ({r1_filter_stats['remaining_fraction']*100:.1f}%)", flush=True)
        # print(f"R2: Filtered {r2_filter_stats['filtered_count']} sequences, "
        #       f"{r2_filter_stats['remaining_count']} remaining ({r2_filter_stats['remaining_fraction']*100:.1f}%)", flush=True)
    
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
    
    r1_matches = []
    r2_matches = []
    
    print("Searching for LTR/linker sequences in filtered reads...", flush=True)
    
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
    
    print(f'Found {len(r1_matches)} R1 sequences and {len(r2_matches)} R2 sequences with matches', flush=True)
    
    if not r1_matches or not r2_matches:
        if not r1_matches and not r2_matches:
            print(f'No matching LTR or linker sequences found in the filtered reads. Check the given inputs.')
        if not r1_matches:
            print(f'No matching LTR sequences found in the filtered reads. Check the given input.')
        if not r2_matches:
            print(f'No matching linker sequences found in the filtered reads. Check the given input.')
        sys.exit()
    
    r1_counts = {}
    for match in r1_matches:
        seq = match['matched_text']
        r1_counts[seq] = r1_counts.get(seq, 0) + 1
    
    r2_counts = {}
    for match in r2_matches:
        seq = match['matched_text']
        r2_counts[seq] = r2_counts.get(seq, 0) + 1
    
    r1_total = len(r1_matches)
    r2_total = len(r2_matches)
    
    r1_freqs = {seq: count/r1_total for seq, count in r1_counts.items()}
    r2_freqs = {seq: count/r2_total for seq, count in r2_counts.items()}
    
    report = []
    report.append("")
    
    report.append(f"Target LTR sequence: {ltr_seq}")
    report.append(f"Target Linker sequence: {linker_seq}")
    report.append(f"Sample size: {sample_size} reads")
    report.append(f"Allowed error rate: {error_rate}")
    report.append(f"Reporting threshold: {threshold}")
    
    if contams:
        report.append(f"R1 sequences remaining after contaminant filtering: {r1_filter_stats['remaining_count']} "
                      f"({r1_filter_stats['remaining_fraction']*100:.1f}%)")
        # report.append(f"R2 sequences remaining after contaminant filtering: {r2_filter_stats['remaining_count']} "
        #               f"({r2_filter_stats['remaining_fraction']*100:.1f}%)")
    
    effective_r1_sample = len(r1_sequences)
    effective_r2_sample = len(r2_sequences)
    
    report.append(f"Sequences with LTR match: {len(r1_matches)} ({(len(r1_matches) / effective_r1_sample) * 100:.2f}%)")
    report.append(f"Sequences with linker match: {len(r2_matches)} ({(len(r2_matches) / effective_r2_sample) * 100:.2f}%)")
    report.append("")
    
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


def generate_consensus(r1_file, r2_file, consensus_length=50, sample_size=5000, threshold=0.8,
                      contams=None, contam_error_rate=0.3, min_seqs_threshold=0.1):
    result = {}
    
    print(f'Sampling {sample_size} reads from input files...')
    r1_sequences = sample_reads(r1_file, sample_size)
    r2_sequences = sample_reads(r2_file, sample_size)
    
    if contams:
        print(f"Filtering contaminants from sampled reads...", flush=True)
        
        r1_sequences, r1_filter_stats = filter_contaminants(
            r1_sequences, contams, contam_error_rate, min_seqs_threshold
        )
        # r2_sequences, r2_filter_stats = filter_contaminants(
        #     r2_sequences, contaminants, contam_error_rate, min_sequences_threshold
        # )
        
        print(f"Removed {r1_filter_stats['filtered_count']} R1 sequences, "
              f"{r1_filter_stats['remaining_count']} sequences remaining ({r1_filter_stats['remaining_fraction']*100:.1f}%)", flush=True)
        # print(f"Removed {r2_filter_stats['filtered_count']} R2 sequences, "
        #       f"{r2_filter_stats['remaining_count']} sequences remaining ({r2_filter_stats['remaining_fraction']*100:.1f}%)", flush=True)
    
    print('Calculating LTR consensus sequence...')
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
    
    print('Calculating linker consensus sequence...')
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
  
def pad_tab(tab, start_at = None, end_at = None):
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
    
    log_exptp = -np.outer(phi, theta)
    
    log_denom = np.log1p(-np.exp(log_exptp))
    
    grad = np.sum(-(1 - x) * phi[:, np.newaxis] + 
                  x * phi[:, np.newaxis] * np.exp(log_exptp - log_denom), axis = 0)
    
    curv = -np.sum(phi[:, np.newaxis]**2 * np.exp(log_exptp - log_denom), axis = 0)
    
    theta = theta - grad / curv
    
    return {"theta": theta, "grad": grad, "info": curv}

def estep(x, theta, phi):
    x = np.asarray(x)
    phi = np.asarray(phi)
    theta = np.asarray(theta)
    
    lambda_mat = np.outer(phi, theta)
    
    log_expnl = -lambda_mat
    expnl = np.exp(np.clip(log_expnl, -500, 500))
    expnl[log_expnl < -500] = 0.0
    
    denom = 1 - expnl
    
    safe_denom = np.maximum(denom, np.finfo(float).eps)
    res = x * lambda_mat / safe_denom
    
    first_term = np.sum(-lambda_mat[x == 0])
    # Account for cases expnl is close to 1
    safe_expnl = np.minimum(expnl, 1.0 - np.finfo(float).eps)
    second_term = np.sum(np.log1p(-safe_expnl[x == 1]))
        
    loglik = first_term + second_term

    return res, loglik

def safer_normalization(preds, mask):
    indices = np.where(mask)[0]
    selected_preds = preds[indices]
    
    if np.any(selected_preds <= np.finfo(float).eps):
        if np.sum(selected_preds) == 0:
            return np.zeros_like(selected_preds)
        return selected_preds / np.sum(selected_preds)
    
    log_preds = np.log(selected_preds)
    max_log = np.max(log_preds)
    exp_shifted = np.exp(log_preds - max_log)
    return exp_shifted / np.sum(exp_shifted)
  
def phi_update_default(obj):
    row_sums = np.sum(obj, axis = 1)
    row_dict = {str(i): row_sums[i] for i in range(len(row_sums))}
    tmp_frame = pad_tab(row_dict)
    
    df = pd.DataFrame({
        'y': tmp_frame['y'],
        'x': tmp_frame['x']
    })
    
    df['y'] = df['y'] + 1e-10
    
    non_zero_indices = df[df['y'] > 1e-5]['x'].values
    
    knot_configs = []
    
    if len(non_zero_indices) <= 3:
        knot_configs.append({
            'knots': [],
            'degree': min(max(1, len(non_zero_indices) - 1), 2)
        })
    elif len(non_zero_indices) <= 10:
        median_idx = int(np.median(non_zero_indices))
        knot_configs.append({
            'knots': [median_idx],
            'degree': 2
        })
    else:
        quantiles = [25, 75] if len(non_zero_indices) < 20 else [20, 50, 80]
        knots = list(np.percentile(non_zero_indices, quantiles).astype(int))
        knot_configs.append({
            'knots': knots,
            'degree': 3
        })
    
    x_range = df['x'].max() - df['x'].min()
    if x_range > 10:
        knot_configs.append({
            'knots': [df['x'].min() + x_range // 3, df['x'].min() + 2 * x_range // 3],
            'degree': 3
        })
    
    unsuccessful_update = False
    for config in knot_configs:
        try:
            knot_formula = f"bs(x, knots={config['knots']}, degree={config['degree']}, include_intercept=True)"
            X = patsy.dmatrix(knot_formula, df)
            
            family = sm.families.Poisson()
            model = sm.GLM(df['y'], X, family = family)
            initial_fit = model.fit(tol = 1e-12, max_iter = 1000)
            
            # Calculate scale using Pearson residuals and residual DF
            mu = np.maximum(initial_fit.mu, 1e-5)
            pearson_resid = initial_fit.resid_pearson
            df_resid = initial_fit.df_resid
            scale = sum(pearson_resid**2) / df_resid
            
            fit = model.fit(scale = scale, tol = 1e-12, max_iter = 1000)
            predicted_values = fit.predict()
            preds = np.maximum(predicted_values, 1e-5)
            break
            
        except Exception as e:
            if config == knot_configs[-1]:
                preds = df['y'].values
                unsuccessful_update = True 
            else:
                continue

    return safer_normalization(preds, tmp_frame['orig']), unsuccessful_update

def maxEM(slmat):
    phi_old = np.sum(slmat, axis = 1)
    phi_old = np.maximum(phi_old, np.finfo(float).eps)
    phi_old = phi_old / np.sum(phi_old)
    
    theta_old = np.sum(slmat, axis = 0)
    theta_old = np.maximum(theta_old, np.finfo(float).eps)

    Y, _ = estep(slmat, theta_old, phi_old)
    theta_old = np.sum(Y, axis = 0)
    
    row_sums = np.sum(Y, axis = 1)
    phi_old, _ = phi_update_default(Y)
    
    phi_min = np.finfo(float).eps
    phi_old = np.maximum(phi_old, phi_min)
    phi_old = phi_old / np.sum(phi_old)
    
    llk_old = float('-inf')
    i = 0
    not_done = True
    
    min_reps = 3
    max_reps = 2000
    max_abs_le = 0.01
    max_rel_le = 1e-6
    
    phi_new = phi_old.copy()
    theta_new = theta_old.copy()
    
    converged = True
    fallback = False
    while not_done:
        mres = mstep(slmat, theta_new, phi_new)
        theta_new = np.maximum(mres["theta"], np.finfo(float).eps)
        
        try:
            Y, llk = estep(slmat, theta_new, phi_new)
            
            if not np.isfinite(llk):
                fail = 1
                theta_new = theta_old.copy()
                Y, llk = estep(slmat, theta_new, phi_new)
            else:
                fail = 0
                if llk > llk_old:
                    llk_old = llk
                    theta_old = theta_new.copy()
                else:
                    fail += 1
                    theta_new = theta_old.copy()
                    Y, _ = estep(slmat, theta_new, phi_new)
        except Exception:
            fail = 1
            theta_new = theta_old.copy()
            Y, llk = estep(slmat, theta_new, phi_new)
        
        phi_new, fallback = phi_update_default(Y)
        phi_new = np.maximum(phi_new, phi_min)
        phi_new = phi_new / np.sum(phi_new)

        try:
            Y, llk = estep(slmat, theta_new, phi_new)
            if np.isfinite(llk) and llk > llk_old:
                llk_old = llk
                phi_old = phi_new.copy()
            else:
                phi_new = phi_old.copy()
                fail += 1
        except Exception:
            phi_new = phi_old.copy()
            fail += 1
        
        adjs = mres["grad"]
        max_abs = np.max(np.abs(adjs))
        denominator = np.maximum(np.abs(theta_old), 1e-10)
        max_rel = np.max(np.abs(adjs / denominator))
        
        i += 1
        not_done = (fail < 2) and ((i < min_reps) or (max_abs > max_abs_le) or (max_rel > max_rel_le))
        
        if i >= max_reps:
            converged = False
            break
          
        if fail >= 2:
            converged = False

    if fallback:
        converged = False  
        
    return {
        'theta': theta_new,
        'phi': phi_new,
        'iter': i,
        'converged': converged
        }

def estimate_abundance(location_indices, lengths, min_length):
    location_indices = np.array(location_indices)
    lengths = np.array(lengths)
    
    if len(location_indices) != len(lengths):
        raise ValueError("Lengths of location_indices and lengths must be equal")
    
    min_length_actual = int(np.min(lengths))
    max_length_actual = int(np.max(lengths))
    
    n_locations = int(np.max(location_indices)) + 1
    length_range = np.arange(min_length_actual, max_length_actual + 1)
    
    slmat = np.zeros((len(length_range), n_locations), dtype = int)
    length_indices = (lengths - min_length_actual).astype(int)
    loc_indices = location_indices.astype(int)
    np.add.at(slmat, (length_indices, loc_indices), 1)
    
    row_sums = slmat.sum(axis=1)
    keep_rows = row_sums > 0

    slmat = slmat[keep_rows, :]

    result = maxEM(slmat)
    
    if result['converged']:
        print(f"Model converged after {result['iter']} iterations")
    else:
        print(f'Model failed to converge; abundance estimates could be unreliable')
    
    result['obs'] = np.bincount(location_indices.astype(int), minlength = n_locations)
    
    result['data'] = {
          'locations': location_indices,
          'lengths': lengths
    }
    
    return result

def sonic_abundance(frag_dict, U3, shift, min_frag_len):
    locations = []
    lengths = []
    
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
        site_count = value['count']
        
        locations.append(site_location)
        lengths.append(frag_length)
    
    if len(locations) > 0:
        try:
            unique_locations = np.unique(locations)
            location_indices = np.array([np.where(unique_locations == loc)[0][0] for loc in locations])
            
            abundance_result = estimate_abundance(location_indices, lengths, min_length = min_frag_len)

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