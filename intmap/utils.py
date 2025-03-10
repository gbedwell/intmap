import os
import gzip
import io
import regex
import numpy as np
import faiss
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict
from datasketch import MinHash
import math
import multiprocessing
from joblib import Parallel, delayed
from operator import itemgetter

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