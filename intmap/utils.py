import os
import gzip
import io
import regex
import numpy as np
import faiss
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict
from joblib import Parallel, delayed
from datasketch import MinHash
import math
from bisect import bisect_left, bisect_right

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
    'chunk_hashes',
    'group_similar_hashes',
    'extract_grouped_entries',
    'TRANS_TABLE',
    'revcomp',
    'collapse_group',
    'final_pass_collapse'
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
    return [int(text) if text.isdigit() else text for text in regex.split(r'(\d+)', chrom)]

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
    
    if strand == '-':
        start = (coordinates['end'] - shift) - 40
        end = (coordinates['end'] - shift) + 15
    else:
        start = (coordinates['start'] + shift) - 15
        end = (coordinates['start'] + shift) + 40
        
    sequence = genome.fetch(chrom, start, end)
    
    if strand == '-':
        sequence = crop.revcomp(sequence)
    return sequence.upper()

def apply_minhash(key, value, string, num_bits, token_size):
    minhash = MinHash(num_perm = num_bits)
    tokens = [string[i:i + token_size] for i in range(len(string) - token_size + 1)]
    
    for token in tokens:
        minhash.update(token.encode('utf8'))
    
    minhash_values = np.array(minhash.hashvalues, dtype = np.uint64)
    
    packed_bits = minhash_values.view(np.uint8)
    
    value['hash'] = packed_bits
    
    return key, value

def set_faiss_threads(nthr):
    if hasattr(faiss, 'omp_set_num_threads'):
        faiss.omp_set_num_threads(nthr)
        return True
    return False

def get_nn(input_dict, num_bits, nthr, len_diff, k):
    set_faiss_threads(nthr)
    
    kmer_groups = defaultdict(dict)
    n_kmers = len_diff + 1
    
    for key, value in input_dict.items():
        seq = value['seq1']
        kmers = sorted([seq[i:i+k] for i in range(len(seq) - k + 1)])
        selected_kmers = kmers[:n_kmers]
        for kmer in selected_kmers:
            kmer_groups[kmer].update({key: value})
    
    all_distances = []
    all_indices = []
    offset = 0
    max_k = max(len(group) for group in kmer_groups.values())
    
    for group in kmer_groups.values():
        minhash_list = [value['hash'] for value in group.values()]
        nb = len(minhash_list)
        if nb == 0:
            continue
            
        dim = num_bits * 64
        db = np.array(minhash_list, dtype='uint8')
        
        hash_index = faiss.IndexBinaryFlat(dim)
        hash_index.add(db)
        
        distances, indices = hash_index.search(db, k = max_k)
        
        indices[indices != -1] += offset
        offset += nb
        
        all_distances.append(distances)
        all_indices.append(indices)
    
    return (np.vstack(all_distances) if all_distances else np.array([]), 
            np.vstack(all_indices) if all_indices else np.array([]))

def chunk_hashes(distances, indices, nthr):
    n = math.ceil(len(distances) / nthr)
    for i in range(0, len(distances), n):
        yield distances[i:i + n], indices[i:i + n]

def group_similar_hashes(hash_distances, hash_indices, nthr, sensitivity, num_bits, input_dict):
    hamming_threshold = (num_bits * 64) * (1 - sensitivity)
    
    # Use input_dict keys to get correct sequence count
    sequence_indices = range(len(input_dict))
    n_sequences = len(sequence_indices)
    
    # Build edges using original sequence indices
    rows, cols = [], []
    mask = hash_distances <= hamming_threshold
    for dist_row, idx_row, row_mask in zip(hash_distances, hash_indices, mask):
        valid_indices = idx_row[row_mask]
        if len(valid_indices) > 0 and idx_row[0] < n_sequences:
            source_idx = idx_row[0]
            valid_indices = valid_indices[valid_indices < n_sequences]
            rows.extend([source_idx] * len(valid_indices))
            cols.extend(valid_indices)
    
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
    
    for component_id in range(n_components):
        component_indices = np.where(labels == component_id)[0]
        if len(component_indices) == 1:
            unique_entries.append(int(component_indices[0]))
        else:
            grouped_entries.append(sorted(map(int, component_indices)))
            
    return unique_entries, grouped_entries

def extract_grouped_entries(groups, input_dict):
    grouped_entries = []
    for group in groups:
        entries = [list(input_dict.values())[index] for index in group]
        grouped_entries.append(entries)

    return grouped_entries

def collapse_group(chrom_strand_tuple, sites, threshold, len_diff):
    sites.sort()
    positions = [s[0] for s in sites]
    position_updates = {}
    
    abundant_sites = [s for s in sites if s[1] >= threshold]
    
    for abundant_pos, abundant_count, abundant_read in abundant_sites:
        left = bisect_left(positions, abundant_pos - len_diff)
        right = bisect_right(positions, abundant_pos + len_diff)
        
        for pos, count, read_name in sites[left:right]:
            if read_name != abundant_read and count < threshold:
                position_updates[read_name] = abundant_pos
    
    return position_updates

def final_pass_collapse(kept_frags, len_diff, nthr, threshold):
    grouped_frags = defaultdict(list)
    for read_name, frag in kept_frags.items():
        key = (frag['chrom'], frag['strand'])
        site_pos = frag['end'] if frag['strand'] == '-' else frag['start']
        grouped_frags[key].append([site_pos, frag['count'], read_name])
    
    results = Parallel(n_jobs=nthr)(
        delayed(collapse_group)(
            key, sites,
            threshold=threshold,
            len_diff=len_diff
        ) for key, sites in grouped_frags.items()
    )
    
    for position_updates in results:
        for read_name, new_pos in position_updates.items():
            if kept_frags[read_name]['strand'] == '-':
                kept_frags[read_name]['end'] = new_pos
            else:
                kept_frags[read_name]['start'] = new_pos
    
    return kept_frags