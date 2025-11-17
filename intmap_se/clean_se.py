import math
from rapidfuzz.distance import Levenshtein
from collections import defaultdict, deque, Counter
import numpy as np
import random
from itertools import groupby
from operator import itemgetter
import random
import hashlib
import faiss
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from intmap.utils import *
import multiprocessing
from joblib import Parallel, delayed
import time
from datasketch import MinHash, LeanMinHash, MinHashLSH
from functools import lru_cache
from pybloom_live import BloomFilter

# def hamming_distance(seq1, seq2):
#     if len(seq1) != len(seq2):
#         return ValueError('Sequences must be of equal length.')
#     return sum(x1 != x2 for x1, x2 in zip(seq1, seq2))

@lru_cache(maxsize=10000)
def seq_similarity(seq1, seq2):
    return 1 - Levenshtein.normalized_distance(seq1, seq2)

def unique_exact_matches(input_dict):
    tmp_exact_kept = defaultdict(list)
    tmp_exact_dup = set()
    
    sort_key = itemgetter('start', 'end', 'ltr_umi', 'linker_umi')
    qual_key = itemgetter('map_qual', 'mean_qual')
    
    for key, entries in input_dict.items():
        entries.sort(key = sort_key)
        
        for _, group_entries in groupby(entries, key = sort_key):
            group = list(group_entries)
            
            if not tmp_exact_dup.intersection(entry['read_name'] for entry in group):
                if len(group) > 1:
                    best_entry = max(group, key=qual_key)
                    best_entry['count'] += sum(
                        entry['count'] for entry in group 
                        if entry is not best_entry
                    )
                    tmp_exact_kept[key].append(best_entry)
                    
                    tmp_exact_dup.update(
                        entry['read_name'] for entry in group 
                        if entry is not best_entry
                    )
                else:
                    tmp_exact_kept[key].extend(group)
    
    return dict(tmp_exact_kept)

def ranged_groupby(entries, tolerance, sort_key = itemgetter('start', 'end')):
    entries = sorted(entries, key=sort_key)
    groups = []
    current_group = []
    last_start, last_end = None, None

    for entry in (e for e in entries):
        start, end = sort_key(entry)
        if current_group and any(
            abs(x) > tolerance for x in (start - last_start, end - last_end)
            ):
            groups.append(current_group)
            current_group = []
        current_group.append(entry)
        last_start, last_end = start, end

    if current_group:
        groups.append(current_group)

    return groups

# def vectorized_hamming_distance(seq1, seq2):
#     if len(seq1) != len(seq2):
#         return ValueError('Sequences must be of equal length.')
#     # Convert strings to arrays of integers (ASCII values)
#     arr1 = np.array([ord(c) for c in seq1])
#     arr2 = np.array([ord(c) for c in seq2])
#     return np.sum(arr1 != arr2)

def vectorized_batch_hamming(reference_umi, umi_list):
    if not umi_list:
        return []
    
    ref_len = len(reference_umi)
    ref_array = np.frombuffer(reference_umi.encode('ascii'), dtype=np.uint8)
    ref_reshaped = ref_array.reshape(1, -1)
    
    umi_lengths = np.array([len(umi) for umi in umi_list])
    if not np.all(umi_lengths == ref_len):
        raise ValueError('Sequences must be of equal length.')
    
    umi_arrays = np.vstack(
        [np.frombuffer(umi.encode('ascii'), dtype=np.uint8) for umi in umi_list]
        )
        
    return np.sum(umi_arrays != ref_reshaped, axis=1)

@lru_cache(maxsize=10000)
def cached_hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must be of equal length.')
    return sum(x1 != x2 for x1, x2 in zip(seq1, seq2))

class OptimizedBKTree:
    def __init__(self, distance_fn, initial_item=None):
        self.distance_fn = distance_fn
        self.tree = {} if initial_item is None else {initial_item: {}}
        
    def add(self, item):
        if not self.tree:
            self.tree[item] = {}
            return
            
        node = next(iter(self.tree))
        while True:
            dist = self.distance_fn(item, node)
            if dist not in self.tree[node]:
                self.tree[node][dist] = item
                self.tree[item] = {}
                break
            node = self.tree[node][dist]
    
    def find(self, item, threshold):
        if not self.tree:
            return []
            
        results = []
        root = next(iter(self.tree))
        
        def search(node, target, max_dist):
            dist = self.distance_fn(target, node)
            if dist <= max_dist:
                results.append((dist, node))
                
            start_dist = max(0, dist - max_dist)
            end_dist = dist + max_dist
            
            for d in range(start_dist, end_dist + 1):
                if d in self.tree[node]:
                    search(self.tree[node][d], target, max_dist)
                    
        search(root, item, threshold)
        return sorted(results)

def cluster_entries_by_umis(entries, threshold, frag_ratio):
    all_n_umis = all(entry.get('ltr_umi') == 'N' and entry.get('linker_umi') == 'N' 
                     for entry in entries)
    
    if all_n_umis:
        return [entries]

    for entry in entries:
        if 'combined_umi' not in entry:
            entry['combined_umi'] = entry['ltr_umi'] + entry['linker_umi']
    
    if len(entries) > 100: 
        entries.sort(key=itemgetter('count', 'map_qual', 'mean_qual'), reverse=True)
        
        first_umi = entries[0]['combined_umi']
        tree = OptimizedBKTree(cached_hamming_distance, first_umi)
        
        clusters = defaultdict(list)
        clusters[0].append(entries[0])
        processed = {first_umi}
        
        for entry in entries[1:]:
            if entry['count'] < (frag_ratio - 1):
                cluster_id = len(clusters)
                clusters[cluster_id].append(entry)
                continue
                
            umi = entry['combined_umi']
            if umi not in processed:
                matches = tree.find(umi, threshold)
                if matches:
                    cluster_id = matches[0][1]
                else:
                    cluster_id = len(clusters)
                    tree.add(umi)
                clusters[cluster_id].append(entry)
                processed.add(umi)
    else:
        clusters = {0: entries}
    
    networks = build_umi_networks(clusters.values(), frag_ratio)
    connected_components = find_connected_components(networks)
    
    entry_clusters = []
    entry_map = {entry['read_name']: entry for entry in entries}
    processed_reads = set()
    
    for component in connected_components:
        entry_clusters.append([entry_map[read_name] for read_name in component])
        processed_reads.update(component)
    
    for read_name, entry in entry_map.items():
        if read_name not in processed_reads:
            entry_clusters.append([entry])
    
    return entry_clusters

def build_umi_networks(umi_clusters, frag_ratio, mem_aware = False):
    networks = defaultdict(list)
    unnetworked = set()
    count_key = itemgetter('count', 'map_qual', 'mean_qual')
    
    for cluster in umi_clusters:
        all_reads = {read['read_name'] for read in cluster}
        networked_set = set()
        processed_set = set()
        sorted_entries = sorted(cluster, key=count_key, reverse=True)
        
        # Creating and storing a variant dictionary for all single-difference UMIs (O(n) + O(1))
        # is much faster than on-the-fly Hamming distance calculations (O(n^2)).
        # However, it's also more memory intensive.
        # The mem_aware flag will toggle between variant creation and Hamming distance, but is not currently used.
        if not mem_aware:
            umi_variants = {}
            
            for entry in sorted_entries:
                umi = entry['combined_umi']
                umi_len = len(umi)
                umi_array = np.array(list(umi))
                
                variants = {umi}
                
                for pos in range(umi_len):
                    original_char = umi[pos]
                    if original_char != 'N':
                        for char in 'ACGT':
                            if char != original_char:
                                variant_array = umi_array.copy()
                                variant_array[pos] = char
                                variant = ''.join(variant_array)
                                variants.add(variant)
                
                umi_variants[umi] = variants
        
        for i, entry in enumerate(sorted_entries):
            if entry['read_name'] in processed_set:
                continue
                
            processed_set.add(entry['read_name'])
            current_count = entry['count']
            current_umi = entry['combined_umi']
            
            candidates = []
            
            for next_entry in sorted_entries[i+1:]:
                next_read = next_entry['read_name']
                if next_read in processed_set:
                    continue
                
                next_count = next_entry['count']
                next_umi = next_entry['combined_umi']
                
                if current_count < (frag_ratio * next_count) - 1:
                    processed_set.add(next_read)
                    continue
                
                if not mem_aware:
                    if next_umi in umi_variants[current_umi]:
                        networks[entry['read_name']].append(next_read)
                        networked_set.add(next_read)
                        networked_set.add(entry['read_name'])
                else:
                    candidates.append(next_entry)
            
            if mem_aware:  
                if candidates:
                    candidate_umis = [c['combined_umi'] for c in candidates]
                    distances = vectorized_batch_hamming(current_umi, candidate_umis)
                    
                    for idx, dist in enumerate(distances):
                        if dist <= 1:
                            next_read = candidates[idx]['read_name']
                            networks[entry['read_name']].append(next_read)
                            networked_set.add(next_read)
            
            if networks[entry['read_name']]:
                networked_set.add(entry['read_name'])
        
        unnetworked.update(all_reads - networked_set)
    
    for read in unnetworked:
        networks[read] = []
        
    return networks

def find_connected_components(networks):
    all_nodes = set(networks.keys()) | {node for nodes in networks.values() for node in nodes}
    node_to_id = {node: i for i, node in enumerate(all_nodes)}
    
    parent = list(range(len(all_nodes)))
    rank = [0] * len(all_nodes)
    
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            if rank[root_x] < rank[root_y]:
                parent[root_x] = root_y
            else:
                parent[root_y] = root_x
                if rank[root_x] == rank[root_y]:
                    rank[root_x] += 1
    
    for node, neighbors in networks.items():
        node_id = node_to_id[node]
        for neighbor in neighbors:
            neighbor_id = node_to_id[neighbor]
            union(node_id, neighbor_id)
    
    components = defaultdict(list)
    for node in all_nodes:
        node_id = node_to_id[node]
        root = find(node_id)
        components[root].append(node)
    
    return sorted(components.values(), key=len, reverse=True)

def unique_fuzzy_matches(input_dict, len_diff, umi_diff, frag_ratio):
    tmp_fuzzy_kept = defaultdict(list)
    tmp_fuzzy_dup = defaultdict(set)
    
    position_key = itemgetter('start', 'end')
    quality_key = itemgetter('count', 'map_qual', 'mean_qual')

    for key, entries in input_dict.items():        
        position_groups = ranged_groupby(entries, len_diff, sort_key = position_key)
        
        for group in position_groups:
            entry_clusters = cluster_entries_by_umis(
                entries = group, 
                threshold = umi_diff,
                frag_ratio = frag_ratio
                )

            for connected_entries in entry_clusters:
                kept_entry = max(connected_entries, key = quality_key)
                kept_read = kept_entry['read_name']
                
                duplicate_reads = set(e['read_name'] for e in connected_entries) - {kept_read}
                kept_entry['count'] += sum(
                    e['count'] for e in connected_entries 
                    if e['read_name'] in duplicate_reads
                    )
                
                tmp_fuzzy_dup[kept_read].update(duplicate_reads)
                tmp_fuzzy_kept[key].append(kept_entry)
                
    tmp_fuzzy_dup = {k: sorted(list(v)) for k, v in tmp_fuzzy_dup.items()}
    return dict(tmp_fuzzy_kept), tmp_fuzzy_dup

def multi_exact_matches(input_dict):
    multi_exact_kept = {}
    
    multi_keys = list(input_dict.keys())
    multi_values = list(input_dict.values())
    
    seq_array = np.array([
        value['ltr_umi'] + value['linker_umi'] + value['seq1']
        for value in multi_values
    ])
    
    unique_indices, inverse_indices, counts = np.unique(seq_array, return_inverse = True, return_counts = True)
    qual_key = itemgetter('map_qual', 'mean_qual')
    
    for seq_idx in range(len(unique_indices)):
        group_indices = np.where(inverse_indices == seq_idx)[0]
        if len(group_indices) == 1:
            idx = group_indices[0]
            multi_exact_kept[multi_keys[idx]] = multi_values[idx]
        else:
            group_entries = [multi_values[i] for i in group_indices]
            best_entry = max(group_entries, key=qual_key)
            best_idx = group_indices[group_entries.index(best_entry)]
            
            best_entry = best_entry.copy()
            best_entry['count'] += sum(
                multi_values[i]['count'] for i in group_indices 
                if multi_values[i] != best_entry
            )
            multi_exact_kept[multi_keys[best_idx]] = best_entry
    
    return multi_exact_kept
    
def verify_sequence_groups_faiss(group, seq_sim, len_diff, nthr):
    set_faiss_threads(nthr)
    
    if len(group) <= 1:
        return [group]
    
    compare_seqs = {}
    for i, entry in enumerate(group):
        seq1 = entry['seq1']
        compare_seqs[i] = seq1
    
    k = 4
    base_to_val = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}
    dim = 4**k
    
    vectors = np.zeros((len(group), dim), dtype=np.float32)
    
    for i, seq in compare_seqs.items():
        kmer_counts = np.zeros(dim, dtype=np.float32)
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j+k]
            if all(base in base_to_val for base in kmer):
                idx = 0
                for base in kmer:
                    idx = idx * 4 + base_to_val.get(base, 0)
                kmer_counts[idx] += 1
        vectors[i] = kmer_counts
    
    norms = np.sqrt(np.sum(vectors**2, axis=1))
    norms[norms == 0] = 1
    vectors = vectors / norms[:, np.newaxis]
    
    index = faiss.IndexFlatIP(dim)
    index.add(vectors)
    
    D, I = index.search(vectors, min(20, len(group)))
    
    rows, cols = [], []
    for i in range(len(group)):
        seq_i = compare_seqs[i]
        len_i = len(seq_i)
        
        for j_idx, j in enumerate(I[i]):
            if i == j:
                continue
                
            seq_j = compare_seqs[j]
            len_j = len(seq_j)
            
            if abs(len_i - len_j) > len_diff:
                continue
                
            if D[i][j_idx] < seq_sim - 0.05:
                continue
                
            sim = seq_similarity(seq_i, seq_j)
            if sim >= seq_sim:
                rows.append(i)
                cols.append(j)
    
    graph = csr_matrix(
        ([1] * len(rows), (rows, cols)),
        shape=(len(group), len(group)),
        dtype=np.int8
    )
    graph = graph.maximum(graph.T)
    
    n_components, labels = connected_components(
        csgraph=graph,
        directed=False,
        return_labels=True
    )
    
    result = []
    for component_id in range(n_components):
        component_indices = np.where(labels == component_id)[0]
        result.append([group[idx] for idx in component_indices])
    
    return result

def verify_sequence_groups(group, seq_sim, len_diff, nthr, min_parallel_size=500):
    if len(group) <= 1:
        return [group]
    
    if len(group) >= min_parallel_size:
        return verify_sequence_groups_faiss(group, seq_sim, len_diff, nthr)
    
    seq1_key = itemgetter('seq1')
    
    compare_seqs = {}
    for i, entry in enumerate(group):
        seq1 = seq1_key(entry)
        compare_seqs[i] = seq1
    
    entry_indices = list(range(len(group)))
    sorted_indices = sorted(entry_indices, key=lambda i: (len(compare_seqs[i]), group[i]['read_name']), reverse=True)
    
    parent = list(range(len(group)))
    rank = [0] * len(group)
    
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y: 
            if rank[root_x] < rank[root_y]:
                parent[root_x] = root_y
            else:
                parent[root_y] = root_x
                if rank[root_x] == rank[root_y]:
                    rank[root_x] += 1
    
    for i in range(len(sorted_indices)):
        idx_i = sorted_indices[i]
        seq_i = compare_seqs[idx_i]
        len_i = len(seq_i)
        
        for j in range(i):
            idx_j = sorted_indices[j]
            seq_j = compare_seqs[idx_j]
            len_j = len(seq_j)
            
            if abs(len_i - len_j) > len_diff:
                continue
            
            sim = seq_similarity(seq_i, seq_j)
            if sim >= seq_sim:
                union(idx_i, idx_j)
    
    subgroups = defaultdict(list)
    for i, idx in enumerate(sorted_indices):
        root = find(idx)
        subgroups[root].append(group[idx])
    
    return list(subgroups.values())
    
def process_group(group, umi_diff, frag_ratio, seq_sim, len_diff, min_parallel_size, nthr):
        kept_reads = {}
        dup_reads = {}
        
        verified_subgroups = verify_sequence_groups(
            group = group, 
            seq_sim = seq_sim, 
            len_diff = len_diff, 
            min_parallel_size = min_parallel_size, 
            nthr = nthr
            )
        
        for subgroup in verified_subgroups:
            entry_clusters = cluster_entries_by_umis(
                entries = subgroup,
                threshold = umi_diff,
                frag_ratio = frag_ratio
                )
            
            for connected_entries in entry_clusters:
                kept_entry = max(connected_entries, key = itemgetter('count', 'map_qual', 'mean_qual'))
                kept_read = kept_entry['read_name']
                
                duplicate_reads = set(e['read_name'] for e in connected_entries) - {kept_read}
                kept_entry['count'] += sum(
                    e['count'] for e in connected_entries 
                    if e['read_name'] in duplicate_reads
                    )
                
                kept_reads[kept_read] = kept_entry
                dup_reads[kept_read] = sorted(list(duplicate_reads))
        
        return kept_reads, dup_reads
    
def multi_fuzzy_matches(groups, umi_diff, frag_ratio, nthr, seq_sim, len_diff, min_parallel_size=500):
    
    sorted_groups = sorted(groups, key=lambda g: (len(g), min(e['read_name'] for e in g) if g else ''), reverse=True)
    large_groups = [g for g in sorted_groups if len(g) >= min_parallel_size]
    regular_groups = [g for g in sorted_groups if len(g) < min_parallel_size]
    
    multi_fuzzy_kept = {}
    multi_fuzzy_dup = {}
    
    if large_groups:        
        for large_group in large_groups:
            kept, dup = process_group(
                group = large_group,
                seq_sim = seq_sim,
                umi_diff = umi_diff,
                frag_ratio = frag_ratio,
                len_diff = len_diff,
                min_parallel_size = min_parallel_size,
                nthr = nthr
                )
            
            multi_fuzzy_kept.update(kept)
            multi_fuzzy_dup.update(dup)
    
    if regular_groups:
        total_groups = len(regular_groups)
        groups_per_batch = max(1, total_groups // nthr)
        batched_groups = [regular_groups[i:i+groups_per_batch] for i in range(0, total_groups, groups_per_batch)]
        
        results = Parallel(n_jobs=nthr)(
            delayed(lambda batch: [process_group(
                group=group,
                seq_sim=seq_sim,
                umi_diff=umi_diff,
                frag_ratio=frag_ratio,
                len_diff=len_diff,
                min_parallel_size=min_parallel_size,
                nthr = nthr
            ) for group in batch])(batch)
            for batch in batched_groups
        )

        for batch_results in results:
            for kept, dup in batch_results:
                multi_fuzzy_kept.update(kept)
                multi_fuzzy_dup.update(dup)
    
    return multi_fuzzy_kept, multi_fuzzy_dup

def stable_hash(s):
    # bits = platform.architecture()[0]
    # n_digits = 16 if bits == '64bit' else 8
    n_digits = 8
    raw_hash = int(hashlib.sha256(s.encode()).hexdigest()[:n_digits], 16)
    return raw_hash

def prefix_weighted_minhash(sequence, num_perm, min_frag_len, token_size):
    minhash = MinHash(num_perm = num_perm, seed = 1)
    for pos, i in enumerate(range(len(sequence) - token_size + 1)):
        kmer = sequence[i:i + token_size]
        weight = 1 if pos < (min_frag_len - token_size + 1) else 0
        for _ in range(weight):
            minhash.update(kmer.encode('utf-8'))
    return LeanMinHash(minhash)

def mm_coarse_cluster(sequences, threshold, num_perm, min_frag_len, token_size):    
    lsh = MinHashLSH(threshold = threshold, num_perm = num_perm)
    
    seq_to_minhash = {}
    
    sorted_seqs = sorted(sequences, key=lambda seq: (len(seq), seq), reverse=True)
    clusters = defaultdict(list)
    
    for seq in sorted_seqs:
        minhash = prefix_weighted_minhash(seq, num_perm, min_frag_len, token_size)
        seq_to_minhash[seq] = minhash
        
        matches = lsh.query(minhash) if len(lsh.keys) > 0 else []
        
        if not matches:
            clusters[seq].append(seq)
            lsh.insert(seq, minhash)
        else:
            sorted_matches = sorted(matches, key=lambda x: (len(x), x))
            best_match = sorted_matches[-1]
            clusters[best_match].append(seq)
            
    return list(clusters.values())

def refine_coarse_clusters(coarse_clusters, seq_sim):
    refined_clusters = []
    for cluster in coarse_clusters:
        if len(cluster) == 1:
            refined_clusters.append([cluster[0]])
            continue
        subclusters = defaultdict(list)
        sorted_seqs = sorted(cluster, key=lambda seq: (len(seq), seq), reverse=True)
        subclusters[sorted_seqs[0]].append(sorted_seqs[0])
        
        for seq in sorted_seqs[1:]:
            seq_len = len(seq)
            best_similarity = 0
            best_key = None
            
            for key in sorted(subclusters.keys(), key=lambda k: (len(k), k)):
                sim_check = seq_similarity(seq, key[:seq_len])
                if sim_check > best_similarity:
                    best_similarity = sim_check
                    best_key = key
            
            if best_similarity >= seq_sim:
                subclusters[best_key].append(seq)
            else:
                subclusters[seq].append(seq)
        
        refined_clusters.extend(subclusters.values())
        
    return refined_clusters

def group_mm_sequences(mm_reads, seq_sim, min_frag_len, num_perm, token_size):
        
    seq_to_entries = {}
    seqs = []
    for entry in mm_reads.values():
        seq = entry['seq1']
        if seq not in seq_to_entries:
            seq_to_entries[seq] = []
        seq_to_entries[seq].append(entry)
    seqs = list(seq_to_entries.keys())

    coarse_clusters = mm_coarse_cluster(
        sequences = seqs,
        threshold = 0.5,
        num_perm = num_perm,
        min_frag_len = min_frag_len,
        token_size = token_size
        )

    refined_seq_clusters = refine_coarse_clusters(
        coarse_clusters = coarse_clusters, 
        seq_sim = seq_sim
        )
        
    del coarse_clusters
    
    final_clusters = []
    for seq_cluster in refined_seq_clusters:
        dict_cluster = []
        for seq in seq_cluster:
            dict_cluster.extend(seq_to_entries[seq])
        final_clusters.append(dict_cluster)
    
    return final_clusters

def build_position_based_index(um_kept_dict, k, nthr):
    um_index = {}
    
    estimated_kmers = sum(1 for read in um_kept_dict.values())
    bloom = BloomFilter(capacity = estimated_kmers * 2, error_rate = 0.01)
    
    stratified_reads = defaultdict(list)
    for read_name, read in um_kept_dict.items():
        key = (read['chrom'], read['strand'])
        stratified_reads[key].append(read)
    
    def process_stratum(stratum_key, reads):
        chrom, strand = stratum_key
        local_um_index = {}
        local_prefix_set = set()
        position_fragment_count = defaultdict(int)
        
        pos_groups = defaultdict(list)
        for read in reads:
            pos = read['end'] if strand == '-' else read['start']
            pos_groups[pos].append(read)
        
        for pos, pos_reads in pos_groups.items():
            sorted_reads = sorted(pos_reads, key=lambda x: (len(x['seq1']), x['seq1']), reverse=False)
            for read in sorted_reads:
                seq = read['seq1']
                prefix = seq[:k]
                if prefix not in local_prefix_set:
                    position_info = (chrom, pos, strand)
                    read_id = read['read_name']
                    local_um_index[prefix] = (position_info, read_id)
                    local_prefix_set.add(prefix)
                position_fragment_count[(chrom, pos, strand)] += 1
    
        return local_um_index, local_prefix_set, position_fragment_count
    
    results = Parallel(n_jobs=nthr)(
        delayed(process_stratum)(key, reads)
        for key, reads in stratified_reads.items()
        )
    
    um_index = defaultdict(list)
    fragment_counts = defaultdict(int)
    for local_index, prefix_set, position_fragment_count in results:
        for prefix, (position_info, read_id) in local_index.items():
            um_index[prefix].append((position_info, read_id))
        for prefix in prefix_set:
            bloom.add(prefix)
        for position, count in position_fragment_count.items():
            fragment_counts[position] += count
    
    return um_index, bloom, fragment_counts

def compare_to_um(mm_group, k, um_index, bloom_filter, um_kept_dict, seq_sim, fragment_counts):
    sorted_mm_reads = sorted(mm_group, key = lambda x: (len(x['seq1']), x['seq1']), reverse = True)
    relocation_info = {}
    
    for mm_read in sorted_mm_reads:
        mm_seq = mm_read['seq1']        
        prefix = mm_seq[:k]
        
        if prefix not in bloom_filter:
            continue
        
        if prefix in um_index:
            positions = um_index[prefix]
            best_similarity = 0
            best_position_info = None
            best_fragment_count = 0
            
            for position_info, um_read_id in positions:
                chrom, pos, strand = position_info
                um_rep_read = um_kept_dict[um_read_id]
                um_seq = um_rep_read['seq1']
            
                if len(um_seq) < len(mm_seq):
                    continue

                common_length = min(len(mm_seq), len(um_seq))
                similarity = seq_similarity(mm_seq[:common_length], um_seq[:common_length])
                
                fragment_count = fragment_counts.get(position_info, 0)
                
                if similarity > best_similarity or (similarity == best_similarity and fragment_count > best_fragment_count):
                    best_similarity = similarity
                    best_position_info = position_info
                    best_fragment_count = fragment_count
                    
                # if best_similarity == 1:
                #     break
                
            if best_similarity >= seq_sim:
                relocation_info[mm_read['read_name']] = best_position_info
        
    if relocation_info:        
        relocated_group = []
        unrelocated_group = []
        for read in mm_group:
            read_name = read['read_name']
            if read_name in relocation_info:
                position_info = relocation_info[read_name]
                chrom, pos, strand = position_info
                relocated_read = read.copy()
                read_len = int(relocated_read['end']) - int(relocated_read['start'])
                relocated_read['chrom'] = str(chrom)
                relocated_read['start'] = int(pos) if strand == '+' else int(pos - read_len)
                relocated_read['end'] = int(pos + read_len) if strand == '+' else int(pos)
                relocated_read['strand'] = str(strand)
                relocated_read['multi'] = 'True - relocated'
                relocated_group.append(relocated_read)
            else:
                unrelocated_group.append(read)
        return (relocated_group, unrelocated_group)
    
    return mm_group

def assign_mm_group(mm_group, mm_group_threshold):
    if len(mm_group) == 1:
        return (mm_group, 1)
    
    if len(mm_group) > 1 and len(mm_group) < mm_group_threshold:
        return ([read for read in mm_group], 2)
    
    positions = np.array([(read['chrom'], read['start'], read['end'], read['strand']) for read in mm_group])
    counts = np.array([read['count'] for read in mm_group])
    
    unique_positions, position_indices, position_counts = np.unique(
        positions, axis=0, return_inverse=True, return_counts=True
        )
    
    total_counts = np.bincount(position_indices, weights=counts)
    frequencies = total_counts / total_counts.sum()
    
    min_freq = 1/(len(unique_positions) * 2)
    valid_positions = unique_positions[frequencies >= min_freq]
    
    if len(valid_positions) > 0:
        position_strings = [f"{p[0]}_{p[1] if p[3] == '+' else p[2]}_{p[3]}" for p in valid_positions]
        seed_string = ','.join(sorted(position_strings)[:2])
        seed = int(hashlib.sha256(seed_string.encode('utf-8')).hexdigest(), 16) % (2**32)
        
        random.seed(seed)
        chosen_pos = random.choices(valid_positions, k=1)[0]
        random.seed()
        
        assigned_group = []
        for read in mm_group:
            read = read.copy()
            chosen_strand = str(chosen_pos[3])
            read_len = len(read['seq1'])
            read['chrom'] = str(chosen_pos[0])
            read['start'] = int(chosen_pos[1]) if chosen_strand == '+' else (int(chosen_pos[2]) - read_len)
            read['end'] = (int(chosen_pos[1]) + read_len) if chosen_strand == '+' else int(chosen_pos[2])
            read['strand'] = chosen_strand
            read['multi'] = 'True - grouped'
            assigned_group.append(read)
            
        return (assigned_group, 0)
    
    return (mm_group, 0)

def batch_compare_to_um(chunk, k, um_index, bloom_filter, um_kept_dict, seq_sim, fragment_counts):
    relocated_reads = {}
    remaining_groups = []
    
    for group in chunk:
        result = compare_to_um(
            mm_group=group,
            k = k,
            um_index = um_index,
            bloom_filter = bloom_filter,
            um_kept_dict = um_kept_dict,
            seq_sim = seq_sim,
            fragment_counts = fragment_counts
            )
        
        if isinstance(result, tuple) and len(result) == 2:
            relocated_group, unrelocated_group = result
            for read in relocated_group:
                relocated_reads[read['read_name']] = read
            if unrelocated_group:
                remaining_groups.append(unrelocated_group)
        else:
            remaining_groups.append(result)
            
    return relocated_reads, remaining_groups

def verify_mm_positions(mm_kept_dict, um_kept_dict, seq_sim, nthr, len_diff, 
                        k, min_frag_len, num_perm, token_size, mm_group_threshold):
    
    n_mm_kept = len(mm_kept_dict)
    
    if not um_kept_dict:
        um_index, bloom_filter, fragment_counts = {}, None, {}
    else:
        um_index, bloom_filter, fragment_counts = build_position_based_index(um_kept_dict, k, nthr)
    
    subgroups = group_mm_sequences(
        mm_reads=mm_kept_dict, 
        seq_sim=seq_sim, 
        min_frag_len=min_frag_len,
        num_perm=num_perm,
        token_size=token_size
        )
    
    print(f'Number of multimapping subgroups: {len(subgroups)}', flush=True)  
    
    nsub = 0
    for group in subgroups:
        nsub += len(group)

    if nsub != len(mm_kept_dict):
        raise Exception(f'Failure in grouping multimapping reads.')
    
    sorted_subgroups = sorted(subgroups, 
                            key=lambda g: (len(g), min(e['read_name'] for e in g) if g else ''), 
                            reverse=True)
    
    relocated_reads = {}
    remaining_groups = []
    
    if um_kept_dict:
        subgroup_chunk_size = max(1, len(sorted_subgroups) // nthr)
        subgroup_chunks = [
            sorted_subgroups[i:i + subgroup_chunk_size] for i in range(0, len(sorted_subgroups), subgroup_chunk_size)
        ]
        
        relocated_results = Parallel(n_jobs=nthr)(
            delayed(batch_compare_to_um)(
                chunk=chunk,
                k=k,
                um_index=um_index,
                bloom_filter=bloom_filter,
                um_kept_dict=um_kept_dict,
                seq_sim=seq_sim,
                fragment_counts = fragment_counts
            ) for chunk in subgroup_chunks
        )

        for chunk_relocated_reads, chunk_remaining_groups in relocated_results:
            relocated_reads.update(chunk_relocated_reads)
            remaining_groups.extend(chunk_remaining_groups)
    else:
        remaining_groups = sorted_subgroups

    n_relocated = len(relocated_reads)
    n_relocated_perc = (n_relocated / n_mm_kept) * 100
    print(f'Number of multimapping reads relocated to uniquely mapped positions: {n_relocated} ({n_relocated_perc:.2f}%)',
            flush=True)

    total_singles = 0
    total_reassigned = 0
    if remaining_groups:
        reassigned_results = Parallel(n_jobs=nthr)(
            delayed(assign_mm_group)(group, mm_group_threshold) for group in remaining_groups
        )

        for group, n_single in reassigned_results:
            if n_single == 1 or n_single == 2:
                total_singles += len(group)
            else:
                total_reassigned += len(group)
            for read in group:
                relocated_reads[read['read_name']] = read
    
    reassigned_perc = (total_reassigned / n_mm_kept) * 100
    ungrouped_perc = (total_singles / n_mm_kept) * 100
    print(f'Number of ungrouped multimapping reads: {total_singles} ({ungrouped_perc:.2f}%)')
    print(f'Number of multimapping reads grouped and reassigned: {total_reassigned} ({reassigned_perc:.2f}%)')

    return relocated_reads