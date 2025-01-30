import math
from rapidfuzz.distance import Levenshtein
from collections import defaultdict, deque, Counter
import numpy as np
import random
from itertools import chain, groupby, combinations
from operator import itemgetter
from joblib import Parallel, delayed
import multiprocessing
import pybktree
import random
import hashlib
import faiss
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from intmap.utils import *

# Calculate Hamming distance
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return ValueError('Sequences must be of equal length.')
    return sum(x1 != x2 for x1, x2 in zip(seq1, seq2))

# Calculate sequence similarity based on maximum sequence length
# and the Levenshtein distance
def seq_similarity(seq1, seq2):
    return 1 - Levenshtein.normalized_distance(seq1, seq2)

def unique_exact_matches(input_dict):
    tmp_exact_kept = defaultdict(list)
    tmp_exact_dup = set()
    
    sort_key = itemgetter('start', 'end', 'ltr_umi', 'linker_umi')
    qual_key = itemgetter('mean_qual')
    
    for key, entries in input_dict.items():
        entries.sort(key = sort_key)
        
        for _, group_entries in groupby(entries, key = sort_key):
            group = list(group_entries)
            
            # Skip if any read in group is already marked as duplicate
            if not tmp_exact_dup.intersection(entry['read_name'] for entry in group):
                if len(group) > 1:
                    best_entry = max(group, key=qual_key)
                    best_entry['count'] += sum(
                        entry['count'] for entry in group 
                        if entry is not best_entry
                    )
                    tmp_exact_kept[key].append(best_entry)
                    
                    # Add non-best entries to duplicates
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

def cluster_entries_by_umis(entries, threshold, frag_ratio):
    if len(entries) > 1000:
        entries.sort(key=itemgetter('count'), reverse = True)
        
        first_umi = entries[0]['ltr_umi'] + entries[0]['linker_umi']
        tree = pybktree.BKTree(hamming_distance, [first_umi])
        
        clusters = defaultdict(list)
        clusters[0].append(entries[0])
        processed = {first_umi}
        
        for entry in entries[1:]:
            umi = entry['ltr_umi'] + entry['linker_umi']
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
    connected_components = traverse_umi_networks(networks, clusters.values())
    
    entry_clusters = []
    entry_map = {entry['read_name']: entry for entry in entries}
    processed_reads = set()
    
    # Add networked reads
    for component in connected_components:
        entry_clusters.append([entry_map[read_name] for read_name in component])
        processed_reads.update(component)
    
    # Add unnetworked reads as single-entry clusters
    for read_name, entry in entry_map.items():
        if read_name not in processed_reads:
            entry_clusters.append([entry])
    
    return entry_clusters

def build_umi_networks(umi_clusters, frag_ratio):
    networks = defaultdict(list)
    unnetworked = set()
    count_key = itemgetter('count')
    
    for cluster in umi_clusters:
        all_reads = {read['read_name'] for read in cluster}
        networked_set = set()
        sorted_entries = sorted(cluster, key=count_key, reverse=True)
        
        for entry in sorted_entries:
            if entry['read_name'] in networked_set:
                continue
                
            networked_set.add(entry['read_name'])  # Add initial entry to network
            to_process = deque([(entry, sorted_entries)])
            
            while to_process:
                current_entry, remaining = to_process.popleft()
                current_read = current_entry['read_name']
                current_umi = current_entry['ltr_umi'] + current_entry['linker_umi']
                current_count = current_entry['count']
                
                for next_entry in remaining:
                    next_read = next_entry['read_name']
                    if next_read in networked_set:
                        continue
                        
                    next_count = next_entry['count']
                    next_umi = next_entry['ltr_umi'] + next_entry['linker_umi']
                    
                    if hamming_distance(current_umi, next_umi) <= 1 and current_count >= (frag_ratio * next_count) - 1:
                        networks[entry['read_name']].append(next_read)
                        networked_set.add(next_read)
                        to_process.append((next_entry, remaining))
        
        unnetworked.update(all_reads - networked_set)
    
    for read in unnetworked:
        networks[read] = []
        
    return networks

def traverse_umi_networks(networks, clusters):
    read_counts = {read['read_name']: read['count'] for cluster in clusters for read in cluster}
    matched = set()
    components = []
    count_getter = read_counts.get
    sorted_nodes = sorted(networks.keys(), key = count_getter, reverse = True)
    
    for node in sorted_nodes:
        if node not in matched:
            component = bfs(node, networks)
            matched.update(component)
            components.append(sorted(component))
            
    return components

def bfs(start_node, network):
    visited = set([start_node])
    queue = deque([start_node])
    
    while queue:
        current = queue.popleft()
        for neighbor in network[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
                
    return visited

def unique_fuzzy_matches(input_dict, len_diff, umi_diff, frag_ratio):
    tmp_fuzzy_kept = defaultdict(list)
    tmp_fuzzy_dup = defaultdict(set)
    
    position_key = itemgetter('start', 'end')
    quality_key = itemgetter('count', 'mean_qual')

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
        value['ltr_umi'] + value['linker_umi'] + value['seq1'] + value['seq2']
        for value in multi_values
    ])
    
    unique_indices, inverse_indices, counts = np.unique(seq_array, return_inverse = True, return_counts = True)
    qual_key = itemgetter('mean_qual')
    
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

def chunk_groups(groups, nthr):
    n = max(1, math.ceil(len(groups) / nthr))
    for i in range(0, len(groups), n):
        yield groups[i:i + n]

def verify_sequence_groups(group, seq_sim):
    if len(group) <= 1:
        return [group]
    
    seq1_key = itemgetter('seq1')
    seq2_key = itemgetter('seq2')
    
    # Pre-compute sequences for all entries
    sequences = [(seq1_key(entry), seq2_key(entry)) for entry in group]
    
    # Pre-compute concatenated sequences when needed
    compare_seqs = {}
    for i, (seq1, seq2) in enumerate(sequences):
        if seq_similarity(seq1, revcomp(seq2)) < seq_sim:
            compare_seqs[i] = seq1 + seq2
        else:
            compare_seqs[i] = seq1
    
    subgroups = [[group[0]]]
    first_seqs = sequences[0]
    
    # Use pre-computed sequences for comparisons
    for i, entry in enumerate(group[1:], 1):
        entry_seqs = sequences[i]
        compare_seq = compare_seqs.get(i, entry_seqs[0])
            
        matched = False
        for subgroup in subgroups:
            ref_idx = group.index(subgroup[0])
            ref_seqs = sequences[ref_idx]
            ref_seq = compare_seqs.get(ref_idx, ref_seqs[0])
                
            if seq_similarity(compare_seq, ref_seq) >= seq_sim:
                subgroup.append(entry)
                matched = True
                break
                
        if not matched:
            subgroups.append([entry])
    
    return subgroups

def multi_fuzzy_matches(groups, umi_diff, frag_ratio, nthr, seq_sim):
    def process_group(group):
        kept_reads = {}
        dup_reads = {}
        
        verified_subgroups = verify_sequence_groups(group, seq_sim)
        
        for subgroup in verified_subgroups:
            entry_clusters = cluster_entries_by_umis(
                entries = subgroup,
                threshold = umi_diff,
                frag_ratio = frag_ratio
                )
            
            for connected_entries in entry_clusters:
                kept_entry = max(connected_entries, key = itemgetter('count', 'mean_qual'))
                kept_read = kept_entry['read_name']
                
                duplicate_reads = set(e['read_name'] for e in connected_entries) - {kept_read}
                kept_entry['count'] += sum(
                    e['count'] for e in connected_entries 
                    if e['read_name'] in duplicate_reads
                    )
                
                kept_reads[kept_read] = kept_entry
                dup_reads[kept_read] = sorted(list(duplicate_reads))
                
        return kept_reads, dup_reads

    results = Parallel(n_jobs = nthr)(delayed(process_group)(group) for group in groups)
    
    multi_fuzzy_kept = {}
    multi_fuzzy_dup = {}
    for kept, dup in results:
        multi_fuzzy_kept.update(kept)
        multi_fuzzy_dup.update(dup)
    
    return multi_fuzzy_kept, multi_fuzzy_dup

def group_mm_sequences(group, seq_sim, len_diff, k=20):
    sorted_group = sorted(group, key=lambda x: len(x['seq1']))
    seqs = [entry['seq1'] for entry in sorted_group]
    
    # Convert sequences to k-mer arrays
    kmer_arrays = []
    for seq in seqs:
        max_k = len(seq) - k + 1
        kmers = np.array([hash(seq[i:i+k]) for i in range(0, min(len_diff, max_k))])
        kmer_arrays.append(kmers)
    
    # Build vectorized index
    unique_kmers = np.unique(np.concatenate(kmer_arrays))
    kmer_to_idx = {kmer: idx for idx, kmer in enumerate(unique_kmers)}
    
    # Create sparse matrix of k-mer occurrences
    n_seqs = len(seqs)
    n_kmers = len(unique_kmers)
    kmer_matrix = np.zeros((n_seqs, n_kmers), dtype=np.int8)
    
    for seq_idx, kmers in enumerate(kmer_arrays):
        for kmer in kmers:
            kmer_matrix[seq_idx, kmer_to_idx[kmer]] = 1
    
    # Fast candidate identification using matrix multiplication
    shared_kmers = kmer_matrix @ kmer_matrix.T
    n_kmers_per_seq = np.array([len(kmers) for kmers in kmer_arrays])
    min_shared = 1
    
    candidates = np.where(shared_kmers >= min_shared)
    
    # Group sequences using identified candidates
    subgroups = []
    processed = set()
    
    for i in range(n_seqs):
        if i in processed:
            continue
            
        current_group = [sorted_group[i]]
        processed.add(i)
        
        candidates_i = candidates[1][candidates[0] == i]
        for j in candidates_i:
            if j not in processed:
                shorter_seq = min(seqs[i], seqs[j], key=len)
                longer_seq = max(seqs[i], seqs[j], key=len)
                max_similarity = max(
                    seq_similarity(shorter_seq, longer_seq[pos:(pos + len(shorter_seq))])
                    for pos in range(len_diff + 1)
                    )
                if max_similarity >= seq_sim:
                    current_group.append(sorted_group[j])
                    processed.add(j)
                    
        subgroups.append(current_group)
    
    return subgroups

def build_position_based_index(um_kept_dict, len_diff, nthr, k=20):
    n_reads = len(um_kept_dict)
    n_kmers = len_diff + 1
    k_mod = k - len_diff
    vectors = np.zeros((n_reads, n_kmers), dtype=np.float32)
    positions = []
    read_names = []
    
    for idx, (read_name, read) in enumerate(um_kept_dict.items()):
        seq = read['seq1']
        kmers = [hash(seq[i:i+k_mod]) for i in range(n_kmers)]
        vectors[idx] = kmers / np.linalg.norm(kmers)
        
        positions.append((read['chrom'], 
                        read['end'] if read['strand'] == '-' else read['start'], 
                        read['strand']))
        read_names.append(read_name)
    
    index = faiss.IndexFlatL2(n_kmers)
    index.add(vectors)
    
    position_array = np.array(positions, dtype=[
        ('chrom', 'U20'), 
        ('pos', np.int32), 
        ('strand', 'U1')
    ])
    
    read_array = np.array(read_names)
    
    return index, position_array, read_array

def compare_to_um(mm_group, um_index, um_positions, um_read_names, 
                    um_kept_dict, seq_sim, k, len_diff):
    longest_mm_read = max(mm_group, key=lambda x: len(x['seq1']))
    longest_mm_seq = longest_mm_read['seq1']
    
    n_kmers = len_diff + 1
    mm_vectors = np.zeros((1, n_kmers), dtype=np.float32)
    
    kmers = [hash(longest_mm_seq[i:i+k]) for i in range(n_kmers)]
    mm_vectors[0] = kmers / np.linalg.norm(kmers)
    
    k_neighbors = 5
    D, I = um_index.search(mm_vectors, k_neighbors)
    
    for dist, um_idx in zip(D[0], I[0]):
        um_pos = um_positions[um_idx]
        um_read = um_kept_dict[str(um_read_names[um_idx])]
        um_seq = um_read['seq1']
        
        if len(longest_mm_seq) >= len(um_seq):
            similarities = [
                seq_similarity(um_seq, longest_mm_seq[i:i+len(um_seq)])
                for i in range(len(longest_mm_seq) - len(um_seq) + 1)
                if i + len(um_seq) <= len(longest_mm_seq)
            ]
            
            if similarities:
                max_similarity = max(similarities)
                if max_similarity >= seq_sim:
                    relocated_group = []
                    for read in mm_group:
                        relocated_read = read.copy()
                        relocated_read['chrom'] = str(um_pos['chrom'])
                        relocated_read['start'] = int(um_pos['pos']) if um_pos['strand'] == '+' else int(um_pos['pos'] - len(longest_mm_seq))
                        relocated_read['end'] = int(um_pos['pos'] + len(longest_mm_seq)) if um_pos['strand'] == '+' else int(um_pos['pos'])
                        relocated_read['strand'] = str(um_pos['strand'])
                        relocated_read['multi'] = 'True - relocated'
                        relocated_group.append(relocated_read)
                    return relocated_group
            
    return None

def assign_mm_group(mm_group):
    if len(mm_group) == 1:
        return(mm_group)
    
    # Vectorized position counting
    positions = np.array([(read['chrom'], read['start'], read['end'], read['strand']) for read in mm_group])
    counts = np.array([read['count'] for read in mm_group])
    
    # Get unique positions and their total counts
    unique_positions, position_indices, position_counts = np.unique(
        positions, axis=0, return_inverse=True, return_counts=True
        )
    
    # Calculate position frequencies
    total_counts = np.bincount(position_indices, weights=counts)
    frequencies = total_counts / total_counts.sum()
    
    # Filter positions meeting minimum frequency
    min_freq = 1/(len(unique_positions) * 2)
    valid_positions = unique_positions[frequencies >= min_freq]
    
    if len(valid_positions) > 0:
        # Deterministic random choice based on position hash
        position_strings = [f"{p[0]}_{p[1]}_{p[2]}_{p[3]}" for p in valid_positions]
        seed_string = ','.join(sorted(position_strings)[:2])
        seed = int(hashlib.sha256(seed_string.encode('utf-8')).hexdigest(), 16) % (2**64)
        
        random.seed(seed)
        chosen_pos = random.choices(valid_positions, k=1)[0]
        
        # Assign all reads to chosen position
        assigned_group = []
        for read in mm_group:
            read = read.copy()
            read['chrom'] = str(chosen_pos[0])
            read['start'] = int(chosen_pos[1])
            read['end'] = int(chosen_pos[2])
            read['strand'] = str(chosen_pos[3])
            read['multi'] = 'True - grouped'
            assigned_group.append(read)
            
        return assigned_group
    
    return mm_group

def verify_mm_positions(mm_kept_dict, um_kept_dict, seq_sim, nthr, len_diff, k):
    print(f'Building unique read index...')
    um_index, um_positions, um_read_names = build_position_based_index(
        um_kept_dict = um_kept_dict, 
        len_diff = len_diff, 
        nthr = nthr, 
        k = k
        )
    
    print(f'Building multimapping read index...')
    mm_reads = list(mm_kept_dict.values())
    subgroups = group_mm_sequences(
        group = mm_reads, 
        seq_sim = seq_sim, 
        len_diff = len_diff, 
        k = k
        )
        
    print(f'Processing multimapping read clusters...')
    relocated_results = Parallel(n_jobs=nthr, prefer="threads")(
        delayed(compare_to_um)(
            mm_group = group, 
            um_index = um_index, 
            um_positions = um_positions, 
            um_read_names = um_read_names, 
            um_kept_dict = um_kept_dict, 
            seq_sim = seq_sim, 
            k = k, 
            len_diff = len_diff
            ) 
        for group in subgroups
    )
    
    relocated_reads = {}
    remaining_groups = []
    
    for i, relocated_group in enumerate(relocated_results):
        if relocated_group:
            for read in relocated_group:
                relocated_reads[read['read_name']] = read
        else:
            remaining_groups.append(subgroups[i])

    if remaining_groups:
        reassigned_results = Parallel(n_jobs=nthr)(
            delayed(assign_mm_group)(group) 
            for group in remaining_groups
            )
        
        for group in reassigned_results:
            for read in group:
                relocated_reads[read['read_name']] = read
    
    mm_kept_dict.update(relocated_reads)
    return mm_kept_dict





