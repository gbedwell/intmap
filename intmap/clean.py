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
        
# def levenshtein_wrapper(a, b):
#     return Levenshtein.distance(a, b)

# def build_um_index(um_kept_dict, prefix_length):
#     index = defaultdict(list)
#     for read_name, um_frag in um_kept_dict.items():
#         prefix = um_frag['seq1'][:prefix_length].encode()
#         index[prefix].append(um_frag)
    
#     unique_prefixes = sorted(set(index.keys()))
#     tree = pybktree.BKTree(levenshtein_wrapper, unique_prefixes)
    
#     return index, tree

# def get_um_candidates(prefix, um_index, um_tree, max_distance):
#     matches = um_tree.find(prefix, max_distance)
#     keys = [key for _, key in matches]    
#     candidates = []
#     for key in keys:
#         candidates.extend(um_index.get(key, []))
    
#     return candidates

def verify_sequence_groups(group, seq_sim):
    n_seqs = len(group)
    if n_seqs <= 1:
        return [group]
    
    seq1_key = itemgetter('seq1')
    seq2_key = itemgetter('seq2')
    
    subgroups = [[group[0]]]
    
    for entry in group[1:]:
        entry_seq1 = seq1_key(entry)
        entry_seq2 = seq2_key(entry)
        
        if seq_similarity(entry_seq1, revcomp(entry_seq2)) >= seq_sim:
            compare_seq = entry_seq1
        else:
            compare_seq = entry_seq1 + entry_seq2
            
        matched = False
        for subgroup in subgroups:
            ref_seq1 = seq1_key(subgroup[0])
            ref_seq2 = seq2_key(subgroup[0])
            
            if seq_similarity(ref_seq1, revcomp(ref_seq2)) >= seq_sim:
                ref_seq = ref_seq1
            else:
                ref_seq = ref_seq1 + ref_seq2
                
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

def build_kmer_index(sequences, k = 15, step = 5):
    kmer_index = defaultdict(set)
    for i, seq in enumerate(sequences):
        for j in range(0, len(seq) - k + 1, step):
            kmer = seq[j:j+k]
            kmer_index[kmer].add(i)
    return kmer_index

def group_mm_sequences(group, seq_sim, k = 15, step = 5):
    sorted_group = sorted(group, key=lambda x: len(x['seq1']))
    
    # Build kmer index for sequences
    seqs = [entry['seq1'] for entry in sorted_group]
    kmer_index = build_kmer_index(seqs, k = k, step = step)
    
    # Pre-compute candidate matches for each sequence
    seq_candidates = {}
    for i, entry in enumerate(sorted_group):
        entry_seq = entry['seq1']
        # Count occurrences of each candidate from k-mer matches
        candidate_counts = Counter()
        for j in range(0, len(entry_seq) - k + 1, step):
            kmer = entry_seq[j:j+k]
            candidate_counts.update(kmer_index[kmer])
        
        # Keep candidates that share enough k-mers
        n_kmers = len(range(0, len(entry_seq) - k + 1, step))
        min_shared = max(1, n_kmers // 4)  # Require ~1/4 of k-mers to be shared
        seq_candidates[i] = {idx for idx, count in candidate_counts.items() if count >= min_shared}
    
    subgroups = [[sorted_group[0]]]
    for i, entry in enumerate(sorted_group[1:], 1):
        entry_seq = entry['seq1']
        matched = False
        
        # Use pre-computed candidates
        candidates = seq_candidates[i]
        
        for subgroup in subgroups:
            ref_idx = sorted_group.index(subgroup[0])
            if ref_idx in candidates:
                ref_seq = subgroup[0]['seq1']
                shorter_seq = min(entry_seq, ref_seq, key=len)
                longer_seq = max(entry_seq, ref_seq, key=len)
                
                max_similarity = max(
                    seq_similarity(shorter_seq, longer_seq[i:i+len(shorter_seq)])
                    for i in range(len(longer_seq) - len(shorter_seq) + 1)
                )

                if max_similarity >= seq_sim:
                    subgroup.append(entry)
                    matched = True
                    break
                    
        if not matched:
            subgroups.append([entry])
    
    return subgroups

def build_position_based_index(um_kept_dict, k = 15, step = 5):
    # Group UM reads by integration site
    position_groups = defaultdict(list)
    for read in um_kept_dict.values():
        if read['strand'] == '+':
            pos_key = (read['chrom'], read['start'], read['strand'])
        else:
            pos_key = (read['chrom'], read['end'], read['strand'])
        position_groups[pos_key].append(read)
    
    # Build k-mer index using one sequence per integration site
    kmer_index = defaultdict(set)
    for pos, reads in position_groups.items():
        rep_seq = max(reads, key = lambda x: len(x['seq1']))['seq1']
        for i in range(0, len(rep_seq) - k + 1, step):
            kmer = rep_seq[i:i+k]
            kmer_index[kmer].add(pos)
            
    return kmer_index, position_groups

def compare_to_um(group, um_index, position_groups, seq_sim, k = 15, step = 5):
    longest_read = max(group, key=lambda x: len(x['seq1']))
    read_seq = longest_read['seq1']
    
    um_candidates = Counter()
    for j in range(0, len(read_seq) - k + 1, step):
        kmer = read_seq[j:j+k]
        um_candidates.update(um_index[kmer])
        
    for pos, count in um_candidates.items():
        um_reads = position_groups[pos]
        um_seq = min(um_reads, key=lambda x: len(x['seq1']))['seq1']
        
        if len(read_seq) <= len(um_seq):
            max_similarity = max(
                seq_similarity(read_seq, um_seq[i:i+len(read_seq)])
                for i in range(len(um_seq) - len(read_seq) + 1)
            )
            
            if max_similarity >= seq_sim:
                relocated_group = []
                for read in group:
                    read = read.copy()
                    read['chrom'] = pos[0]
                    read['start'] = pos[1] if pos[2] == '+' else (pos[1] - len(read['seq1']))
                    read['end'] = (pos[1] + len(read['seq1'])) if pos[2] == '+' else pos[1]
                    read['strand'] = pos[2]
                    read['multi'] = 'True - relocated'
                    relocated_group.append(read)
                return relocated_group
    return None

def assign_mm_group(group):
    positions = defaultdict(int)
    for read in group:
        pos_key = (read['chrom'], read['start'], read['end'], read['strand'])
        positions[pos_key] += read['count']
    
    total_counts = sum(positions.values())
    min_freq = 1/(len(positions) * 2)
    filtered_positions = {
        pos: count for pos, count in positions.items() 
        if count/total_counts >= min_freq
        }
    
    seed_string = ','.join(str(pos) for pos in sorted(filtered_positions.keys())[:2])
    seed = int(hashlib.sha256(seed_string.encode('utf-8')).hexdigest(), 16) % (2**64)
    
    random.seed(seed)
    chosen_pos = random.choices(
        population=list(filtered_positions.keys()),
        k=1
        )[0]
    
    assigned_group = []
    for read in group:
        read = read.copy()
        read['chrom'], read['start'], read['end'], read['strand'] = chosen_pos
        read['multi'] = 'True - grouped'
        assigned_group.append(read)
        
    return assigned_group

def verify_mm_positions(mm_kept_dict, um_kept_dict, seq_sim, nthr, k = 15, step = 5):
    # Build position-based index for UM reads
    um_index, position_groups = build_position_based_index(um_kept_dict, k = k, step = step)
    
    # Group similar MM reads
    mm_reads = list(mm_kept_dict.values())
    subgroups = group_mm_sequences(mm_reads, seq_sim, k = k, step = step)
    
    # Compare MM groups to UM reads
    relocated_results = Parallel(n_jobs=nthr)(
        delayed(compare_to_um)(group, um_index, position_groups, seq_sim, k, step)
        for group in subgroups
        )
    
    relocated_reads = {}
    remaining_groups = []
    
    for i, result in enumerate(relocated_results):
        if result:
            for read in result:
                relocated_reads[read['read_name']] = read
        else:
            remaining_groups.append(subgroups[i])

    if remaining_groups:
        reassigned_results = Parallel(n_jobs=nthr)(
            delayed(assign_mm_group)(group) 
            for consensus_id, group in enumerate(remaining_groups)
        )
        
        for group in reassigned_results:
            for read in group:
                relocated_reads[read['read_name']] = read
    
    mm_kept_dict.update(relocated_reads)
    return mm_kept_dict





