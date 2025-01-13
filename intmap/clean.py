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
from bisect import bisect_left, bisect_right
import random
import hashlib
import faiss

def revcomp(seq):
    trans = str.maketrans("acgtumrwsykvhdbnACGTUMRWSYKVHDBN", "TGCAAKYWSRMBDHVNTGCAAKYWSRMBDHVN")
    seq_list = list(seq)
    seq_list.reverse()
    rev_seq = ''.join(seq_list)
    revcomp_seq = str.translate(rev_seq, trans)
    return revcomp_seq

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
    tmp_exact_kept = {}
    tmp_exact_dup = set()

    for key, entries in input_dict.items():
        entries.sort(key = itemgetter('start', 'end', 'ltr_umi', 'linker_umi'))
        
        for group_key, group_entries in groupby(entries, key=itemgetter('start', 'end', 'ltr_umi', 'linker_umi')):
            group = list(group_entries)

            if any(entry['read_name'] in tmp_exact_dup for entry in group):
                continue
            
            if len(group) > 1:
                best_entry = max(group, key = itemgetter('mean_qual'))
                best_entry['count'] += sum(entry['count'] for entry in group if entry != best_entry)
                tmp_exact_kept.setdefault(key, []).append(best_entry)

                tmp_exact_dup.update(entry['read_name'] for entry in group if entry != best_entry)
            else:
                tmp_exact_kept.setdefault(key, []).extend(group)

    return tmp_exact_kept

def ranged_groupby(entries, tolerance, sort_key = lambda x: (x['start'], x['end'])):
    entries = sorted(entries, key = sort_key)
    groups = []
    current_group = []
    last_start, last_end = None, None

    for entry in entries:
        start, end = sort_key(entry)
        if (current_group and
            (abs(start - last_start) > tolerance and
            abs(end - last_end) > tolerance)):
            groups.append(current_group)
            current_group = []
        current_group.append(entry)
        last_start, last_end = start, end

    if current_group:
        groups.append(current_group)

    return groups

def one_hot_encoding(umi, encoding_map, char_dict):
    indices = np.array([char_dict[char] for char in umi])
    one_hot = encoding_map[indices].flatten()
    
    return one_hot.astype(np.float32)

def cluster_umis(umis, threshold):
    encoding_map = np.array([
        [1, 0, 0, 0, 0, 0, 0, 0],  # A
        [0, 1, 0, 0, 0, 0, 0, 0],  # T
        [0, 0, 1, 0, 0, 0, 0, 0],  # G
        [0, 0, 0, 1, 0, 0, 0, 0],  # C
        [0, 0, 0, 0, 1, 0, 0, 0]   # N
        ])
        
    char_to_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 4}
    
    encoded_umis = np.array(
        [one_hot_encoding(umi = umi, encoding_map = encoding_map, char_dict = char_to_index)
            for umi in umis]
        )

    dim = len(encoded_umis[0])
    index = faiss.IndexFlatL2(dim)
    index.add(encoded_umis)
    
    # Use a graph to catch edge cases where A is like B and B is like C, but A is not like C.
    # These UMIs should probably all be grouped together and processed in the directional UMI network.
    graph = defaultdict(set)
    for i, umi in enumerate(umis):
        distances, indexes = index.search(encoded_umis[i:i+1], min(100, len(umis)))
        for idx, dist in zip(indexes[0], distances[0]):
            if dist <= threshold:
                graph[i].add(idx)
                graph[idx].add(i)

    visited = set()
    clusters = []

    def iterative_dfs(start_node):
        stack = [start_node]
        cluster = []
        visited.add(start_node)
        while stack:
            node = stack.pop()
            cluster.append(umis[node])
            for neighbor in graph[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        return cluster

    for i in range(len(umis)):
        if i not in visited:
            cluster = iterative_dfs(i)
            clusters.append(cluster)

    return clusters

def cluster_entries_by_umis(entries, threshold):
    combined_umis = [(entry['ltr_umi'] + entry['linker_umi']) for entry in entries]

    combined_clusters = cluster_umis(combined_umis, threshold)

    cluster_map = defaultdict(set)
    for idx, entry in enumerate(entries):
        combined_umi = entry['ltr_umi'] + entry['linker_umi']
        for cluster_idx, cluster in enumerate(combined_clusters):
            if combined_umi in cluster:
                cluster_map[cluster_idx].add(idx)

    umi_clusters = []
    visited_entries = set()

    for cluster_idx, cluster_indices in cluster_map.items():
        new_indices = cluster_indices - visited_entries
        if new_indices:
            umi_clusters.append([entries[idx] for idx in new_indices])
            visited_entries.update(new_indices)

    isolated_entries = set(range(len(entries))) - visited_entries
    for idx in isolated_entries:
        umi_clusters.append([entries[idx]])

    return umi_clusters

def build_umi_networks(umi_clusters, frag_ratio):
    networks = defaultdict(list)
    unnetworked = set()
    
    for cluster in umi_clusters:
        all_reads = {read['read_name'] for read in cluster}
        networked_set = set()
        for entry1, entry2 in combinations(cluster, 2):
            count1 = entry1['count']
            count2 = entry2['count']

            umi1 = entry1['ltr_umi'] + entry1['linker_umi']
            umi2 = entry2['ltr_umi'] + entry2['linker_umi']
            
            read1 = entry1['read_name']
            read2 = entry2['read_name']
            
            hd = hamming_distance(umi1, umi2)

            if hd <= 1:
                if count1 >= (frag_ratio * count2) - 1:
                    networks[read1].append(read2)
                    networked_set.add(read1)
                    networked_set.add(read2)
                elif count2 >= (frag_ratio * count1) - 1:
                    networks[read2].append(read1)
                    networked_set.add(read2)
                    networked_set.add(read1)
                    
        unnetworked.update(all_reads - networked_set)
        
    for read in unnetworked:
        networks[read] = []

    return networks

def traverse_umi_networks(networks, clusters):
    read_counts = {read['read_name']: read['count'] for cluster in clusters for read in cluster}
    matched = set()
    components = []
    sorted_nodes = sorted(networks.keys(), key = lambda x: read_counts[x], reverse = True)   
    for node in sorted_nodes:
        if node not in matched: 
            component = list(bfs(node, networks))
            component.sort()
            matched.update(component)
            components.append(component)

    return components         

def bfs(start_node, network):
    visited = set()
    queue = deque([start_node])
    queued = set([start_node])

    while queue:
        current_node = queue.popleft()
        visited.add(current_node)
        
        neighbors = network.get(current_node, [])
        for neighbor in neighbors:
            if neighbor not in visited and neighbor not in queued:
                queue.append(neighbor)
                queued.add(neighbor)
                
    return visited

def unique_fuzzy_matches(input_dict, len_diff, umi_diff, frag_ratio):
    tmp_fuzzy_kept = {}
    tmp_fuzzy_dup = {}

    for key, entries in input_dict.items():
        position_groups = ranged_groupby(entries, len_diff, sort_key=lambda x: (x['start'], x['end']))
        for group in position_groups:
            umi_clusters = cluster_entries_by_umis(entries = group, threshold = (umi_diff * 2))
            umi_networks = build_umi_networks(umi_clusters, frag_ratio)
            connected_reads = traverse_umi_networks(networks = umi_networks, clusters = umi_clusters)

            clustered_entries = {}
            for group in umi_clusters:
                for entry in group:
                    clustered_entries[entry['read_name']] = entry

            for branch in connected_reads:
                connected_entries = [clustered_entries[read_name] for read_name in branch if read_name in clustered_entries]
                kept_entry = max(connected_entries, key=lambda x: (x['count'], x['mean_qual']))
                kept_read = kept_entry['read_name']
                
                # TO-DO: Add error handling if kept_read already exists in kept_fuzzy_dup
                
                for entry in connected_entries:
                    if entry['read_name'] != kept_read:
                        kept_entry['count'] += entry['count']
                        tmp_fuzzy_dup.setdefault(kept_read, set()).add(entry['read_name'])
                
                tmp_fuzzy_kept[key] = tmp_fuzzy_kept.get(key, []) + [kept_entry]
            
    tmp_fuzzy_dup = {key: list(value) for key, value in tmp_fuzzy_dup.items()}            
    return tmp_fuzzy_kept, tmp_fuzzy_dup

def multi_exact_matches(input_dict):
    multi_exact_kept = {}

    seq_list = [
        value['ltr_umi'] + value['linker_umi'] + value['seq1'] + value['seq2']
        for value in input_dict.values()
        ]
    multi_keys = list(input_dict.keys())
    multi_values = list(input_dict.values())

    seq_array = np.array(seq_list)
    seq_elements, seq_counts = np.unique(seq_array, return_counts = True)

    unique_sequences = set(seq_elements[seq_counts == 1])
    non_unique_sequences = set(seq_elements[seq_counts > 1])

    for i, seq in enumerate(seq_array):
        if seq in unique_sequences:
            multi_exact_kept[multi_keys[i]] = multi_values[i]

    grouped_non_unique_indices = {}
    for i, seq in enumerate(seq_array):
        if seq in non_unique_sequences:
            grouped_non_unique_indices.setdefault(seq, []).append(i)

    for indices in grouped_non_unique_indices.values():
        tmp_keys = [multi_keys[i] for i in indices]
        tmp_values = [multi_values[i] for i in indices]
        
        best_entry = max(tmp_values, key=lambda x: x['mean_qual'])
        best_idx = tmp_values.index(best_entry)

        best_entry['count'] += sum(entry['count'] for entry in tmp_values if entry != best_entry)
        multi_exact_kept[tmp_keys[best_idx]] = best_entry

    return multi_exact_kept

def chunk_groups(groups, nthr):
    n = max(1, math.ceil(len(groups) / nthr))
    for i in range(0, len(groups), n):
        yield groups[i:i + n]
        
def levenshtein_wrapper(a, b):
    return Levenshtein.distance(a, b)

def build_um_index(um_kept_dict, prefix_length):
    index = defaultdict(list)
    for read_name, um_frag in um_kept_dict.items():
        prefix = um_frag['seq1'][:prefix_length]
        index[prefix].append(um_frag)
        
    tree = pybktree.BKTree(levenshtein_wrapper, list(index.keys()))
    return index, tree

def get_um_candidates(prefix, um_index, um_tree, max_distance):
    matches = um_tree.find(prefix, max_distance)
    keys = [key for _, key in matches]    
    candidates = []
    for key in keys:
        candidates.extend(um_index.get(key, []))
    
    return candidates

def verify_mm_positions(mm_frag, candidate_um_frags, seq_sim):
    mm_seq1 = mm_frag['seq1']
    mm_seq2 = mm_frag['seq2']
    mm_strand = mm_frag['strand']

    if seq_similarity(mm_seq1, revcomp(mm_seq2)) < seq_sim:
        return mm_frag, False

    for um_frag in candidate_um_frags:
        um_seq1 = um_frag['seq1']
        um_strand = um_frag['strand']

        if len(mm_seq1) > len(um_seq1):
            continue
        
        um_seq1_truncated = um_seq1[:len(mm_seq1)]
        similarity = seq_similarity(mm_seq1, um_seq1_truncated)

        if similarity >= seq_sim:
            mm_frag['start'] = um_frag['start'] if um_strand == '+' else (um_frag['end'] - len(mm_seq1))
            mm_frag['end'] = (um_frag['start'] + len(mm_seq1)) if um_strand == '+' else um_frag['end']
            mm_frag['multi'] = 'True - relocated'
            return mm_frag, True

    return mm_frag, False

def process_mm_fuzzy_group(group, umi_diff, frag_ratio, um_index, um_tree, seq_sim, 
                            prefix_length, unique_coordinates, max_distance):
    tmp_fuzzy_kept = {}
    tmp_fuzzy_dup = {}
    
    umi_clusters = cluster_entries_by_umis(entries = group, threshold = (umi_diff * 2))
    umi_networks = build_umi_networks(umi_clusters, frag_ratio)
    connected_reads = traverse_umi_networks(networks = umi_networks, clusters = umi_clusters)
    
    clustered_entries = {}
    for umi_group in umi_clusters:
        for entry in umi_group:
            clustered_entries[entry['read_name']] = entry
            
    for branch in connected_reads:
        connected_entries = [clustered_entries[read_name] for read_name in branch if read_name in clustered_entries]
        kept_entry = max(connected_entries, key=lambda x: (x['count'], x['mean_qual']))
        kept_read = kept_entry['read_name']

        for entry in connected_entries:
            if entry['read_name'] != kept_read:
                kept_entry['count'] += entry['count']
                tmp_fuzzy_dup.setdefault(kept_read, set()).add(entry['read_name'])

        tmp_fuzzy_kept[kept_entry['read_name']] = kept_entry
    
    tmp_fuzzy_dup = {key: list(value) for key, value in tmp_fuzzy_dup.items()}
    
    relocated_reads = {}
    prefix_cache = set()
    match_found = False
    for read_name, mm_frag in tmp_fuzzy_kept.items():
        if not match_found:
            prefix = mm_frag['seq1'][:prefix_length]
            if prefix not in prefix_cache:
                candidate_um_frags = get_um_candidates(
                    prefix = prefix, 
                    um_index = um_index, 
                    um_tree = um_tree,
                    max_distance = max_distance
                    )
                checked_read, match_found = verify_mm_positions(mm_frag, candidate_um_frags, seq_sim)
                relocated_reads[read_name] = checked_read
                prefix_cache.add(prefix)
            else:
                relocated_reads[read_name] = mm_frag
        else:
            relocated_reads[read_name] = mm_frag

    grouped_reads = {}
    # target_positions = {}
    # seen_keys = set()
    best_read_name = None
    for read_name, read_info in relocated_reads.items():
        if read_info['multi'] == 'True - relocated':
            best_read_name = read_name
            break
        
        # current_chrom = read_info['chrom']
        # current_pos = math.floor((read_info['start'] + read_info['end']) / 2)
        # current_key = current_chrom + '_' + str(current_pos)
        # if current_key not in seen_keys:
        #     if current_chrom in unique_coordinates:
        #         unique_positions = unique_coordinates[current_chrom][0]
        #         start = current_pos - win_len # win_len has been removed as a variable
        #         end = current_pos + win_len
        #         start_index = bisect_left(unique_positions, start)
        #         end_index = bisect_right(unique_positions, end)
        #         density = end_index - start_index
        #         target_positions[read_name] = max(density, 1e-6)
        #         seen_keys.add(current_key)
        #     else:
        #         target_positions[read_name] = 1e-6
        #         seen_keys.add(current_key)

    if best_read_name is None:
        read_names = list(relocated_reads.keys())
        # target_densities = list(target_positions.values())
        # total = sum(target_densities)
        # normalized_densities = [value / total for value in target_densities]
        seed_string = ','.join(read_names)
        seed = int(hashlib.sha256(seed_string.encode('utf-8')).hexdigest(), 16) % (2**64)
        random.seed(seed)
        
        best_read_name = random.choices(
            population = read_names,
            # weights = normalized_densities,
            k = 1
            )[0]

    best_read = relocated_reads[best_read_name]
    best_start, best_end, best_strand = best_read['start'], best_read['end'], best_read['strand']
    
    for read_name, read_info in relocated_reads.items():
        if read_info['multi'] != 'True - relocated':
            read_info['start'] = best_start
            read_info['end'] = best_end
            read_info['strand'] = best_strand
            read_info['multi'] = 'True - grouped'
            
    prefixes = [read_info['seq1'][:int(prefix_length * 1.5)] for read_info in relocated_reads.values()]
    transposed_prefixes = list(zip(*prefixes))
    prefix_index = ''.join(Counter(pos).most_common(1)[0][0] for pos in transposed_prefixes)
    for read_name, read_info in relocated_reads.items():
        read_info['seq1_index'] = prefix_index

    return relocated_reads, tmp_fuzzy_dup

def chunk_groups(groups, nthr):
    groups = sorted(groups, key=len, reverse=True)
    chunks = [[] for _ in range(nthr)]
    
    for i, group in enumerate(groups):
        chunks[i % nthr].append(group)
    
    return chunks
        
def process_mm_fuzzy_chunk(chunk, umi_diff, frag_ratio, um_index, um_tree, seq_sim, 
                            prefix_length, unique_coordinates, max_distance):
    tmp_fuzzy_kept = {}
    tmp_fuzzy_dup = {}

    for group in chunk:
        kept, dup = process_mm_fuzzy_group(
            group = group, 
            umi_diff = umi_diff, 
            frag_ratio = frag_ratio,
            um_index = um_index,
            um_tree = um_tree,
            seq_sim = seq_sim, 
            prefix_length = prefix_length,
            unique_coordinates = unique_coordinates,
            max_distance = max_distance
            )
        tmp_fuzzy_kept.update(kept)
        for key, value in dup.items():
            if key in tmp_fuzzy_dup:
                tmp_fuzzy_dup[key].extend(value)
            else:
                tmp_fuzzy_dup[key] = value
        
    return tmp_fuzzy_kept, tmp_fuzzy_dup

def multi_fuzzy_matches(groups, umi_diff, frag_ratio, nthr, 
                        seq_sim, um_kept_dict, prefix_length,
                        max_distance):
    
    um_index, um_tree = build_um_index(um_kept_dict, prefix_length)
    
    unique_coordinates = {}
    for um_read in um_kept_dict.values():
        position = um_read['start'] if um_read['strand'] == '+' else um_read['end']
        chrom = um_read['chrom']
        if chrom not in unique_coordinates:
            unique_coordinates[chrom] = (set(), 0)
        if position not in unique_coordinates[chrom][0]:
            unique_coordinates[chrom][0].add(position)
            unique_coordinates[chrom] = (unique_coordinates[chrom][0], unique_coordinates[chrom][1] + 1)

    # Sort the positions for each chromosome
    for chrom in unique_coordinates:
        positions, count = unique_coordinates[chrom]
        unique_coordinates[chrom] = (sorted(positions), count) 
        
    group_chunks = list(chunk_groups(groups = groups, nthr = nthr))

    # Process each chunk
    results = Parallel(n_jobs=nthr)(
        delayed(process_mm_fuzzy_chunk)(
            chunk = chunk, 
            umi_diff = umi_diff, 
            frag_ratio = frag_ratio,
            um_index = um_index,
            um_tree = um_tree,
            seq_sim = seq_sim, 
            prefix_length = prefix_length,
            unique_coordinates = unique_coordinates,
            max_distance = max_distance
            ) for chunk in group_chunks
        )

    # Combine fuzzy matches
    multi_fuzzy_kept = {}
    multi_fuzzy_dup = {}
    
    for kept, dup in results:
        multi_fuzzy_kept.update(kept)
        for key, value in dup.items():
            if key in multi_fuzzy_dup:
                multi_fuzzy_dup[key].extend(value)
            else:
                multi_fuzzy_dup[key] = value

    return multi_fuzzy_kept, multi_fuzzy_dup

def verify_mm_positions_final_pass(input_dict, max_distance):
    seq1_indices = list(set(read_info['seq1_index'] for read_info in input_dict.values()))
    tree = pybktree.BKTree(levenshtein_wrapper, seq1_indices)

    similar_index_groups = defaultdict(set)
    for index in seq1_indices:
        similar_indexes = [match[1] for match in tree.find(index, max_distance)]
        similar_index_groups[index].update(similar_indexes)

    grouped_reads = defaultdict(list)
    visited = set()
    for read_name, read_info in input_dict.items():
        if read_name in visited:
            continue
        read_index = read_info['seq1_index']
        for group_index, group_members in similar_index_groups.items():
            if read_index in group_members:
                grouped_reads[group_index].append(read_info)
                visited.add(read_name)
                break
            
    for group_index, group_reads in grouped_reads.items():
        longest_read = max(group_reads, key=lambda x: (len(x['seq1']), x['mean_qual']))
        longest_len = len(longest_read['seq1'])
        longest_position = (longest_read['start'], longest_read['end'], longest_read['strand'])

        for read_info in group_reads:
            if (read_info['start'], read_info['end'], read_info['strand']) != longest_position:
                input_dict[read_info['read_name']]['start'] = longest_read['start']
                input_dict[read_info['read_name']]['end'] = longest_read['end']
                input_dict[read_info['read_name']]['strand'] = longest_read['strand']
                input_dict[read_info['read_name']]['multi'] = 'True - grouped'

    return input_dict