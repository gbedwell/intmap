import math
from rapidfuzz.distance import Levenshtein
from collections import defaultdict
import numpy as np
import random
from itertools import chain
from joblib import Parallel, delayed
import multiprocessing

# Calculate Hamming distance
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return ValueError('Sequences must be of equal length.')
    return sum(x1 != x2 for x1, x2 in zip(seq1, seq2))

# Calculate sequence similarity based on maximum sequence length
# and the Levenshtein distance
def seq_similarity(seq1, seq2):
    return 1 - Levenshtein.normalized_distance(seq1, seq2)

# Define unique exact matching function
# For each key in input_dict, count the number of identical fragments.
# Identical fragments are based on start and end position and UMIs. 
def unique_exact_matches(input_dict):
    tmp_exact_kept = {}
    tmp_exact_dup = set()  # Store already processed reads to avoid duplicates
    
    for key, entries in input_dict.items():
        # Sort entries by 'start' and 'end' positions
        entries.sort(key=lambda x: (x['start'], x['end']))
        
        i = 0
        while i < len(entries):
            current_entry = entries[i]
            current_start = current_entry['start']
            current_end = current_entry['end']
            current_umi1 = current_entry['ltr_umi']
            current_umi2 = current_entry['linker_umi']
            
            if current_entry['read'] in tmp_exact_dup:
                i += 1
                continue

            # Group all entries with the same start, end, and UMI values
            group = [current_entry]
            
            j = i + 1
            while j < len(entries):
                next_entry = entries[j]
                
                # Compare start, end, and UMI values
                if (next_entry['start'] == current_start and 
                    next_entry['end'] == current_end and
                    next_entry['ltr_umi'] == current_umi1 and 
                    next_entry['linker_umi'] == current_umi2):
                    
                    group.append(next_entry)
                else:
                    break
                j += 1
            
            # Determine the entry with the highest quality in the group
            if len(group) > 1:
                best_entry = max(group, key=lambda x: x['mean_qual'])
                total_count = sum(entry['count'] for entry in group if entry != best_entry)
                best_entry['count'] += total_count
                tmp_exact_kept[key] = tmp_exact_kept.get(key, []) + [best_entry]
                
                # Mark all other entries as duplicates
                for entry in group:
                    if entry != best_entry:
                        tmp_exact_dup.add(entry['read'])
            else:
                tmp_exact_kept[key] = tmp_exact_kept.get(key, []) + group
            
            # Move to the next group
            i = j
    
    return tmp_exact_kept


def unique_fuzzy_matches(input_dict, len_diff, umi_diff, frag_ratio):
    tmp_fuzzy_kept = {}
    tmp_fuzzy_dup = {}
    tmp_dup_set = set()  
    
    for key, entries in input_dict.items():
        # Sort entries by 'start' and 'end' positions
        entries.sort(key=lambda x: (x['start'], x['end']))
        
        i = 0
        while i < len(entries):
            current_entry = entries[i]
            current_start = current_entry['start']
            current_end = current_entry['end']
            current_umi1 = current_entry['ltr_umi']
            current_umi2 = current_entry['linker_umi']
            current_count = current_entry['count']
            
            if current_entry['read'] in tmp_dup_set:
                i += 1
                continue
            
            # Initialize an empty group and populate it with entries that match the fuzzy criteria, 
            # not including current_entry
            group = []
            
            k = i - 1
            while k >= 0:
                past_entry = entries[k]
                past_start = past_entry['start']
                past_end = past_entry['end']
                past_umi1 = past_entry['ltr_umi']
                past_umi2 = past_entry['linker_umi']

                # Check the length threshold for start and end positions
                if (abs((current_end - current_start) - (past_end - past_start)) <= 2 * len_diff and 
                    abs(past_start - current_start) <= len_diff and 
                    abs(past_end - current_end) <= len_diff):
                    
                    # Check the UMI difference using Hamming distance
                    if (hamming_distance(current_umi1, past_umi1) <= umi_diff and
                        hamming_distance(current_umi2, past_umi2) <= umi_diff):
                        group.append(past_entry)
                else:
                    break
                k -= 1
            
            j = i + 1
            while j < len(entries):
                next_entry = entries[j]
                next_start = next_entry['start']
                next_end = next_entry['end']
                next_umi1 = next_entry['ltr_umi']
                next_umi2 = next_entry['linker_umi']

                # Check the length threshold for start and end positions
                if (abs((current_end - current_start) - (next_end - next_start)) <= 2 * len_diff and 
                    abs(next_start - current_start) <= len_diff and 
                    abs(next_end - current_end) <= len_diff):
                    
                    # Check the UMI difference using Hamming distance
                    if (hamming_distance(current_umi1, next_umi1) <= umi_diff and
                        hamming_distance(current_umi2, next_umi2) <= umi_diff):
                        group.append(next_entry)
                else:
                    break
                j += 1
                
            # If group is not empty, find the maximum count in the group (excluding current_entry)
            if len(group) > 0:
                # Get the entry with the maximum count
                max_count_entry = max(group, key=lambda x: (x['count'], x['mean_qual']))
                max_count = max_count_entry['count']                            

                # Check the frag_ratio condition
                if (max_count >= (frag_ratio * current_count) + 1):
                    # If the highest count is greater than the threshold, mark current_entry as a duplicate
                    tmp_dup_set.add(current_entry['read'])

                    # Check if max_count_entry is already in tmp_fuzzy_kept
                    if key in tmp_fuzzy_kept:
                        for kept_entry in tmp_fuzzy_kept[key]:
                            if kept_entry['read'] == max_count_entry['read']:
                                # Update the count of max_count_entry in tmp_fuzzy_kept
                                kept_entry['count'] += current_count
                                break
                    else:
                        # Update the count of max_count_entry in input_dict
                        for entry in input_dict[key]:
                            if entry['read'] == max_count_entry['read']:
                                entry['count'] += current_count
                                break

                    kept_read = max_count_entry['read']
                    dup_read = current_entry['read']
                    if kept_read not in tmp_fuzzy_dup:
                        tmp_fuzzy_dup[kept_read] = [dup_read]
                    else:
                        tmp_fuzzy_dup[kept_read].append(dup_read)

                    if dup_read in tmp_fuzzy_dup:
                        tmp_fuzzy_dup[kept_read] += tmp_fuzzy_dup[dup_read]
                        del tmp_fuzzy_dup[dup_read]
            else:
                # If the group is empty (no matching entries), keep the current entry
                tmp_fuzzy_kept[key] = tmp_fuzzy_kept.get(key, []) + [current_entry]
            
            # Move to the next entry
            i += 1
    
    return tmp_fuzzy_kept, tmp_fuzzy_dup

def multi_exact_matches(input_dict):
    multi_exact_kept = {}

    seq_list = []
    multi_keys = []
    multi_values = []
    for key, value in input_dict.items():
        seq_list.append(value['ltr_umi'] + value['linker_umi'] + 
                        value['seq1'] + value['seq2'])
        multi_keys.append(key)
        multi_values.append(value)
        
    seq_array = np.array(seq_list)
    seq_elements, seq_counts = np.unique(seq_array, return_counts = True)
    multi_unique_set = {val for val, count in zip(seq_elements, seq_counts) if count == 1}
    multi_unique_indices = [index for index, value in enumerate(seq_array)
                            if value in multi_unique_set]
    multi_non_unique_set = {val for val, count in zip(seq_elements, seq_counts) if count > 1}
    multi_non_unique_indices = [index for index, value in enumerate(seq_array) 
                                if value in multi_non_unique_set]

    for i in multi_unique_indices:
        multi_exact_kept.update({multi_keys[i]: multi_values[i]})
        
    grouped_non_unique_indices = {}
    for index, value in enumerate(seq_array):
        if value in multi_non_unique_set:
            value = str(value)
            if value not in grouped_non_unique_indices:
                grouped_non_unique_indices[value] = []
            grouped_non_unique_indices[value].append(index)

    for group in grouped_non_unique_indices:
        tmp_keys = [multi_keys[i] for i in grouped_non_unique_indices[group]]
        tmp_values = [multi_values[i] for i in grouped_non_unique_indices[group]]
        best_entry = max(tmp_values, key=lambda x: x['mean_qual'])
        best_idx = tmp_values.index(best_entry)
        total_count = sum(entry['count'] for entry in tmp_values if entry != best_entry)
        best_entry['count'] += total_count
        multi_exact_kept.update({tmp_keys[best_idx]: best_entry})        

    return multi_exact_kept

def multimap_comparisons(index_list, threshold, greater_than=False):
    graph = defaultdict(set)

    for (i, j), score in index_list:
        if ((greater_than and score > threshold) or 
            (not greater_than and score <= threshold)):
            graph[i].add(j)
            graph[j].add(i)

    visited = set()
    components = []
    
    def iterative_dfs(start_node):
        stack = [start_node]
        component = set()
        
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                component.add(node)
                
                for neighbor in graph[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
                        
        return sorted(component)
    
    for node in graph:
        if node not in visited:
            component = iterative_dfs(node)
            components.append(component)
            
    return components

def chunk_groups(groups, nthr):
    n = max(1, math.ceil(len(groups) / nthr))
    for i in range(0, len(groups), n):
        yield groups[i:i + n]

def process_mm_fuzzy_group(group, input_dict, umi_diff, frag_ratio, sim_threshold):
    len_group = 0
    tmp_dict = {}
    # Create temporary dictionary holding fragment information
    for index in group:
        items = list(input_dict.items())[index]
        tmp_dict.update({items[0]: items[1]})
        len_group += 1

    # Make key and value lists for easy index-based access
    group_values = []
    group_keys = []
    for key, value in tmp_dict.items():
        group_values.append(value)
        group_keys.append(key)
    group_indices = list(range(len(group_values)))
    kept_indices = []

    # Concatenate the sequenced fragments for each read
    seqs = [(entry['seq1'] + entry['seq2']) for entry in group_values]
    
    # Calculate pairwise sequence similarities for each sequence in the group
    n = len(seqs)
    sim_idx = []
    for i in range(n):
        for j in range(i + 1, n):
            obs_seq_sim = seq_similarity(seqs[i], seqs[j])
            sim_idx.append(((i, j), obs_seq_sim))

    # Filter for similarity scores > sim_threshold
    # Add indices for reads that don't satisfy the similarity threshold to kept_indices
    sim_groups = multimap_comparisons(index_list = sim_idx, threshold = sim_threshold, greater_than = True)
    sim_flat = [item for sublist in sim_groups for item in sublist]
    kept_indices.extend(i for i in group_indices if i not in sim_flat)

    # For the remaining groups, check UMI similarity
    if sim_groups:
        for gg in sim_groups:
            tmp_values = [group_values[i] for i in gg]
            tmp_ltr_umi = [entry['ltr_umi'] for entry in tmp_values]
            tmp_linker_umi = [entry['linker_umi'] for entry in tmp_values]

            ltr_umi_idx = []
            linker_umi_idx = []
            n = len(tmp_values)
            for i in range(n):
                for j in range(i + 1, n):
                    ltr_umi_dist = hamming_distance(tmp_ltr_umi[i], tmp_ltr_umi[j])
                    ltr_umi_idx.append(((i, j), ltr_umi_dist))
                    linker_umi_dist = hamming_distance(tmp_linker_umi[i], tmp_linker_umi[j])
                    linker_umi_idx.append(((i, j), linker_umi_dist))
            
            ltr_groups = multimap_comparisons(index_list = ltr_umi_idx, 
                                            threshold = umi_diff, 
                                            greater_than = False)
            linker_groups = multimap_comparisons(index_list = linker_umi_idx, 
                                                threshold = umi_diff, 
                                                greater_than = False)
                        
            # Add indices for reads with unmatched UMIs to kept_indices
            ltr_flat = [item for sublist in ltr_groups for item in sublist]
            linker_flat = [item for sublist in linker_groups for item in sublist]
            ltr_set = set(ltr_flat)
            linker_set = set(linker_flat)
            umi_common = list(ltr_set.intersection(linker_set))
            umi_uncommon = ltr_set.symmetric_difference(linker_set)
            kept_indices.extend(umi_uncommon)

            # Remove unmatched indices from groups
            for sublist in ltr_groups:
                ltr_pruned = [[value for value in sublist if value in umi_common]]

            for sublist in linker_groups:
                linker_pruned = [[value for value in sublist if value in umi_common]]
                
            # Only proceed with groups common to both ltr and linker lists
            # Groups that are not common are discarded entirely
            common_groups = [
                [item for item in ltr_pruned[i] if item in linker_pruned[i]] 
                for i in range(len(ltr_pruned))
                ]

            # For the remaining groups, compare each count value to the maximum count value.
            # If max count >= (frag_ratio * current count) + 1, current read is called
            # a duplicate of max count. Their counts are aggregated together.            
            if len(common_groups[0]) > 0:
                for gg in common_groups:
                    pruned_values = [group_values[i] for i in gg]
                    ref_entry = max(pruned_values, key=lambda x: (x['count'], x['mean_qual']))
                    ref_idx = gg[pruned_values.index(ref_entry)]
                    dup_entries = [entry for entry in pruned_values
                                    if ref_entry['count'] >= ((frag_ratio * entry['count']) + 1)]
                    dup_count = sum(entry['count'] for entry in dup_entries)
                    group_values[ref_idx]['count'] += dup_count
                    
                    ref_key = group_keys[ref_idx]
                    dup_keys = [group_keys[gg[pruned_values.index(entry)]] for entry in dup_entries]
                    tmp_fuzzy_dup = {}
                    if len(dup_keys) > 0:
                        if ref_key in tmp_fuzzy_dup:
                            tmp_fuzzy_dup[ref_key].extend(dup_keys)
                        else:
                            tmp_fuzzy_dup[ref_key] = dup_keys

                    non_dup_entries = [entry for entry in pruned_values 
                                        if ref_entry['count'] < ((frag_ratio * entry['count']) + 1)]
                    non_dup_indices = [gg[pruned_values.index(entry)] for entry in non_dup_entries]
                    kept_indices.extend(non_dup_indices)
            else:
                tmp_fuzzy_dup = {}
    else:
        tmp_fuzzy_dup = {}
    
    tmp_fuzzy_kept = {}
    for idx in kept_indices:
        tmp_fuzzy_kept.update({group_keys[idx]: group_values[idx]})
            
    return tmp_fuzzy_kept, tmp_fuzzy_dup

def process_mm_fuzzy_chunk(chunk, input_dict, umi_diff, frag_ratio, sim_threshold):
    tmp_fuzzy_kept = {}
    tmp_fuzzy_dup = {}

    for group in chunk:
        kept, dup = process_mm_fuzzy_group(
            group, 
            input_dict, 
            umi_diff, 
            frag_ratio, 
            sim_threshold)
        tmp_fuzzy_kept.update(kept)
        for key, value in dup.items():
            if key in tmp_fuzzy_dup:
                tmp_fuzzy_dup[key].extend(value)
            else:
                tmp_fuzzy_dup[key] = value

    return tmp_fuzzy_kept, tmp_fuzzy_dup

# Find fuzzy multimapping duplicates
def multi_fuzzy_matches(groups, input_dict, umi_diff, frag_ratio, sim_threshold, nthr):
    group_chunks = list(chunk_groups(groups, nthr))

    # Process each chunk
    results = Parallel(n_jobs=nthr)(
        delayed(process_mm_fuzzy_chunk)(
            chunk, input_dict, umi_diff, frag_ratio, sim_threshold
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

def build_um_index(um_kept_dict, prefix_length):
    index = defaultdict(list)
    for read_name, um_frag in um_kept_dict.items():
        prefix = um_frag['seq1'][:prefix_length]
        index[prefix].append(um_frag)
    return index

def verify_mm_positions(mm_frag, candidate_um_frags, seq_sim):
    mm_seq1 = mm_frag['seq1']
    mm_strand = mm_frag['strand']

    for um_frag in candidate_um_frags:
        um_seq1 = um_frag['seq1']
        um_strand = um_frag['strand']

        if len(mm_seq1) >= len(um_seq1):
            continue
        
        um_seq1_truncated = um_seq1[:len(mm_seq1)]
        similarity = seq_similarity(mm_seq1, um_seq1_truncated)

        if similarity > seq_sim:
            mm_frag['start'] = um_frag['start'] if um_strand == '+' else (um_frag['end'] - len(mm_seq1))
            mm_frag['end'] = (um_frag['start'] + len(mm_seq1)) if um_strand == '+' else um_frag['end']
            mm_frag['multi'] = 'True - relocated'
            return mm_frag, 1

    return mm_frag, 0

def verify_chunk_positions(mm_chunk, um_index, seq_sim, prefix_length):
    relocated_count = 0
    modified_chunk = []

    for read_name, mm_frag in mm_chunk.items():
        prefix = mm_frag['seq1'][:prefix_length]
        candidate_um_frags = um_index.get(prefix, [])

        modified_mm_frag, relocated = verify_mm_positions(mm_frag, candidate_um_frags, seq_sim)
        modified_chunk.append({read_name: modified_mm_frag})
        relocated_count += relocated

    return modified_chunk, relocated_count

def verify_mm_parallel(mm_kept_dict, um_kept_dict, seq_sim, nthr, chunk_size, prefix_length):
    um_index = build_um_index(um_kept_dict, prefix_length)

    mm_items = list(mm_kept_dict.items())
    mm_chunks = [dict(mm_items[i:i + chunk_size]) for i in range(0, len(mm_items), chunk_size)]

    results = Parallel(n_jobs=nthr)(
        delayed(verify_chunk_positions)(mm_chunk, um_index, seq_sim, prefix_length) for mm_chunk in mm_chunks
        )

    modified_mm_list = [mm for chunk, _ in results for mm in chunk]
    modified_mm_dict = {}
    for mm_frag in modified_mm_list:
        for key, value in mm_frag.items():
            modified_mm_dict[key] = value
    total_relocated = sum(count for _, count in results)

    return modified_mm_dict, total_relocated

def seq_sim_trimmed(seq1, seq2, sim_threshold):
        min_length = min(len(seq1), len(seq2))
        trimmed_seq1 = seq1[:min_length]
        trimmed_seq2 = seq2[:min_length]
        sim_score = seq_similarity(trimmed_seq1, trimmed_seq2)
        hamming_dist = hamming_distance(trimmed_seq1[:10], trimmed_seq2[:10])
        return sim_score if (sim_score > sim_threshold and hamming_dist <= 2) else 0

def process_mm_frag_group(group, input_dict, sim_threshold, unique_dict, win_len):
    len_group = 0
    tmp_dict = {}
    # Create temporary dictionary holding fragment information
    for index in group:
        items = list(input_dict.items())[index]
        tmp_dict.update({items[0]: items[1]})
        len_group += 1

    # Make key and value lists for easy index-based access
    group_values = []
    group_keys = []
    for key, value in tmp_dict.items():
        group_values.append(value)
        group_keys.append(key)
    group_indices = list(range(len(group_values)))
    kept_indices = []

    seqs = [entry['seq1'] for entry in group_values]
    
    # Calculate pairwise sequence similarities for each sequence in the group
    n = len(seqs)
    sim_idx = []
    for i in range(n):
        for j in range(i + 1, n):
            obs_seq_sim = seq_sim_trimmed(seqs[i], seqs[j], sim_threshold)
            sim_idx.append(((i, j), obs_seq_sim))

    sim_groups = multimap_comparisons(index_list = sim_idx, threshold = sim_threshold, greater_than = True)
    
    grouped_indices = set(item for sublist in sim_groups for item in sublist)
    ungrouped_indices = [i for i in group_indices if i not in grouped_indices]
    
    grouped_results = {}
    for group in sim_groups:
        group_keys_subset = [group_keys[idx] for idx in group]
        group_values_subset = [group_values[idx] for idx in group]
        group_reads = {key: value for key, value in zip(group_keys_subset, group_values_subset)}
        
        target_positions = {}
        best_seq_id = None
        for seq_id, read_data in group_reads.items():
            if read_data['multi'] == 'True - relocated':
                best_seq_id = seq_id
                break
            
            current_pos = read_data['start'] if read_data['strand'] == '+' else read_data['end']
            unique_positions = [
                um_read['start'] if um_read['strand'] == '+' else int(um_read['end'])
                for um_read in unique_dict.values()
                ]

            # Count unique reads within +/- win_len bases
            density = sum(1 for pos in unique_positions if (current_pos - win_len) <= pos <= (current_pos + win_len))
            target_positions[seq_id] = density

        # Determine the representative read (highest density)
        
        if best_seq_id is None:
            best_seq_id = max(target_positions, key = target_positions.get)
        best_read = group_reads[best_seq_id]
        best_start, best_end, best_strand = best_read['start'], best_read['end'], best_read['strand']

        # Update all reads in the group to the representative read's coordinates
        for seq_id, read_data in group_reads.items():
            if read_data['multi'] != 'True - relocated':
                read_data['start'] = best_start if best_strand == '+' else (best_end - len(read_data['seq1']))
                read_data['end'] = (best_start + len(read_data['seq1'])) if best_strand == '+' else best_end
                read_data['strand'] = best_strand
                read_data['multi'] = 'True - grouped'

        grouped_results.update(group_reads)

    ungrouped_results = {group_keys[idx]: group_values[idx] for idx in ungrouped_indices}

    combined_results = {**grouped_results, **ungrouped_results}

    return combined_results

def process_mm_frag_chunk(chunk, input_dict, sim_threshold, unique_dict, win_len):
    tmp_mm_kept = {}

    for group in chunk:
        kept = process_mm_frag_group(
            group,
            input_dict,
            sim_threshold,
            unique_dict,
            win_len)
        tmp_mm_kept.update(kept)

    return tmp_mm_kept

def multi_frag_matches(groups, input_dict, sim_threshold, nthr, unique_dict, win_len):
    group_chunks = list(chunk_groups(groups, nthr))

    results = Parallel(n_jobs=nthr)(
        delayed(process_mm_frag_chunk)(
            chunk, input_dict, sim_threshold, unique_dict, win_len
            ) for chunk in group_chunks
        )

    multi_frag_kept = {}
    for kept in results:
        multi_frag_kept.update(kept)

    return multi_frag_kept