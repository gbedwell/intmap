import math
from Levenshtein import distance as lev_dist
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
    max_len = max(len(seq1), len(seq2))
    similarity = 1 - (lev_dist(seq1, seq2) / max_len)
    return similarity

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
            
            # Initialize an empty group and populate it with only entries that match the fuzzy criteria, 
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
            if group:
                # Get the entry with the maximum count
                max_count_entry = max(group, key=lambda x: (x['count'], x['mean_qual']))
                max_count = max_count_entry['count']                            

                # Check the frag_ratio condition
                if (max_count >= (frag_ratio * current_count) + 1):
                    # If the highest count is greater than the threshold, mark current_entry as a duplicate
                    tmp_dup_set.add(current_entry['read'])

                    # Check if max_count_entry is already in tmp_fuzzy_kept
                    found_in_kept = False
                    if key in tmp_fuzzy_kept:
                        for kept_entry in tmp_fuzzy_kept[key]:
                            if kept_entry['read'] == max_count_entry['read']:
                                # Update the count of max_count_entry in tmp_fuzzy_kept
                                kept_entry['count'] += current_count
                                found_in_kept = True
                                break
                    if not found_in_kept:
                        # Update the count of max_count_entry in input_dict
                        for entry in input_dict[key]:
                            if entry['read'] == max_count_entry['read']:
                                entry['count'] += current_count
                                break  # Stop once we've updated the entry

                    kept_read = max_count_entry['read']
                    if kept_read not in tmp_fuzzy_dup:
                        tmp_fuzzy_dup[kept_read] = []
                    tmp_fuzzy_dup[kept_read].append(current_entry['read'])

                else:
                    # Otherwise, keep current_entry
                    tmp_fuzzy_kept[key] = tmp_fuzzy_kept.get(key, []) + [current_entry]
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
        if ((greater_than and score >= threshold) or 
            (not greater_than and score <= threshold)):
            graph[i].add(j)
            graph[j].add(i)

    visited = set()
    components = []
    
    def iterative_dfs(start_node):
        stack = [start_node]
        component = []
        
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                component.append(node)
                
                for neighbor in graph[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
                        
        return component
    
    for node in graph:
        if node not in visited:
            component = iterative_dfs(node)
            components.append(component)
            
    return components

def chunk_groups(groups, nthr):
    n = math.ceil(len(groups) / nthr)
    for i in range(0, len(groups), n):
        yield groups[i:i + n]

def process_mm_fuzzy_group(group, input_dict, umi_diff, frag_ratio, sim_threshold):
    tmp_dict = {}
    # Create temporary dictionary holding fragment information
    for index in group:
        items = list(input_dict.items())[index]
        tmp_dict.update({items[0]: items[1]})

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
            seq_sim = seq_similarity(seqs[i], seqs[j])
            sim_idx.append(((i, j), seq_sim))
            
    # Filter for similarity scores >= sim_threshold
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
                ltr_pruned = [[sublist[i] for i in umi_common]]
            for sublist in linker_groups:
                linker_pruned = [[sublist[i] for i in umi_common]]
                
            # Only proceed with groups common to both ltr and linker lists
            # Groups that are not common are discarded entirely
            common_groups = [
                [item for item in ltr_pruned[i] if item in linker_pruned[i]] 
                for i in range(len(ltr_pruned))
                ]
            
            # For the remaining groups, compare each count value to the maximum count value.
            # If max count >= (frag_ratio * current count) + 1, current read is called
            # a duplicate of max count. Their counts are aggregated together.
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