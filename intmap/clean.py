import math
from rapidfuzz.distance import Levenshtein
from collections import defaultdict, deque
import numpy as np
import random
from itertools import groupby
from operator import itemgetter
import pybktree
import random
import hashlib
import faiss
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from intmap.utils import *
import multiprocessing
from joblib import Parallel, delayed
import time

# Calculate Hamming distance
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return ValueError('Sequences must be of equal length.')
    return sum(x1 != x2 for x1, x2 in zip(seq1, seq2))

# Calculate sequence similarity based on maximum sequence length
# and the Levenshtein distance
def seq_similarity(seq1, seq2):
    return 1 - Levenshtein.normalized_distance(seq1, seq2)

def vectorized_seq_similarity(windows, target):
    # Convert windows to list of strings
    window_strings = [''.join(w[:len(target)]) for w in windows]
    target_string = ''.join(target)
    
    # Process all windows at once using rapidfuzz
    distances = [Levenshtein.normalized_distance(w, target_string) for w in window_strings]
    return 1 - np.array(distances)

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

def verify_sequence_groups(group, seq_sim):
    if len(group) <= 1:
        return [group]
    
    seq1_key = itemgetter('seq1')
    seq2_key = itemgetter('seq2')
    
    # Pre-compute sequences for all entries
    sequences = [(seq1_key(entry), seq2_key(entry)) for entry in group]
    
    compare_seqs = {}
    for i, (seq1, seq2) in enumerate(sequences):
        if seq_similarity(seq1, revcomp(seq2)) < seq_sim:
            compare_seqs[i] = seq1 + seq2
        else:
            compare_seqs[i] = seq1
    
    subgroups = [[group[0]]]
    first_seqs = sequences[0]
    
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

    results = Parallel(n_jobs=nthr)(
        delayed(process_group)(group) for group in groups
        )
    
    multi_fuzzy_kept = {}
    multi_fuzzy_dup = {}
    for kept, dup in results:
        multi_fuzzy_kept.update(kept)
        multi_fuzzy_dup.update(dup)
    
    return multi_fuzzy_kept, multi_fuzzy_dup

def stable_hash(s):
    # bits = platform.architecture()[0]
    # n_digits = 16 if bits == '64bit' else 8
    n_digits = 8
    raw_hash = int(hashlib.sha256(s.encode()).hexdigest()[:n_digits], 16)
    return raw_hash

def generate_mm_kmers(seqs, k, max_k=25):
    for seq in seqs:
        max_k = min((len(seq) - k + 1), max_k)
        yield np.array([stable_hash(seq[i:i+k]) for i in range(0, max_k)])

def batch_kmer_gen(seq_batch, k):
    return list(generate_mm_kmers(seq_batch, k))

### TO-DO: VERIFY group_mm_sequences MORE RIGOROUSLY ###
def group_mm_sequences(group, seq_sim, len_diff, min_frag_len, k, nthr):
    batch_size = min(1000, len(group))
    sorted_group = sorted(group, key=lambda x: len(x['seq1']))
    seqs = [entry['seq1'] for entry in sorted_group]
    n_seqs = len(seqs)
    set_faiss_threads(nthr)
    # Generate k-mers
    # kmer_arrays = list(generate_mm_kmers(seqs, min_k))
    kmer_arrays = []
    for i in range(0, len(seqs), batch_size):
        batch_seqs = seqs[i:min(i + batch_size, len(seqs))]
        batch_kmers = Parallel(n_jobs=nthr)(
            delayed(batch_kmer_gen)([seq], k) for seq in batch_seqs
        )
        kmer_arrays.extend([kmers[0] for kmers in batch_kmers if kmers])
    
    # Build vectorized index
    unique_kmers = np.unique(np.concatenate(kmer_arrays))
    kmer_to_idx = {kmer: idx for idx, kmer in enumerate(unique_kmers)}
    
    # Create fixed-size vectors for FAISS
    vectors = np.zeros((n_seqs, len(unique_kmers)), dtype=np.float32)
    for idx, kmers in enumerate(kmer_arrays):
        for kmer in kmers:
            vectors[idx, kmer_to_idx[kmer]] = 1

    # Normalize vectors
    norms = np.linalg.norm(vectors, axis=1)
    vectors = vectors / norms[:, np.newaxis]
    
    # Build index and compute distances
    min_points_per_cluster = 50
    if n_seqs < min_points_per_cluster:
        index = faiss.IndexFlatL2(len(unique_kmers))
        index.add(vectors)
    else:
        # Dynamically allocate cluster number
        # Use k-mer variance as a metric of diversity
        # For binary data, max_var = p(1-p) = 0.5 * 0.5
        max_variance = 0.25
        # Normalize variance to be between 0 and 1
        norm_variance = np.var(vectors, axis=0).sum() / (len(unique_kmers) * max_variance)
        # Scale cluster number gradually with sequence count
        seq_scale = np.log10(n_seqs)
        cluster_factor = int(norm_variance * seq_scale * min_points_per_cluster)
        nlist = min(4096, max(1, cluster_factor))
        quantizer = faiss.IndexFlatL2(len(unique_kmers))
        index = faiss.IndexIVFFlat(quantizer, len(unique_kmers), nlist)
        index.train(vectors)
        index.add(vectors)
        index.nprobe = max(nlist // 2, 2) if nlist > 1 else 1
    
    # Build connectivity graph using batched range search
    # adj_matrix = csr_matrix((n_seqs, n_seqs), dtype=bool)
    rows = []
    cols = []
    
    for i in range(0, n_seqs, batch_size):
        batch_end = min(i + batch_size, n_seqs)
        batch_vectors = vectors[i:batch_end]
        

        k_nn = min(100 * 10 ** (1 + n_seqs // 100000), 100000)
        D, I = index.search(batch_vectors, k_nn)
        
        batch_seqs = np.array(seqs[i:batch_end])
        for j in range(batch_end - i):
            matches = I[j][D[j] < ((1 - seq_sim) * 5)]
            
            if len(matches) > 0:
                match_seqs = np.array([seqs[m] for m in matches])
                similarities = np.array([seq_similarity(batch_seqs[j], match_seq) for match_seq in match_seqs])
                valid_matches = matches[similarities >= seq_sim]
                
                # Add valid matches to graph
                valid_matches = valid_matches[valid_matches > (i + j)]
                if len(valid_matches) > 0:
                    rows.extend([i+j] * len(valid_matches))
                    cols.extend(valid_matches)
                    rows.extend(valid_matches)
                    cols.extend([i+j] * len(valid_matches))
    
    adj_matrix = csr_matrix(
        (np.ones(len(rows), dtype=bool), (rows, cols)), shape=(n_seqs, n_seqs)
    )
    
    # Find connected components
    n_components, labels = connected_components(adj_matrix, directed=False)
    
    # Group by components
    subgroups = []
    for i in range(n_components):
        component_indices = np.where(labels == i)[0]
        subgroups.append([sorted_group[idx] for idx in component_indices])
    
    return subgroups

def build_position_based_index(um_kept_dict, len_diff, nthr, k):
    n_reads = len(um_kept_dict)
    n_kmers = len_diff + 1
    vectors = np.zeros((n_reads, n_kmers), dtype=np.float32)
    positions = []
    read_names = []
    
    for idx, (read_name, read) in enumerate(um_kept_dict.items()):
        seq = read['seq1']
        kmers = [stable_hash(''.join(seq[i:i+k])) for i in range(n_kmers)]
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

### TO-DO: VERIFY compare_to_um MORE RIGOROUSLY ###
def compare_to_um(mm_group, um_index, um_positions, um_read_names, 
                    um_kept_dict, seq_sim, k, len_diff, mm_clone_threshold):
    
    if len(mm_group) < len(um_read_names) * mm_clone_threshold:
        return None
    
    longest_mm_read = max(mm_group, key=lambda x: len(x['seq1']))
    longest_mm_seq = longest_mm_read['seq1']
    longest_mm_len = len(longest_mm_seq)
    
    # Generate k-mers from representative sequence to compare to uniquely mapping k-mers
    n_kmers = len_diff + 1
    mm_vectors = np.zeros((1, n_kmers), dtype=np.float32)
    kmers = [stable_hash(''.join(longest_mm_seq[i:i+k])) for i in range(n_kmers)]
    mm_vectors[0] = kmers / np.linalg.norm(kmers)
    
    longest_mm_seq = np.array(list(longest_mm_read['seq1']))
    
    k_neighbors = 5
    distance_threshold = (1 - seq_sim) * 4
    D, I = um_index.search(mm_vectors, k_neighbors)
    
    # Batch process nearest neighbors
    valid_indices = I[0][D[0] < distance_threshold]
    if len(valid_indices) == 0:
        return None
    
    um_reads_batch = np.array([um_kept_dict[str(name)] for name in um_read_names[valid_indices]])
    um_seqs_batch = np.array([np.array(list(read['seq1'])) for read in um_reads_batch], dtype=object)
    um_positions_batch = um_positions[valid_indices]

    valid_length_mask = np.array([len(seq) >= longest_mm_len for seq in um_seqs_batch])
    if not np.any(valid_length_mask):
        return None

    valid_seqs = um_seqs_batch[valid_length_mask]
    max_um_len = np.max([len(seq) for seq in valid_seqs])

    # Create padded sequence array
    padded_seqs = np.full((len(valid_seqs), max_um_len), '-', dtype=str)
    for i, seq in enumerate(valid_seqs):
        padded_seqs[i, :len(seq)] = seq

    # Create windows array
    n_windows = len(padded_seqs[0]) - longest_mm_len + 1
    all_windows = np.zeros((len(valid_seqs), n_windows, longest_mm_len), dtype=str)
    for i, seq in enumerate(padded_seqs):
        all_windows[i] = np.lib.stride_tricks.sliding_window_view(seq, longest_mm_len)

    # Calculate similarities using vectorized operations
    similarities = vectorized_seq_similarity(
        all_windows.reshape(-1, longest_mm_len), 
        longest_mm_seq
    )
    mean_similarities = similarities.reshape(len(valid_seqs), -1).mean(axis=1)

    mask = mean_similarities >= seq_sim
    if not np.any(mask):
        return None

    best_match_idx = np.argmax(mean_similarities * mask)
    best_position = um_positions_batch[best_match_idx]
    
    relocated_group = []
    for read in mm_group:
        relocated_read = read.copy()
        read_len = int(relocated_read['end'] - relocated_read['start'])
        relocated_read['chrom'] = str(best_position['chrom'])
        relocated_read['start'] = int(best_position['pos']) if best_position['strand'] == '+' else int(best_position['pos'] - read_len)
        relocated_read['end'] = int(best_position['pos'] + read_len) if best_position['strand'] == '+' else int(best_position['pos'])
        relocated_read['strand'] = str(best_position['strand'])
        relocated_read['multi'] = 'True - relocated'
        relocated_group.append(relocated_read)
    return relocated_group

def assign_mm_group(mm_group):
    if len(mm_group) == 1:
        return (mm_group, 1)
    
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
        position_strings = [f"{p[0]}_{p[1] if p[3] == '+' else p[2]}_{p[3]}" for p in valid_positions]
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
            
        return (assigned_group, 0)
    
    return (mm_group, 0)

def batch_compare_to_um(chunk, um_index, um_positions, um_read_names,
                        um_kept_dict, seq_sim, k, len_diff, mm_clone_threshold):
    return [compare_to_um(
        mm_group=group,
        um_index=um_index,
        um_positions=um_positions,
        um_read_names=um_read_names,
        um_kept_dict=um_kept_dict,
        seq_sim=seq_sim,
        k=k,
        len_diff=len_diff,
        mm_clone_threshold=mm_clone_threshold
    ) for group in chunk]
        
def verify_mm_positions(mm_kept_dict, um_kept_dict, seq_sim, nthr, len_diff, 
                        k, min_frag_len, mm_clone_threshold):
    
    um_index, um_positions, um_read_names = build_position_based_index(
        um_kept_dict = um_kept_dict, 
        len_diff = len_diff, 
        nthr = nthr, 
        k = k
        )
    
    mm_reads = list(mm_kept_dict.values())
    subgroups = group_mm_sequences(
        group = mm_reads, 
        seq_sim = seq_sim, 
        len_diff = len_diff, 
        k = k, 
        min_frag_len = min_frag_len,
        nthr = nthr
        )
    print(f'Number of multimapping subgroups: {len(subgroups)}', flush = True)  
    
    sorted_subgroups = sorted(subgroups, key=len, reverse=True)
    subgroup_chunk_size = max(1, len(sorted_subgroups) // nthr)
    subgroup_chunks = [
        sorted_subgroups[i:i + subgroup_chunk_size] for i in range(0, len(sorted_subgroups), subgroup_chunk_size)
        ]
    relocated_results = Parallel(n_jobs=nthr)(
        delayed(batch_compare_to_um)(
            chunk = chunk,
            um_index = um_index, 
            um_positions = um_positions, 
            um_read_names = um_read_names,
            um_kept_dict = um_kept_dict, 
            seq_sim = seq_sim, 
            k = k, 
            len_diff = len_diff,
            mm_clone_threshold = mm_clone_threshold
            ) for chunk in subgroup_chunks
        )
    
    relocated_results = [item for sublist in relocated_results for item in sublist]
    
    relocated_reads = {}
    remaining_groups = []
    
    for i, relocated_group in enumerate(relocated_results):
        if relocated_group:
            for read in relocated_group:
                relocated_reads[read['read_name']] = read
        else:
            remaining_groups.append(subgroups[i])

    print(f'Number of multimapping reads relocated to uniquely mapping positions: {len(relocated_reads)}',
            flush = True)

    if remaining_groups:
        reassigned_results = Parallel(n_jobs=nthr)(
            delayed(assign_mm_group)(group) for group in remaining_groups
            )
        
        total_singles = 0
        total_reassigned = 0
        for group, n_single in reassigned_results:
            total_singles += n_single
            for read in group:
                total_reassigned += 1
                relocated_reads[read['read_name']] = read
    
    print(f'Number of multimapping reads grouped and reassigned: {total_reassigned - total_singles}')
    print(f'Number of ungrouped multimapping reads: {total_singles}')
    
    mm_kept_dict.update(relocated_reads)
    return mm_kept_dict