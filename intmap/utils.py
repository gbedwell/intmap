import os
import gzip
import io
import regex
import numpy as np
import faiss
import glob
import subprocess
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict, Counter
import math
from joblib import Parallel, delayed
import pandas as pd
import statsmodels.api as sm
import patsy
import ruptures as rpt
import pyranges as pr
import mmh3
import parasail

__all__ = [
    'zipped',
    'ensure_open',
    'open_file', 
    'natural_key',
    'check_genome_compression',
    'fetch_sequence',
    'is_misprimed',
    'set_faiss_threads',
    'apply_minhash',
    'get_nn',
    'group_similar_hashes',
    'extract_grouped_entries',
    'trans_table',
    'revcomp',
    'collapse_group',
    'final_pass_collapse',
    'sonic_abundance',
    'calculate_consensus',
    'check_consensus',
    'define_consensus',
    'diversity_metrics',
    'compare_annotations'
]

trans_table = bytes.maketrans(b"acgtumrwsykvhdbnACGTUMRWSYKVHDBN", 
                              b"TGCAAKYWSRMBDHVNTGCAAKYWSRMBDHVN")

def revcomp(seq):
    return seq.translate(trans_table)[::-1]

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
        raise Exception("Invalid file input or file is closed")

class gzipParallel:
    def __init__(self, filename, mode, nthr = 1):
        self.filename = filename
        self.mode = mode
        self.nthr = nthr
        self.proc = None
        self.raw_out = None
        self.fobj = None

    def __enter__(self):
        if 'r' in self.mode:
            self.proc = subprocess.Popen(
                ['bgzip', '-d', '-c', '-@', str(self.nthr), self.filename],
                stdout = subprocess.PIPE,
                bufsize = 8 * 1024 * 1024
            )
            self.fobj = io.TextIOWrapper(self.proc.stdout, encoding = 'utf-8')
        else:
            self.raw_out = open(self.filename, 'wb', buffering = 8 * 1024 * 1024)
            self.proc = subprocess.Popen(
                ['bgzip', '-c', '-@', str(self.nthr)],
                stdin = subprocess.PIPE,
                stdout = self.raw_out,
                bufsize = 8 * 1024 * 1024
            )
            self.fobj = io.TextIOWrapper(self.proc.stdin, encoding = 'utf-8')
        return self.fobj

    def __exit__(self, *args):
        self.fobj.flush()
        self.fobj.close()
        self.proc.wait()
        if self.raw_out:
            self.raw_out.close()

def open_file(filename, mode, nthr = 1):
    if zipped(filename):
        if nthr > 1:
            return gzipParallel(filename, mode, nthr)
        else:
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
        start = max(((coordinates['end'] - shift) - 50), 0)
        end = min((coordinates['end'] - shift) + 50, chr_len)
    else:
        start = max(((coordinates['start'] + shift) - 50), 0)
        end = min(((coordinates['start'] + shift) + 50), chr_len)
        
    sequence = genome.fetch(chrom, start, end)
    
    if strand == '-':
        sequence = revcomp(sequence)
    return sequence.upper()


def is_misprimed(seq, primer, matrix, seed_len = 6, max_seed_error = 0.2, min_match_frac = 0.7,
                 gap_open = 10, gap_extend = 2):
    def aln_score(query, target, seed_at_5prime):
        result = parasail.sw_trace_striped_16(query, target, gap_open, gap_extend, matrix)
        if result.score <= 0:
            return False

        q_aln = result.traceback.query
        t_aln = result.traceback.ref
        q_aligned_len = sum(1 for c in q_aln if c != '-')
        q_end = result.end_query
        q_start = q_end - q_aligned_len + 1
        qlen = len(query)

        seed_error = 0
        match = 0

        def in_seed(pos):
            return pos < seed_len if seed_at_5prime else pos >= qlen - seed_len

        for i in range(q_start):
            if in_seed(i):
                seed_error += 1

        q_pos = q_start
        for k in range(len(q_aln)):
            qch = q_aln[k]
            tch = t_aln[k]
            if qch == tch:
                match += 1
            elif in_seed(q_pos):
                seed_error += 1
            if qch != '-':
                q_pos += 1

        for i in range(q_end + 1, qlen):
            if in_seed(i):
                seed_error += 1

        return (seed_error / seed_len <= max_seed_error) and (match >= len(query) * min_match_frac)

    return aln_score(primer, seq, seed_at_5prime = False) or \
           aln_score(revcomp(primer), seq, seed_at_5prime = True)


def get_minimizers(string, k, w):
    n_kmers = len(string) - k + 1
    if n_kmers <= 0:
        return []
    if n_kmers < w:
        return [string[i:i + k].encode('utf8') for i in range(n_kmers)]

    kmer_hashes = []
    for i in range(n_kmers):
        h = mmh3.hash(string[i:i+k].encode('utf8'), seed = 1, signed=False)
        kmer_hashes.append(h)

    minimizer_positions = set()
    for start in range(n_kmers - w + 1):
        min_pos = start
        for j in range(start + 1, start + w):
            if kmer_hashes[j] < kmer_hashes[min_pos]:
                min_pos = j
        minimizer_positions.add(min_pos)

    return [string[pos:pos + k].encode('utf8') for pos in sorted(minimizer_positions)]

def apply_minhash(key, value, string, num_perm, token_size, b_bits, minimizer_win):
    hashvalues = np.full(num_perm, np.iinfo(np.uint32).max, dtype = np.uint32)

    if minimizer_win is not None:
        tokens = get_minimizers(string, token_size, minimizer_win)
    else:
        tokens = [string[i:i + token_size].encode('utf8')
                  for i in range(len(string) - token_size + 1)]

    for token in tokens:
        for perm_idx in range(num_perm):
            h = mmh3.hash(token, seed=perm_idx, signed=False)
            if h < hashvalues[perm_idx]:
                hashvalues[perm_idx] = h

    if b_bits == 1:
        bits = (hashvalues & 1).astype(np.uint8)
        packed = np.packbits(bits, bitorder = 'little')
    else:
        mask = (1 << b_bits) - 1
        truncated = (hashvalues & mask).astype(np.uint32)
        total_bits = num_perm * b_bits
        n_bytes = (total_bits + 7) // 8
        packed = np.zeros(n_bytes, dtype = np.uint8)
        for idx in range(num_perm):
            bit_offset = idx * b_bits
            byte_idx = bit_offset // 8
            bit_idx = bit_offset % 8
            val = int(truncated[idx])
            bits_remaining = b_bits
            val_shift = 0
            while bits_remaining > 0:
                bits_in_this_byte = min(8 - bit_idx, bits_remaining)
                packed[byte_idx] |= ((val >> val_shift) << bit_idx) & 0xFF
                val_shift += bits_in_this_byte
                bits_remaining -= bits_in_this_byte
                byte_idx += 1
                bit_idx = 0

    value['hash'] = packed
    return key, value

def set_faiss_threads(nthr):
    if hasattr(faiss, 'omp_set_num_threads'):
        faiss.omp_set_num_threads(nthr)
        return True
    return False

def get_nn(input_dict, num_perm, nthr, b_bits):
    set_faiss_threads(nthr)
    n = len(input_dict)
    nn_k = min(n, 1000)

    total_bits = num_perm * b_bits
    dim = (total_bits + 7) // 8
    n_bits = dim * 8

    db = np.zeros((n, dim), dtype='uint8')
    for i, (_, value) in enumerate(input_dict.items()):
        db[i] = value['hash']

    if n > 5000:
        n_list = min(int(np.sqrt(n)), 256)
        n_search = math.ceil(n_list * 0.125)
        quantizer = faiss.IndexBinaryFlat(n_bits)
        index = faiss.IndexBinaryIVF(quantizer, n_bits, n_list)
        index.train(db)
        index.add(db)
        index.nprobe = min(n_list, n_search)
    else:
        index = faiss.IndexBinaryFlat(n_bits)
        index.add(db)

    distances, indices = index.search(db, k = nn_k)
    return distances, indices

def group_similar_hashes(hash_distances, hash_indices, nthr, similarity, num_perm, input_dict, b_bits):
    hamming_threshold = (num_perm * b_bits) * (1 - similarity) * 0.5
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
        shape = (n_sequences, n_sequences),
        dtype = np.int8
    )
    graph = graph.maximum(graph.T)
    
    n_components, labels = connected_components(
        csgraph = graph,
        directed = False,
        return_labels = True
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

def collapse_group(chrom_strand_tuple, pos_counts, cluster_win, min_count, 
                    abundant_fc, read_mapping):
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
            nearby_mask = np.abs(valid_positions - abundant_positions[i]) <= cluster_win
            if any(nearby_mask):
                valid_counts = abundant_counts[:i][valid_positions_mask]
                valid_dup_counts = abundant_dup_counts[:i][valid_positions_mask]
                
                max_nearby_count = np.max(valid_counts[nearby_mask])
                max_nearby_dup_count = np.max(valid_dup_counts[nearby_mask])
                
                if max_nearby_count >= abundant_counts[i] * abundant_fc:
                    keep_abundant[i] = False
                    collapsed.add(i)
                elif max_nearby_dup_count >= abundant_dup_counts[i] * abundant_fc:
                    keep_abundant[i] = False
                    collapsed.add(i)
    
    final_abundant_positions = abundant_positions[keep_abundant]
    
    read_updates = {}
    for abundant_pos in final_abundant_positions:
        collapse_mask = np.abs(positions - abundant_pos) <= cluster_win
        for pos in positions[collapse_mask]:
            for read_name in read_mapping[(chrom_strand_tuple[0], pos, chrom_strand_tuple[1])]:
                read_updates[read_name] = abundant_pos
    
    return read_updates

def consolidate_new_dups(kept_dict):
      fragment_data = []
      read_names = []

      for read_name, frag in kept_dict.items():
          fragment_data.append((frag['chrom'], frag['start'], frag['end'], frag['strand'], frag['ltr_umi'], frag['linker_umi']))
          read_names.append(read_name)

      fragment_data = np.array(fragment_data, dtype=[('chrom', 'U25'), ('start', 'i8'), ('end', 'i8'), ('strand', 'U1'), ('ltr_umi', 'U50'), ('linker_umi', 'U50')])
      
      read_names = np.array(read_names)

      _, inverse_indices, counts = np.unique(fragment_data, return_inverse = True, return_counts = True)

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
    
def final_pass_collapse(kept_frags, cluster_win, nthr, min_count, abundant_fc, umi_dedup_fp, low_confidence_fc):
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
                    
    results = Parallel(n_jobs = 1)(
        delayed(collapse_group)(
            chrom_strand_tuple = key, 
            pos_counts = pos, 
            min_count = min_count, 
            abundant_fc = abundant_fc, 
            cluster_win = cluster_win, 
            read_mapping = read_mapping)
        for key, pos in grouped_data.items()
    )
    
    for read_updates in results:
        for read_name, pos in read_updates.items():
            if kept_frags[read_name]['strand'] == '-':
                kept_frags[read_name]['end'] = pos
            else:
                kept_frags[read_name]['start'] = pos
    
    reads_to_remove = consolidate_new_dups(kept_frags)

    read_mapping = defaultdict(list)
    for read_name, frag in kept_frags.items():
        if read_name not in reads_to_remove:
            pos = frag['end'] if frag['strand'] == '-' else frag['start']
            read_mapping[(frag['chrom'], pos, frag['strand'])].append(read_name)

    if umi_dedup_fp:
        all_n_umis = all(
            frag['ltr_umi'] == 'N' and frag['linker_umi'] == 'N'
            for frag in kept_frags.values()
        )
        
        if not all_n_umis:
            pos_dup_map = defaultdict(float)
            for read_name, frag in kept_frags.items():
                pos = frag['end'] if frag['strand'] == '-' else frag['start']
                pos_dup_map[(frag['chrom'], pos, frag['strand'])] += frag.get('count', 0)

            umi_pos_reads = defaultdict(lambda: defaultdict(list))
            for pos_key, read_names in read_mapping.items():
                for read_name in read_names:
                    umi = kept_frags[read_name]['ltr_umi'] + kept_frags[read_name]['linker_umi']
                    umi_pos_reads[umi][pos_key].append(read_name)

            def pos_rank(pos_key):
                reads = [read for read in read_mapping[pos_key] if read not in reads_to_remove]
                if not reads:
                    return (0, 0, 0)
                best = max(reads, key = lambda x: (kept_frags[x].get('map_qual', 0), kept_frags[x].get('mean_qual', 0)))
                return (pos_dup_map[pos_key], kept_frags[best].get('map_qual', 0), kept_frags[best].get('mean_qual', 0))

            for umi, pos_read_map in umi_pos_reads.items():
                if len(pos_read_map) > 1:
                    best_pos = max(pos_read_map, key = pos_rank)
                    for pos_key, read_names in pos_read_map.items():
                        if pos_key != best_pos:
                            reads_to_remove.update(read_names)
        
    removed_frags = {}
    for read_name in reads_to_remove:
        removed_frags[read_name] = kept_frags[read_name]
        del kept_frags[read_name]
    
    if len(kept_frags) > 0:
        post_collapse_counts = defaultdict(float)
        for read_name, frag in kept_frags.items():
            pos = frag['end'] if frag['strand'] == '-' else frag['start']
            post_collapse_counts[(frag['chrom'], pos, frag['strand'])] += frag.get('count', 0)
            
        position_dup_counts = list(post_collapse_counts.values())

        total_position_counts = sum(position_dup_counts)
        average_position_count = total_position_counts / len(position_dup_counts)
        print(f"Average number of observations per position: {average_position_count:.2f}", flush = True)

        lc_threshold = math.floor(average_position_count / low_confidence_fc)
        print(f"Low-confidence threshold: {lc_threshold}", flush = True)
        
        low_confidence_positions = {
            pos_key for pos_key, dup_count in post_collapse_counts.items()
            if dup_count < lc_threshold
        }

        low_confidence_reads = []
        for read_name, frag in kept_frags.items():
            pos = frag['end'] if frag['strand'] == '-' else frag['start']
            position_key = (frag['chrom'], pos, frag['strand'])
            if position_key in low_confidence_positions:
                low_confidence_reads.append(read_name)
        
        if low_confidence_reads:
            for read_name in low_confidence_reads:
                removed_frags[read_name] = kept_frags[read_name]
                del kept_frags[read_name]
        
        if len(reads_to_remove) > 0 or len(low_confidence_reads) > 0:
            print(f"Removing {len(reads_to_remove) + len(low_confidence_reads)} additional reads in the final pass.", flush = True)
    
    return kept_frags, removed_frags
  
def pad_tab(tab, start_at = None, end_at = None):
    tab_keys = [int(k) for k in tab.keys()]
    
    if start_at is None:
        start_at = 0
    if end_at is None:
        end_at = max(tab_keys) - 1 + 10
    
    new_tab_names = list(range(start_at, end_at + 1))
    
    if not all(k in new_tab_names for k in tab_keys):
        raise Exception("Invalid lengths - perhaps some are 0")
    
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
    
    curv = -np.sum(phi[:, np.newaxis] ** 2 * np.exp(log_exptp - log_denom), axis = 0)
    
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
            scale = sum(pearson_resid ** 2) / df_resid
            
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
        raise Exception("Lengths of location_indices and lengths must be equal")
    
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
      
def calculate_consensus(file, consensus_length, base_threshold, rc, cpt_penalty, peak_threshold):
    seqs = []
    with open_file(f'{file}', 'rt') as ff:
        line_num = 0
        get_seq = (
            (lambda line: line.strip()[-consensus_length:].rjust(consensus_length, 'N'))
            if rc else
            (lambda line: line.strip()[:consensus_length].ljust(consensus_length, 'N'))
        )

        for line in ff:
            if line_num % 4 == 1:
                seq = get_seq(line)
                seqs.append(seq)
            line_num += 1

    seqs = np.array([list(seq) for seq in seqs])

    consensus = []
    fracs = []
    for col in seqs.T:
        counts = Counter(col)
        total = sum(counts.values())
        most_common_base, freq = max(counts.items(), key=lambda x: x[1])
        fracs.append(freq / total)
        if freq / total >= base_threshold:
            consensus.append(most_common_base)
        else:
            consensus.append('N')

    consensus = ''.join(consensus)
    fracs = np.array(fracs)

    cpt_algo = rpt.Pelt(model = "l2", min_size = 3, jump = 1).fit(fracs)
    breakpoints = cpt_algo.predict(pen = cpt_penalty)
    starts = [0] + breakpoints[:-1]
    ends = breakpoints

    seg_means = np.array([fracs[start:end].mean() for start, end in zip(starts, ends)])
    peak = seg_means.max()
    keep = seg_means >= (peak_threshold * peak)

    est_bases = []
    for (start, end), k in zip(zip(starts, ends), keep):
        if k:
            est_bases.extend(list(consensus[start:end]))
        else:
            est_bases.extend(['N'] * (end - start))
    est_consensus = ''.join(est_bases)

    if rc:
        consensus = revcomp(consensus)
        est_consensus = revcomp(est_consensus)

    return({'raw': consensus, 'estimated': est_consensus})

      
def define_consensus(r1, r2, consensus_length, sample_size, base_threshold, 
                     out_nm, seqtk_path, processed_directory, long_read, cpt_penalty, 
                     peak_threshold):
  
    if r2 and long_read:
        raise Exception("r2 and long_read cannot both be true.")
      
    print(f'Sampling sequences from input files...\n', flush = True)
    
    file_prefix = f'{processed_directory}/{out_nm}_samp_'

    r1_samp = (
        f'{seqtk_path} sample '
        f'-s100 {r1} {sample_size} > '
        f'{file_prefix}R1.fq'
    )
    
    subprocess.call(r1_samp, shell = True)
    
    if r2:
        r2_samp = (
            f'{seqtk_path} sample '
            f'-s100 {r2} {sample_size} > '
            f'{file_prefix}R2.fq'
        )
        
        subprocess.call(r2_samp, shell = True)

    print(f'Calculating consensus sequences...', flush = True)
    r1_consensus = calculate_consensus(
        file = f'{file_prefix}R1.fq',
        consensus_length = consensus_length, 
        base_threshold = base_threshold,
        rc = False,
        cpt_penalty = cpt_penalty,
        peak_threshold = peak_threshold
    )
    
    if long_read:
        r2_consensus = calculate_consensus(
            file = f'{file_prefix}R1.fq',
            consensus_length = consensus_length, 
            base_threshold = base_threshold,
            rc = True,
            cpt_penalty = cpt_penalty,
            peak_threshold = peak_threshold
        )
    elif r2:
        r2_consensus = calculate_consensus(
            file = f'{file_prefix}R2.fq',
            consensus_length = consensus_length, 
            base_threshold = base_threshold,
            rc = False,
            cpt_penalty = cpt_penalty,
            peak_threshold = peak_threshold
        )
    else:
        r2_consensus = None

    for tmp_file in glob.glob(f'{file_prefix}*.fq'):
        os.remove(tmp_file)
        
    return r1_consensus, r2_consensus
  
def consensus_check_internal(file, target_seq, consensus_threshold, consensus_error_rate,
                             cut_path, nthr, processed_directory, rc):
    if rc:
        consensus_trim = (
            f'{cut_path} '
            f'-e {consensus_error_rate} '
            f'-j {nthr} '
            f"--info-file {os.path.join(processed_directory, 'consensus_info.txt')} "
            f'--discard-untrimmed '
            f'-a \"{target_seq};min_overlap={len(target_seq)}\" '
            f"--quiet -o {os.path.join(processed_directory, 'consensus_crop.fq.gz')} "
            f'{file}'
        )
    else:
        consensus_trim = (
            f'{cut_path} '
            f'-e {consensus_error_rate} '
            f'-j {nthr} '
            f"--info-file {os.path.join(processed_directory, 'consensus_info.txt')} "
            f'--discard-untrimmed '
            f'-g \"{target_seq};min_overlap={len(target_seq)};rightmost\" '
            f"--quiet -o {os.path.join(processed_directory, 'consensus_crop.fq.gz')} "
            f'{file}'
        )

    subprocess.call(consensus_trim, shell = True)
    
    contam_data = pd.read_csv(
        f"{os.path.join(processed_directory, 'consensus_info.txt')}", 
        sep='\t',
        header=None,
        names=[0,1,2,3,4,5,6],
        dtype={0: str, 1: int, 4: str, 5: str, 6: str}, 
        usecols=[0,1,2,3,4,5,6],
        on_bad_lines='skip'
    )
    sample_size = len(contam_data)
    contam_data = contam_data[contam_data.iloc[:, 1] != -1]
    
    if len(contam_data) == 0:
        raise Exception("No sequences matches found. Check input sequence and/or data path.")
    
    matches = contam_data[5].tolist()
    match_counts = contam_data[5].value_counts().to_dict()
    total_match = len(matches)
    
    print(f'{total_match} matches found in {sample_size} reads ({total_match / sample_size * 100:.2f}%)', flush = True)

    match_fractions = {key: value / total_match for key, value in match_counts.items()}
    match_fractions = {key: value for key, value in match_fractions.items() if value >= consensus_threshold}
    match_fractions['Other'] = 1 - sum(match_fractions.values())
    
    if rc:
        adjoining_seqs = contam_data[4].dropna().str[-15:].tolist()
    else:
        adjoining_seqs = contam_data[6].dropna().str[:15].tolist()

    seq_counts =  pd.Series(adjoining_seqs).value_counts().to_dict()
    total_seq = len(adjoining_seqs)
    seq_fractions = {key: value / total_seq for key, value in seq_counts.items()}
    overrep = {key: value for key, value in seq_fractions.items() if value >= consensus_threshold}
    
    for pattern in [
        os.path.join(processed_directory, 'consensus_crop.fq.gz'),
        os.path.join(processed_directory, 'consensus_info.txt')
        ]:
        for tmp_file in glob.glob(pattern):
            os.remove(tmp_file)
    
    return match_fractions, overrep


def check_consensus(r1, r2, ltr, linker, sample_size, consensus_threshold, consensus_error_rate, 
                    cut_path, nthr, processed_directory, long_read, seqtk_path, out_nm):
  
    if r2 and long_read:
        raise Exception("r2 and long_read cannot both be true.")
      
    if long_read:
        linker = revcomp(linker)
    
    print(f'Sampling {sample_size} reads from input files...\n', flush = True)
    
    file_prefix = f'{processed_directory}/{out_nm}_samp_'

    r1_samp = (
        f'{seqtk_path} sample '
        f'-s100 {r1} {sample_size} > '
        f'{file_prefix}R1.fq'
    )
    
    subprocess.call(r1_samp, shell = True)
    
    if r2:
        r2_samp = (
            f'{seqtk_path} sample '
            f'-s100 {r2} {sample_size} > '
            f'{file_prefix}R2.fq'
        )
    
        subprocess.call(r2_samp, shell = True)
        
    print('Checking integrant-end consensus...', flush = True)
    ltr_matches, ltr_overrep = consensus_check_internal(
        file = f'{file_prefix}R1.fq', target_seq = ltr, consensus_threshold = consensus_threshold,
        consensus_error_rate = consensus_error_rate, cut_path = cut_path, 
        nthr = nthr, processed_directory = processed_directory, rc = False
    )
    
    if long_read:
        print('\n')
        print('Checking linker-end consensus...', flush = True)
        linker_matches, linker_overrep = consensus_check_internal(
            file = f'{file_prefix}R1.fq', target_seq = linker, consensus_threshold = consensus_threshold, 
            consensus_error_rate = consensus_error_rate, cut_path = cut_path, 
            nthr = nthr, processed_directory = processed_directory, rc = True
        )
    elif r2:
        print('\n')
        print('Checking linker-end consensus...', flush = True)
        linker_matches, linker_overrep = consensus_check_internal(
            file = f'{file_prefix}R2.fq', target_seq = linker, consensus_threshold = consensus_threshold, 
            consensus_error_rate = consensus_error_rate, cut_path = cut_path, 
            nthr = nthr, processed_directory = processed_directory, rc = False
        )
    else:
        linker_matches, linker_overrep = None, None
    
    print('\n')
    print('The following integrant-end matches have been found:', flush = True)
    for seq, fraction in ltr_matches.items():
        print(f'{seq}\t{fraction:.4f}', flush = True)

    print('\n')
    print('The following adjoining sequences are overrepresented in integrant-end reads:', flush = True)
    for seq, fraction in ltr_overrep.items():
        print(f'{seq}\t{fraction:.4f}', flush = True)

    if linker_matches and linker_overrep:
        print('\n')
        print('The following linker-end matches have been found:', flush = True)
        for seq, fraction in linker_matches.items():
            print(f'{seq}\t{fraction:.4f}', flush = True)
        
        print('\n')
        print('The following adjoining sequences are overrepresented in linker-end reads:', flush = True)
        for seq, fraction in linker_overrep.items():
            print(f'{seq}\t{fraction:.4f}', flush = True)
            
    for tmp_file in glob.glob(f'{file_prefix}*.fq'):
        os.remove(tmp_file)
            
def diversity_metrics(site_counts):
    arr = np.asarray(site_counts)
    total = arr.sum()
    if total == 0:
        return {'shannon': 0, 'pielou': 0, 'simpson': 0}
    p = arr[arr > 0] / total
    shannon = -np.sum(p * np.log(p))
    pielou = shannon / np.log(p.size) if p.size > 1 else 0
    simpson = 1 - np.sum(p ** 2)
    return {
        'shannon': shannon,
        'pielou': pielou,
        'simpson': simpson
    }

def compare_annotations(sites, features, out_nm, processed_directory, mle, peaks = False):    
    anno_nm = os.path.splitext(os.path.basename(features))[0]
    feats = pr.read_bed(features)
    nearest = sites.k_nearest(feats, k = 1, strandedness = False, apply_strand_suffix = False)
    overlap = pr.count_overlaps({'Overlap': feats}, sites, strandedness = False)
    total = len(sites) if not mle else sum(sites.Count)
    
    if not mle:
        n_ov = np.count_nonzero(overlap.Overlap > 0)
    else:
        n_ov = overlap.Count[overlap.Overlap > 0].sum()
    print(f"Fraction {'IS' if not peaks else 'peak'} overlap with {anno_nm}: {round(n_ov / total, 4)}")
    
    nearest = nearest.as_df()
    site_cols = [c for c in ['Chromosome', 'Start', 'End', 'Name', 'Count', 'Strand'] if c in nearest.columns]
    feat_cols = [c for c in ['Name_b', 'Start_b', 'End_b', 'Strand_b', 'Distance'] if c in nearest.columns]
    nearest = nearest[site_cols + feat_cols]
    nearest = nearest.rename(
        columns = {
            'Chromosome': 'chrom', 'Start': 'site_start', 'End': 'site_end', 'Name': 'site_name', 'Count': 'count',
            'Strand': 'site_strand', 'Name_b': 'feat_name', 'Start_b': 'feat_start', 'End_b': 'feat_end',
            'Strand_b': 'feat_strand', 'Distance': 'distance'
        }
    )
    
    overlap.to_bed(
        path = os.path.join(processed_directory, f'{out_nm}_{anno_nm}_overlap.bed.gz'), 
        keep = True, 
        compression = 'gzip'
    )
    
    nearest.to_csv(
        path_or_buf = os.path.join(processed_directory, f'{out_nm}_{anno_nm}_nearest.bed.gz'), 
        sep = '\t', 
        header = True, 
        compression = 'gzip'
    )