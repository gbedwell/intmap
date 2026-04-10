from intmap.utils import get_nn, apply_minhash, group_similar_hashes
from intmap.utils import collapse_group
from intmap.utils import final_pass_collapse
from intmap.utils import consolidate_new_dups
import random
import numpy as np
from datasketch import MinHash
import math
from .conftest import load_test_data
  
def generate_test_data(n_sequences = 100, seq_length = 100, num_perm = 128, token_size = 6, b_bits = 32, win_size = None):
    if n_sequences % 5 != 0:
        raise ValueError("n_sequences must be divisible by 5.")
      
    if seq_length % 2 != 0:
        raise ValueError("seq_length must be divisible by 2.")
      
    if num_perm % 8 != 0:
        raise ValueError("num_perm must be divisible by 8.")
    
    test_dict = {}
    bases = ['A', 'T', 'G', 'C']
    
    def mutate_sequence(seq, mut_rate):
        n_mutations = int(len(seq) * mut_rate)
        seq = list(seq)
        for _ in range(n_mutations):
            mutation_type = np.random.choice(
                ['substitution', 'insertion', 'deletion'], 
                p=[0.9, 0.05, 0.05]
            )
            pos = np.random.randint(0, len(seq))
            
            if mutation_type == 'substitution':
                seq[pos] = np.random.choice([b for b in bases if b != seq[pos]])
            elif mutation_type == 'insertion':
                seq.insert(pos, np.random.choice(bases))
            elif mutation_type == 'deletion' and len(seq) > 1:
                seq.pop(pos)
        
        return ''.join(seq)
    
    np.random.seed(1234)
  
    for i in range(0, n_sequences, 5):
        base_seq = ''.join(np.random.choice(bases, seq_length))
        
        sequences = [
            base_seq,
            mutate_sequence(seq = base_seq, mut_rate = 0.01),
            mutate_sequence(seq = base_seq, mut_rate = 0.025),
            mutate_sequence(seq = base_seq, mut_rate = 0.05),
            mutate_sequence(seq = base_seq, mut_rate = 0.10)
          ]
        
        for j, seq in enumerate(sequences[:min(5, n_sequences - i)]):
            key = f'read_{i + j}'
            value = {'seq1': seq[:seq_length // 2], 'seq2': seq[seq_length - (seq_length // 2):]}
            key, value = apply_minhash(key, value, seq, num_perm, token_size = token_size, b_bits = b_bits, minimizer_win = win_size)
            test_dict[key] = value
            
    return test_dict
    
def test_group_similar_hashes():
    n_sequences = 100
    num_perm = 128
    seq_sim = 0.95
    seq_len1 = 200
    seq_len2 = 5000
    token_size1 = 4
    token_size2 = 8
    hash_sim1 = round(seq_sim ** token_size1 / (2 - seq_sim ** token_size1), 2) - 0.05
    hash_sim2 = round(seq_sim ** token_size2 / (2 - seq_sim ** token_size2), 2) - 0.05
    b_bits = 32
    win_size1 = token_size1
    win_size2 = 2 * token_size2
    total_bits = num_perm * b_bits

    input_dict1 = generate_test_data(
        n_sequences = n_sequences, 
        seq_length = seq_len1, 
        num_perm = num_perm,
        token_size = token_size1,
        b_bits = b_bits,
        win_size = win_size1
      )
    
    distances1, indices1 = get_nn(
        input_dict = input_dict1, 
        num_perm = num_perm, 
        nthr = 1,
        b_bits = b_bits
      )
    
    assert distances1.shape[0] == len(input_dict1)
    assert distances1.shape[1] == len(input_dict1)
    assert indices1.shape == distances1.shape
    
    unique_entries1, grouped_entries1 = group_similar_hashes(
        hash_distances = distances1,
        hash_indices = indices1,
        nthr = 1,
        similarity = hash_sim1,
        num_perm = num_perm,
        input_dict = input_dict1,
        b_bits = b_bits
      )
    
    print(grouped_entries1)
        
    all_grouped1 = set([idx for group in grouped_entries1 for idx in group])
    assert not set(unique_entries1).intersection(all_grouped1)
    
    total_sequences1 = len(unique_entries1) + sum(len(group) for group in grouped_entries1)
    assert total_sequences1 == n_sequences
    
    assert len(grouped_entries1) == 20
    
    first_entry1 = []
    second_entry1 = []
    for grouping in grouped_entries1:
        assert isinstance(grouping, list)
        assert len(grouping) >= 4 and len(grouping) <= 5
        first_entry1.append(grouping[0])
        second_entry1.append(grouping[1])
        
    assert len(first_entry1) == len(list(set(first_entry1)))
    assert len(second_entry1) == len(list(set(second_entry1)))
    assert first_entry1[0] == 0 and first_entry1[19] == 95
    assert second_entry1[0] == 1 and second_entry1[19] == 96
    for idx in range(len(first_entry1) - 1):
        assert first_entry1[idx + 1] - first_entry1[idx] == 5
    for idx in range(len(second_entry1) - 1):
        assert second_entry1[idx + 1] - second_entry1[idx] == 5

    input_dict2 = generate_test_data(
        n_sequences = n_sequences, 
        seq_length = seq_len2, 
        num_perm = num_perm,
        token_size = token_size2,
        b_bits = b_bits,
        win_size = win_size2
      )
    
    distances2, indices2 = get_nn(
        input_dict = input_dict2, 
        num_perm = num_perm, 
        nthr = 1,
        b_bits = b_bits
      )
    
    assert distances2.shape[0] == len(input_dict2)
    assert distances2.shape[1] == len(input_dict2)
    assert indices2.shape == distances2.shape
    
    unique_entries2, grouped_entries2 = group_similar_hashes(
        hash_distances = distances2,
        hash_indices = indices2,
        nthr = 1,
        similarity = hash_sim2,
        num_perm = num_perm,
        input_dict = input_dict2,
        b_bits = b_bits
      )
    
    print(grouped_entries2)
    
    all_grouped2 = set([idx for group in grouped_entries2 for idx in group])
    assert not set(unique_entries2).intersection(all_grouped2)
    
    total_sequences2 = len(unique_entries2) + sum(len(group) for group in grouped_entries2)
    assert total_sequences2 == n_sequences
    
    assert len(grouped_entries2) == 20
    
    first_entry2 = []
    second_entry2 = []
    for grouping in grouped_entries2:
        assert isinstance(grouping, list)
        assert len(grouping) >= 4 and len(grouping) <= 5
        first_entry2.append(grouping[0])
        second_entry2.append(grouping[1])
        
    assert len(first_entry2) == len(list(set(first_entry2)))
    assert len(second_entry2) == len(list(set(second_entry2)))
    assert first_entry2[0] == 0 and first_entry2[19] == 95
    assert second_entry2[0] == 1 and second_entry2[19] == 96
    for idx in range(len(first_entry2) - 1):
        assert first_entry2[idx + 1] - first_entry2[idx] == 5
    for idx in range(len(second_entry2) - 1):
        assert second_entry2[idx + 1] - second_entry2[idx] == 5
    
def test_collapse_group():
    pos_counts = [
        (100, 80, 0),  # (position, count) - Abundant
        (102, 5, 0),   # Will collapse to 100
        (105, 3, 0),   # Will collapse to 100
        (150, 75, 0),  # Abundant
        (153, 4, 0),   # Will collapse to 150
        (200, 2, 0),   # Below min_count threshold
        (201, 1, 0)
    ]
    
    read_mapping = {
        ("chr1", 100, "+"): ["read1"],
        ("chr1", 102, "+"): ["read2"],
        ("chr1", 105, "+"): ["read3"],
        ("chr1", 150, "+"): ["read4"],
        ("chr1", 153, "+"): ["read5"],
        ("chr1", 200, "+"): ["read6"],
        ("chr1", 201, "+"): ["read6"]
    }

    expected = {
        "read1": 100,
        "read2": 100,
        "read3": 100,
        "read4": 150,
        "read5": 150
    }
    
    result = collapse_group(
        chrom_strand_tuple=("chr1", "+"),
        pos_counts=pos_counts,
        cluster_win=5,
        min_count=10,
        abundant_fc=2,
        read_mapping=read_mapping
    )
    
    result = {k: int(v) for k, v in result.items()}
    
    assert result == expected
    
def test_consolidate_new_dups():
    kept_dict = {
        'read1': {'chrom': 'chr1', 'start': 100, 'end': 200, 'strand': '+', 'count': 10, 'map_qual': 30, 'mean_qual': 25, 'ltr_umi': 'N', 'linker_umi': 'N'},
        'read2': {'chrom': 'chr1', 'start': 100, 'end': 200, 'strand': '+', 'count': 5, 'map_qual': 20, 'mean_qual': 20, 'ltr_umi': 'N', 'linker_umi': 'N'},
        'read3': {'chrom': 'chr1', 'start': 300, 'end': 400, 'strand': '-', 'count': 8, 'map_qual': 25, 'mean_qual': 22, 'ltr_umi': 'N', 'linker_umi': 'N'},
        'read4': {'chrom': 'chr1', 'start': 300, 'end': 400, 'strand': '-', 'count': 6, 'map_qual': 15, 'mean_qual': 18, 'ltr_umi': 'N', 'linker_umi': 'N'},
        'read5': {'chrom': 'chr2', 'start': 500, 'end': 600, 'strand': '+', 'count': 12, 'map_qual': 35, 'mean_qual': 30, 'ltr_umi': 'N', 'linker_umi': 'N'}
    }

    reads_to_remove = consolidate_new_dups(kept_dict)

    assert 'read2' in reads_to_remove  # Duplicate of read1
    assert 'read4' in reads_to_remove  # Duplicate of read3
    assert 'read1' not in reads_to_remove  # Best read for the first group
    assert 'read3' not in reads_to_remove  # Best read for the second group
    assert 'read5' not in reads_to_remove  # Unique read

    assert kept_dict['read1']['count'] == 15  # read1 + read2
    assert kept_dict['read3']['count'] == 14  # read3 + read4
    assert kept_dict['read5']['count'] == 12  # No change
