from intmap.utils import get_nn, apply_minhash, group_similar_hashes
from intmap.utils import collapse_group
from intmap.utils import final_pass_collapse
import random
import numpy as np
from datasketch import MinHash
import math
from .conftest import load_test_data

def generate_test_data(n_sequences = 100, seq_length = 100, num_perm = 128):
    test_dict = {}
    bases = ['A', 'T', 'G', 'C']
    
    # Create groups of similar sequences
    for i in range(0, n_sequences, 3):
        # Generate base sequence
        base_seq = ''.join(np.random.choice(bases, seq_length))
        
        # Create original and two similar variants
        sequences = [
            base_seq,  # Original
            base_seq[:90] + ''.join(np.random.choice(bases, 10)),  # 90% similar
            base_seq[:80] + ''.join(np.random.choice(bases, 20))   # 80% similar
            ]
        
        # Add sequences to test dict
        for j, seq in enumerate(sequences[:min(3, n_sequences-i)]):
            key = f'read_{i+j}'
            value = {'seq1': seq, 'seq2': ''.join(np.random.choice(bases, seq_length))}
            key, value = apply_minhash(key, value, seq, num_perm, token_size=4)
            test_dict[key] = value
            
    return test_dict

def test_get_nn():
    len_diff = 5
    input_dict = generate_test_data(n_sequences = 100, seq_length = 100, num_perm = 128)
    distances, indices = get_nn(input_dict, num_perm = 128, nthr = 1, len_diff = len_diff)
    
    max_rows = len(input_dict) * 1 # Number of iterations in get_nn()
    
    assert distances.shape[0] <= max_rows  # Cannot exceed max possible rows
    assert indices.shape == distances.shape  # Shapes must match
    assert distances.shape[1] <= len(input_dict)  # Each group can have at most len(input_dict) neighbors
    
def test_group_similar_hashes():
    n_sequences = 100
    input_dict = generate_test_data(
        n_sequences = n_sequences, 
        seq_length = 100, 
        num_perm = 128
        )
    
    distances, indices = get_nn(
        input_dict = input_dict, 
        num_perm = 128, 
        nthr = 1, 
        len_diff = 5
        )
    
    unique_entries, grouped_entries = group_similar_hashes(
        hash_distances = distances,
        hash_indices = indices,
        nthr = 1,
        similarity = 0.85,
        num_perm = 128,
        input_dict = input_dict
        )
    
    # Verify output types
    assert isinstance(unique_entries, list)
    assert isinstance(grouped_entries, list)
    assert all(isinstance(x, (int, np.integer)) for x in unique_entries)
    assert all(isinstance(group, list) for group in grouped_entries)
    
    # Verify no overlap between unique and grouped entries
    all_grouped = set([idx for group in grouped_entries for idx in group])
    assert not set(unique_entries).intersection(all_grouped)
    
    # Verify group properties
    for group in grouped_entries:
        assert len(group) > 1  # Groups must contain multiple entries
        assert len(set(group)) == len(group)  # No duplicates
    
    # Verify total number of sequences is preserved
    total_sequences = len(unique_entries) + sum(len(group) for group in grouped_entries)
    assert total_sequences == n_sequences
    
def test_collapse_group():
    # Create test data matching the new function signature
    pos_counts = [
        (100, 80, 0),  # (position, count) - Abundant
        (102, 5, 0),   # Will collapse to 100
        (105, 3, 0),   # Will collapse to 100
        (150, 75, 0),  # Abundant
        (153, 4, 0),   # Will collapse to 150
        (200, 2, 0),   # Below min_count threshold
    ]
    
    # Create read mapping dictionary
    read_mapping = {
        ("chr1", 100, "+"): ["read1"],
        ("chr1", 102, "+"): ["read2"],
        ("chr1", 105, "+"): ["read3"],
        ("chr1", 150, "+"): ["read4"],
        ("chr1", 153, "+"): ["read5"],
        ("chr1", 200, "+"): ["read6"]
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
        len_diff=5,
        min_count=10,
        count_fc=2,
        read_mapping=read_mapping
    )
    
    result = {k: int(v) for k, v in result.items()}
    
    assert result == expected