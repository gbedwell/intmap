from intmap.utils import get_nn, apply_minhash, group_similar_hashes
from intmap.utils import collapse_group
from intmap.utils import final_pass_collapse
import random
import numpy as np
from datasketch import MinHash
import math
from .conftest import load_test_data

def generate_test_data(n_sequences = 100, seq_length = 100, num_bits = 128):
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
            key, value = apply_minhash(key, value, seq, num_bits, token_size=4)
            test_dict[key] = value
            
    return test_dict

def test_get_nn():
    len_diff = 5
    input_dict = generate_test_data(n_sequences = 100, seq_length = 100, num_bits = 128)
    distances, indices = get_nn(input_dict, num_bits = 128, nthr = 1, len_diff = len_diff, k = 8)
    
    n_kmers = len_diff + 1
    max_rows = len(input_dict) * n_kmers
    
    assert distances.shape[0] <= max_rows  # Cannot exceed max possible rows
    assert indices.shape == distances.shape  # Shapes must match
    assert distances.shape[1] <= len(input_dict)  # Each group can have at most len(input_dict) neighbors
    
def test_group_similar_hashes():
    n_sequences = 100
    input_dict = generate_test_data(
        n_sequences = n_sequences, 
        seq_length = 100, 
        num_bits = 128
        )
    
    distances, indices = get_nn(
        input_dict = input_dict, 
        num_bits = 128, 
        nthr = 1, 
        len_diff = 5,
        k = 8
        )
    
    unique_entries, grouped_entries = group_similar_hashes(
        hash_distances = distances,
        hash_indices = indices,
        nthr = 1,
        sensitivity = 0.85,
        num_bits = 128,
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
    sites = [
        [100, 80, "read1"],  # Abundant
        [102, 5, "read2"],   # Collapse to read1
        [105, 3, "read3"],   # Collapse to read1
        [150, 75, "read4"],  # Abundant
        [153, 4, "read5"],   # Collapse to read4
        [200, 2, "read6"],   # Too far to collapse
    ]
    
    expected = {
        "read2": 100,
        "read3": 100,
        "read5": 150
    }
    
    result = collapse_group(
        chrom_strand_tuple=("chr1", "+"),
        sites=sites,
        threshold=10,
        len_diff=5
    )
    
    assert result == expected

def test_final_pass_collapse():
    input_data, expected = load_test_data('utils', 'final_pass_collapse')
    result = final_pass_collapse(input_data, len_diff=5, nthr=1)
    assert result == expected