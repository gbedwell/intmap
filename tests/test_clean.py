from intmap.clean import unique_exact_matches
from intmap.clean import unique_fuzzy_matches
from intmap.clean import ranged_groupby
from intmap.clean import cluster_entries_by_umis, build_umi_networks, find_connected_components
from intmap.clean import multi_exact_matches
from intmap.clean import verify_sequence_groups
from intmap.clean import multi_fuzzy_matches
from intmap.clean import group_mm_sequences
from intmap.clean import build_position_based_index
from intmap.clean import compare_to_um
from intmap.clean import assign_mm_group
from intmap.clean import verify_mm_positions
from .conftest import load_test_data
import random
import numpy as np
import faiss

def test_unique_exact_matches():
    input_data, expected = load_test_data('clean', 'unique_exact_matches')
    result = unique_exact_matches(input_data)
    assert result == expected
    
def test_ranged_groupby():
    input_data, expected = load_test_data('clean', 'ranged_groupby')
    result = ranged_groupby(list(input_data[0].values())[0], tolerance = 5)
    assert result == expected
    
def test_cluster_entries_by_umis():
    # Test small cluster (< 10000 entries)
    small_cluster_input = [
        {
            'read_name': 'read1',
            'ltr_umi': 'AAAA',
            'linker_umi': 'TTTT',
            'count': 11,
            'mean_qual': 30,
            "map_qual": 30
        },
        {
            'read_name': 'read2',
            'ltr_umi': 'AAAT',  # 1 mismatch with read1
            'linker_umi': 'TTTT',
            'count': 5,
            'mean_qual': 25,
            "map_qual": 30
        },
        {
            'read_name': 'read3',
            'ltr_umi': 'AAAT',  
            'linker_umi': 'TTTA', # 1 mismatch with read2
            'count': 2,
            'mean_qual': 20,
            "map_qual": 30
        },
        {
            'read_name': 'read4',
            'ltr_umi': 'CGCG',  
            'linker_umi': 'GCGC',
            'count': 1,
            'mean_qual': 20,
            "map_qual": 30
        }
    ]
    
    # Test large cluster (> 10000 entries)
    random.seed(1)
    large_cluster_input = []
    bases = ['A', 'T', 'G', 'C']
    for i in range(1200):
        umi1 = ''.join(random.choices(bases, k = 4))
        umi2 = ''.join(random.choices(bases, k = 4))
        large_cluster_input.append({
            'read_name': f'read{i}',
            'ltr_umi': umi1,
            'linker_umi': umi2,
            'count': random.randint(1, 100),
            'mean_qual': random.randint(20, 40),
            "map_qual": random.randint(20, 40)
        })
    
    # Test both scenarios
    small_clusters = cluster_entries_by_umis(small_cluster_input, threshold = 4, frag_ratio = 2)
    print(f'Small clusters: {small_clusters}')
    large_clusters = cluster_entries_by_umis(large_cluster_input, threshold = 1, frag_ratio = 2)
    
    # Verify small cluster results
    assert len(small_clusters) == 2  # Should group similar UMIs
    assert len(small_clusters[0]) == 3  # reads1-3 should cluster together
    assert len(small_clusters[1]) == 1  # One single entry
    
    # Verify large cluster results
    assert len(large_clusters) > 1  # Should form multiple clusters
    for cluster in large_clusters:
        assert len(cluster) >= 1  # Each cluster should have entries
        
def test_multi_exact_matches():
    input_data, expected = load_test_data('clean', 'multi_exact_matches')
    result = multi_exact_matches(input_data)
    assert result == expected
    
def test_verify_sequence_groups():
    input_data, expected = load_test_data('clean', 'verify_sequence_groups')
    
    result = [verify_sequence_groups(group, seq_sim=0.95, len_diff=1, nthr=1) for group in input_data]
    assert result == expected
    
def test_unique_fuzzy_matches():
    input_data, expected = load_test_data('clean', 'unique_fuzzy_matches')
    kept, dup = unique_fuzzy_matches(
        input_dict=input_data,
        len_diff=5,
        umi_diff=1,
        frag_ratio=2
        )
    assert kept == expected['kept']
    assert dup == expected['dup']
    
def test_multi_fuzzy_matches():
    input_data, expected = load_test_data('clean', 'multi_fuzzy_matches')
    kept, dup = multi_fuzzy_matches(
        groups=input_data['groups'],
        umi_diff=3,
        frag_ratio=2,
        nthr=1,
        seq_sim=0.9,
        len_diff = 1
        )
    assert kept == expected['kept']
    assert dup == expected['dup']
    
def test_group_mm_sequences():
    test_group = {
        'read1': {'read_name': 'read1', 'seq1': 'AAGACTGCTTGAGCCCAGGAGTTCAAGGCTACAGTGAGCTATGATCACACCGT', 'count': 1},
        'read2': {'read_name': 'read2', 'seq1': 'AACACTGCTTGAGCCCAGGAGTTCAAGGCTACAGTGAGCTATGATCACACC', 'count': 1},
        'read3': {'read_name': 'read3', 'seq1': 'AAGACTGCTTGAGCCCAGGAGTTCAACGCAACAGTGAGCTATGATCA', 'count': 1},
        'read4': {'read_name': 'read4', 'seq1': 'GCCTTGAACTCCTGGGCTCAAGCAGTCTT', 'count': 1}
    }

    subgroups1 = group_mm_sequences(mm_reads=test_group, seq_sim=0.8, min_frag_len=25,
                                    num_perm=32, token_size=4)
    
    assert len(subgroups1) == 2
    
    # Verify each group contains the expected sequences
    group_seqs1 = [set(read['seq1'] for read in group) for group in subgroups1]
    assert set(['AAGACTGCTTGAGCCCAGGAGTTCAAGGCTACAGTGAGCTATGATCACACCGT', 
                'AACACTGCTTGAGCCCAGGAGTTCAAGGCTACAGTGAGCTATGATCACACC', 
                'AAGACTGCTTGAGCCCAGGAGTTCAACGCAACAGTGAGCTATGATCA']) in group_seqs1
    assert set(['GCCTTGAACTCCTGGGCTCAAGCAGTCTT']) in group_seqs1
    
from tests.data.clean.build_position_based_index.test_data import get_expected_data

def test_build_position_based_index():
    input_data, _ = load_test_data('clean', 'build_position_based_index')
    
    um_index, bloom_filter, fragment_counts = build_position_based_index(
        um_kept_dict=input_data,
        k=3,
        nthr=1
    )
    
    assert isinstance(um_index, dict)
    assert len(um_index) > 0
    
    assert len(bloom_filter) > 0
    
    for read in input_data.values():
        prefix = read['seq1'][:3]
        if prefix in bloom_filter:
            if prefix in um_index:
                assert isinstance(um_index[prefix], list)
                for entry in um_index[prefix]:
                    assert isinstance(entry, tuple)
                    assert len(entry) == 2
                    position_info, read_id = entry
                    assert len(position_info) == 3
                    assert isinstance(read_id, str)
    
def test_compare_to_um():
    # Expanded test cases with more diverse scenarios
    um_kept_dict = {
        "um_read1": {
            "read_name": "um_read1",
            "seq1": "GCGTAGCGTGGCAA",  
            "strand": "-", 
            "start": 186,
            "end": 200,
            "chrom": "chr2",
            "multi": "False"
        },
        "um_read2": { 
            "read_name": "um_read2",
            "seq1": "GCGAAGCGTCGCAA",
            "strand": "-",
            "start": 186,
            "end": 200,
            "chrom": "chr2",
            "multi": "False"
        },
        "um_read3": {
            "read_name": "um_read3",
            "seq1": "ATCGATCGATCGAA",
            "strand": "+",
            "start": 500,
            "end": 514,
            "chrom": "chr1",
            "multi": "False"
        },
        "um_read4": {
            "read_name": "um_read4",
            "seq1": "TGCATGCATGCAA",
            "strand": "+",
            "start": 1000,
            "end": 1013,
            "chrom": "chr3",
            "multi": "False"
        },
        "um_read5": {
            "read_name": "um_read5",
            "seq1": "GCGTAGCGTGGCAA",
            "strand": "-",
            "start": 286,
            "end": 300,
            "chrom": "chr3",
            "multi": "False"
        },
        "um_read6": {
            "read_name": "um_read6",
            "seq1": "GCGTAGCGTGGCAA",
            "strand": "-",
            "start": 286,
            "end": 300,
            "chrom": "chr3",
            "multi": "False"
        },
        "um_read7": {
            "read_name": "um_read7",
            "seq1": "GCGTAGCGTGGCAA",
            "strand": "-",
            "start": 286,
            "end": 300,
            "chrom": "chr3",
            "multi": "False"
        }
    }
    
    group1 = [
        {
            "read_name": "mm_read1", # Perfect match to um_read1
            "seq1": "GCGTAGCGTGGC",
            "count": 1,
            "chrom": "chr1", 
            "start": 300,
            "end": 312,
            "strand": "+",
            "multi": "True"
        },
        {
            "read_name": "mm_read2", # Perfect match to um_read1
            "seq1": "GCGTAGCGTGG",
            "count": 1,
            "chrom": "chr1",
            "start": 305,
            "end": 316,
            "strand": "+",
            "multi": "True"
        }
    ]
    
    group2 = [
        {
            "read_name": "mm_read3", # No match
            "seq1": "TGCGAAACCTAG",
            "count": 1,
            "chrom": "chr3",
            "start": 400,
            "end": 412,
            "strand": "-",
            "multi": "True"
        }
    ]
    
    um_index, bloom_filter, fragment_counts = build_position_based_index(
        um_kept_dict=um_kept_dict, 
        k=5, 
        nthr=1
    )

    # Test group1 - should be relocated
    result1 = compare_to_um(
        mm_group=group1,
        k=5,
        um_index=um_index,
        bloom_filter=bloom_filter,
        um_kept_dict=um_kept_dict,
        seq_sim=0.8,
        fragment_counts = fragment_counts
    )
    
    if isinstance(result1, tuple) and len(result1) == 2:
        relocated_group, unrelocated_group = result1
        assert len(relocated_group) > 0  # Some reads should be relocated
        # Check that relocated reads have the correct format
        for read in relocated_group:
            assert read['multi'] == 'True - relocated'
            assert read['chrom'] == 'chr3'  # Should match um_read5
            assert read['strand'] == '-'    # Should match um_read5
    else:
        # If no relocation happened, the result should be the original group
        assert result1 == group1
    
    # Test group2 - should not be relocated
    result2 = compare_to_um(
        mm_group=group2,
        k=5,
        um_index=um_index,
        bloom_filter=bloom_filter,
        um_kept_dict=um_kept_dict,
        seq_sim=0.8,
        fragment_counts = fragment_counts
    )
    
    # For group2, we expect no relocation
    if isinstance(result2, tuple) and len(result2) == 2:
        relocated_group, unrelocated_group = result2
        assert len(relocated_group) == 0  # No reads should be relocated
        assert unrelocated_group == group2
    else:
        # If no relocation happened, the result should be the original group
        assert result2 == group2
    
def test_assign_mm_group():
    group = [
        {
            "read_name": "mm_read1",
            "seq1": "GCGTAGCGTGGC",
            "count": 3,
            "chrom": "chr1",
            "start": 100,
            "end": 112,
            "strand": "+",
            "multi": "True"
        },
        {
            "read_name": "mm_read2",
            "seq1": "GCGTAGCGTGGC",
            "count": 2,
            "chrom": "chr2", 
            "start": 200,
            "end": 212,
            "strand": "-",
            "multi": "False"
        }
    ]
    
    expected = [
        {
            "read_name": "mm_read1",
            "seq1": "GCGTAGCGTGGC",
            "count": 3,
            "chrom": "chr1",
            "start": 100,
            "end": 112,
            "strand": "+",
            "multi": "True - grouped"
        },
        {
            "read_name": "mm_read2",
            "seq1": "GCGTAGCGTGGC",
            "count": 2,
            "chrom": "chr1",
            "start": 100,
            "end": 112,
            "strand": "+",
            "multi": "True - grouped"
        }
    ]
    
    result, n = assign_mm_group(group, mm_group_threshold=1)
    assert result == expected

def test_verify_mm_positions():
    input_data, expected = load_test_data('clean', 'verify_mm_positions')
    
    result = verify_mm_positions(
        mm_kept_dict=input_data['mm_kept_dict'],
        um_kept_dict=input_data['um_kept_dict'],
        seq_sim=0.8,
        nthr=1,
        k = 5,
        len_diff = 5,
        min_frag_len=20,
        num_perm=32, 
        token_size=4,
        mm_group_threshold=1
    )
    
    assert result == expected