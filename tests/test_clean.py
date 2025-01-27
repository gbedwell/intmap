from intmap.clean import unique_exact_matches
from intmap.clean import unique_fuzzy_matches
from intmap.clean import ranged_groupby
from intmap.clean import cluster_entries_by_umis, build_umi_networks, traverse_umi_networks, bfs
from intmap.clean import multi_exact_matches
from intmap.clean import verify_sequence_groups
from intmap.clean import multi_fuzzy_matches
from intmap.clean import build_kmer_index
from intmap.clean import group_mm_sequences
from intmap.clean import build_position_based_index
from intmap.clean import compare_to_um
from intmap.clean import assign_mm_group
from intmap.clean import verify_mm_positions
from .conftest import load_test_data
import random

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
            'mean_qual': 30
        },
        {
            'read_name': 'read2',
            'ltr_umi': 'AAAT',  # 1 mismatch with read1
            'linker_umi': 'TTTT',
            'count': 5,
            'mean_qual': 25
        },
        {
            'read_name': 'read3',
            'ltr_umi': 'AAAT',  
            'linker_umi': 'TTTA', # 1 mismatch with read2
            'count': 2,
            'mean_qual': 20
        },
        {
            'read_name': 'read4',
            'ltr_umi': 'CGCG',  
            'linker_umi': 'GCGC',
            'count': 1,
            'mean_qual': 20
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
            'mean_qual': random.randint(20, 40)
        })
    
    # Test both scenarios
    small_clusters = cluster_entries_by_umis(small_cluster_input, threshold = 4, frag_ratio = 2)
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
    
    result = [verify_sequence_groups(group, seq_sim=0.95) for group in input_data]
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
        seq_sim=0.9
        )
    assert kept == expected['kept']
    assert dup == expected['dup']
    
def test_build_kmer_index():
    sequences = [
        "AAATGCGTAGCGTGGC",
        "TGCGTAGCGTGGCGAT",
        "GCGTAGCGTGGC"
    ]
    
    kmer_index = build_kmer_index(sequences, k=5, step=2)
    kmer_dict = dict(kmer_index)
    
    expected = {
        'AAATG': {0}, 
        'ATGCG': {0}, 
        'GCGTA': {0, 2}, 
        'GTAGC': {0, 2}, 
        'AGCGT': {0, 2}, 
        'CGTGG': {0, 2}, 
        'TGCGT': {1}, 
        'CGTAG': {1}, 
        'TAGCG': {1}, 
        'GCGTG': {1}, 
        'GTGGC': {1}, 
        'GGCGA': {1}
    }
    
    assert kmer_dict == expected
    
    # Test with different parameters
    kmer_index = build_kmer_index(sequences, k=3, step=1)
    assert "AAA" in kmer_index
    assert "TGC" in kmer_index
    assert len(kmer_index["GCG"]) == 3
    
def test_group_mm_sequences():
    test_group = [
        {'read_name': 'read1', 'seq1': 'AAATGCGTAGCGTGGC', 'count': 1},
        {'read_name': 'read2', 'seq1': 'TGCGTAGCGTGGCGAT', 'count': 1},
        {'read_name': 'read3', 'seq1': 'GCGTAGCGTGGC', 'count': 1},
        {'read_name': 'read4', 'seq1': 'TACGTACGTACGT', 'count': 1}
    ]

    subgroups1 = group_mm_sequences(test_group, seq_sim=0.8, k=5, step=1)

    assert len(subgroups1) == 2
    
    # Verify each group contains the expected sequences
    group_seqs1 = [[read['seq1'] for read in group] for group in subgroups1]
    assert ['GCGTAGCGTGGC', 'AAATGCGTAGCGTGGC', 'TGCGTAGCGTGGCGAT'] in group_seqs1
    assert ['TACGTACGTACGT'] in group_seqs1
    
    subgroups2 = group_mm_sequences(test_group, seq_sim=0.8, k=5, step=2)

    # read2 and read3 share no exact k-mer matches when step=2.
    assert len(subgroups2) == 3
    
    # Verify each group contains the expected sequences
    group_seqs2 = [[read['seq1'] for read in group] for group in subgroups2]
    assert ['GCGTAGCGTGGC', 'AAATGCGTAGCGTGGC'] in group_seqs2
    assert ['TACGTACGTACGT'] in group_seqs2
    assert ['TGCGTAGCGTGGCGAT'] in group_seqs2
    
def test_build_position_based_index():
    input_data, expected = load_test_data('clean', 'build_position_based_index')
    
    # Convert expected kmer_index arrays to sets of tuples
    expected_kmer_index = {
        kmer: {tuple(pos) for pos in positions} 
        for kmer, positions in expected['kmer_index'].items()
    }
    
    kmer_index, position_groups = build_position_based_index(
        um_kept_dict=input_data,
        k=5,
        step=2
    )
    
    # Convert position tuples to strings for comparison
    position_dict = {
        f"{pos[0]}_{pos[1]}_{pos[2]}": [read['read_name'] for read in reads]
        for pos, reads in position_groups.items()
    }
    
    # Convert numeric values in kmer_index to match expected format
    kmer_dict = {
        kmer: {(str(pos[0]), pos[1], str(pos[2])) for pos in positions}
        for kmer, positions in kmer_index.items()
    }
    
    assert kmer_dict == expected_kmer_index
    assert position_dict == expected['position_groups']
    
def test_compare_to_um():
    group = [
        {
            "read_name": "mm_read1",
            "seq1": "GCGTAGCGTGGC",
            "count": 1,
            "chrom": "chr1",
            "start": 300,
            "end": 312,
            "strand": "+"
        }
    ]
    
    um_index = {
        "GCGTA": {("chr2", 200, "-")},
        "GTAGC": {("chr2", 200, "-")},
        "AGCGT": {("chr2", 200, "-")},
        "CGTGG": {("chr2", 200, "-")},
        "TGGCT": {("chr2", 200, "-")}
    }
    
    position_groups = {
        ("chr2", 200, "-"): [
            {
                "read_name": "um_read1",
                "seq1": "GCGTAGCGTGGC",
                "strand": "-",
                "start": 185,
                "end": 200
            }
        ]
    }
    
    expected = [
        {
            "read_name": "mm_read1",
            "seq1": "GCGTAGCGTGGC",
            "count": 1,
            "chrom": "chr2",
            "start": 188,
            "end": 200,
            "strand": "-",
            "multi": "True - relocated"
        }
    ]
    
    result = compare_to_um(
        group=group,
        um_index=um_index,
        position_groups=position_groups,
        seq_sim=0.8,
        k=5,
        step=2
    )
    
    assert result == expected
    
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
    
    result = assign_mm_group(group)
    assert result == expected
    
def test_verify_mm_positions():
    input_data, expected = load_test_data('clean', 'verify_mm_positions')
    
    result = verify_mm_positions(
        mm_kept_dict=input_data['mm_kept_dict'],
        um_kept_dict=input_data['um_kept_dict'],
        seq_sim=0.8,
        nthr=1,
        k = 5,
        step = 2
    )
    
    assert result == expected