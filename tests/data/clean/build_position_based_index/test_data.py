import numpy as np
import faiss
import hashlib

def stable_hash(s):
    # bits = platform.architecture()[0]
    # n_digits = 16 if bits == '64bit' else 8
    n_digits = 8
    raw_hash = int(hashlib.sha256(s.encode()).hexdigest()[:n_digits], 16)
    return raw_hash

def generate_new_test_data(input_data):
    from intmap.clean import build_position_based_index
    
    index, positions, read_names = build_position_based_index(
        um_kept_dict=input_data,
        k=3,
        len_diff=5,
        nthr=1
    )
    
    return index, positions, read_names

def get_expected_data():
    import numpy as np
    import faiss

    n_kmers = 6

    # Create 2 vectors with dimension 6
    vectors = np.array([
        [np.float32(0.6837285), np.float32(0.46537253), np.float32(0.35371256), np.float32(0.21648666), np.float32(0.2978466), np.float32(0.23505755)],
        [np.float32(0.573622), np.float32(0.45269674), np.float32(0.26834884), np.float32(0.13978006), np.float32(0.21314712), np.float32(0.573622)],
    ], dtype=np.float32)

    # Create position array
    positions = np.array([
        ('chr1', 100, '+'),
        ('chr2', 200, '-'),
    ], dtype=[('chrom', 'U20'), ('pos', np.int32), ('strand', 'U1')])

    # Create read names array
    read_names = np.array([
        'read1',
        'read3',
    ])

    # Build expected FAISS index
    index = faiss.IndexFlatL2(6)
    index.add(vectors)

    return index, positions, read_names

def update_test_data(input_data):
    """Update the get_expected_data function with current implementation output"""
    from intmap.clean import build_position_based_index
    
    # Get the current output
    index, positions, read_names = build_position_based_index(
        um_kept_dict=input_data,
        k=3,
        len_diff=5,
        nthr=1
    )
    
    # Print code that would reproduce this data
    print("def get_expected_data():")
    print("    import numpy as np")
    print("    import faiss")
    print("    ")
    print(f"    n_kmers = {index.d}")
    print("    ")
    
    # Print code to recreate vectors
    print(f"    # Create {index.ntotal} vectors with dimension {index.d}")
    print("    vectors = np.array([")
    for i in range(index.ntotal):
        vec = index.reconstruct(i)
        print(f"        {list(vec)},")
    print("    ], dtype=np.float32)")
    print("    ")
    
    # Print code to recreate positions
    print("    # Create position array")
    print("    positions = np.array([")
    for i in range(len(positions)):
        print(f"        ('{positions[i]['chrom']}', {positions[i]['pos']}, '{positions[i]['strand']}'),")
    print("    ], dtype=[('chrom', 'U20'), ('pos', np.int32), ('strand', 'U1')])")
    print("    ")
    
    # Print code to recreate read names
    print("    # Create read names array")
    print("    read_names = np.array([")
    for name in read_names:
        print(f"        '{name}',")
    print("    ])")
    print("    ")
    
    # Print code to recreate index
    print("    # Build expected FAISS index")
    print(f"    index = faiss.IndexFlatL2({index.d})")
    print("    index.add(vectors)")
    print("    ")
    print("    return index, positions, read_names")