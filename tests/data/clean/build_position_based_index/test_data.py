import numpy as np
import faiss

def get_expected_data():
    n_reads = 3
    n_kmers = 6  # len_diff + 1
    
    # Create expected vectors
    vectors = np.zeros((n_reads, n_kmers), dtype=np.float32)
    vectors = vectors / np.linalg.norm(vectors)
    
    # Create position array
    positions = np.array([
        ('chr1', 100, '+'),
        ('chr1', 100, '+'),
        ('chr2', 200, '-')
    ], dtype=[('chrom', 'U20'), ('pos', np.int32), ('strand', 'U1')])
    
    # Create read names array
    read_names = np.array(['read1', 'read2', 'read3'])
    
    # Build expected FAISS index
    index = faiss.IndexFlatL2(n_kmers)
    index.add(vectors)
    
    return index, positions, read_names