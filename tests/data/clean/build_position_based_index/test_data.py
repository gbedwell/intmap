import numpy as np
import faiss
import hashlib

def stable_hash(s):
    # bits = platform.architecture()[0]
    # n_digits = 16 if bits == '64bit' else 8
    n_digits = 8
    raw_hash = int(hashlib.sha256(s.encode()).hexdigest()[:n_digits], 16)
    return raw_hash

def get_expected_data():
    n_reads = 3
    k = 3
    len_diff = 5
    n_kmers = len_diff + 1
    
    # Create expected vectors
    vectors = np.zeros((n_reads, n_kmers), dtype=np.float32)
    
    # Generate actual k-mers and their hashes for each sequence
    seq1 = "AAATGCGTAGCGTGGC"
    seq2 = "AAATGCGTAGCGTGGCT"
    seq3 = "GCGTAGCGTGGC"
    
    vectors[0] = [stable_hash(seq1[i:i+k]) for i in range(n_kmers)]
    vectors[1] = [stable_hash(seq2[i:i+k]) for i in range(n_kmers)]
    vectors[2] = [stable_hash(seq3[i:i+k]) for i in range(n_kmers)]
    
    # Normalize vectors
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    
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