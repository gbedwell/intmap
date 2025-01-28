def get_expected_data():
    kmer_index = {
        hash('AAATG'): {('chr1', 100, '+')},
        hash('ATGCG'): {('chr1', 100, '+')},
        hash('GCGTA'): {('chr2', 200, '-'), ('chr1', 100, '+')},
        hash('GTAGC'): {('chr2', 200, '-'), ('chr1', 100, '+')},
        hash('AGCGT'): {('chr2', 200, '-'), ('chr1', 100, '+')},
        hash('CGTGG'): {('chr2', 200, '-'), ('chr1', 100, '+')},
        hash('TGGCT'): {('chr1', 100, '+')}
    }
    
    position_groups = {
        ('chr1', 100, '+'): [
            {
                'read_name': 'read1',
                'chrom': 'chr1',
                'start': 100,
                'end': 115,
                'strand': '+',
                'seq1': 'AAATGCGTAGCGTGGC'
            },
            {
                'read_name': 'read2',
                'chrom': 'chr1',
                'start': 100,
                'end': 116,
                'strand': '+',
                'seq1': 'AAATGCGTAGCGTGGCT'
            }
        ],
        ('chr2', 200, '-'): [
            {
                'read_name': 'read3',
                'chrom': 'chr2',
                'end': 200,
                'start': 185,
                'strand': '-',
                'seq1': 'GCGTAGCGTGGC'
            }
        ]
    }
    
    return kmer_index, position_groups