from intmap.crop import (
    compile_patterns,
    process_reads_parallel,
    init_params
    )
import regex
import math
import gzip
import os
import numpy as np

TEST_DATA_DIR = "tests/data/crop"
os.makedirs(TEST_DATA_DIR, exist_ok = True)

def create_test_fastq():
    ltr5 = "GCGCGC"
    ltr3 = "TGTGTG"
    ltr3_prime = "TGTGTA"
    linker5 = "ATATAT" 
    linker3 = "ACACTC"
    
    # R1 with LTR
    r1_reads = [
        "@read1",
        f"{ltr5}{ltr3}AAAAGTCTAGCTAG",
        "+",
        "IIIIIIIIIIIIIII",
        "@read2",
        f"{ltr5}{ltr3_prime}AAAAGTCTAGCTAG",
        "+", 
        "IIIIIIIIIIIIIIIII"
    ]
    
    # R2 with linker
    r2_reads = [
        "@read1",
        f"{linker5}{linker3}CCCCGTCTAGCTAG",
        "+",
        "IIIIIIIIIIIIIII",
        "@read2", 
        f"{linker5}{linker3}CCCCGTCTAGCTAG",
        "+",
        "IIIIIIIIIIIIIII"
    ]
    
    with gzip.open(f'{TEST_DATA_DIR}/test_R1.fq.gz', 'wt') as f1, gzip.open(f'{TEST_DATA_DIR}/test_R2.fq.gz', 'wt') as f2:
        f1.write('\n'.join(r1_reads))
        f2.write('\n'.join(r2_reads))

def test_crop_perfect_and_mismatch():
    # Set up test parameters
    class Args:
        nm = 'test'
        r1 = f'{TEST_DATA_DIR}/test_R1.fq.gz'
        r2 = f'{TEST_DATA_DIR}/test_R2.fq.gz'
        min_qual = 20
        ltr3 = "TGTGTG"
        linker3 = "ACACTC"
        ltr5 = "GCGCGC"
        linker5 = "ATATAT"
        ltr3_error_rate = 0.2
        linker3_error_rate = 0.2
        ltr5_error_rate = 0.2
        linker5_error_rate = 0.2
        min_frag_len = 10
        read_length = 100
        ltr_umi_len = 0
        ltr_umi_offset = 0
        linker_umi_len = 0
        linker_umi_offset = 0
        ltr_umi_pattern = None
        linker_umi_pattern = None
        c = None
        nthr = 1
        crop_chunk_size = 10

    create_test_fastq()
    
    args = Args()
    params = init_params(args)
    patterns = compile_patterns(args.ltr3, args.linker3, args.ltr5, args.linker5,
                                args.ltr3_error_rate, args.linker3_error_rate,
                                args.ltr5_error_rate, args.linker5_error_rate)
    is_zipped = True
    out_nm = 'out_nm'
    with gzip.open(args.r1) as f1, gzip.open(args.r2) as f2:
        chunk1 = [line.decode() for line in f1.readlines()]
        chunk2 = [line.decode() for line in f2.readlines()]
        cropped_reads1, cropped_reads2 = process_reads_parallel(chunk1, chunk2, patterns, params, is_zipped, out_nm)
    
    # Verify results
    assert len(cropped_reads1) == 2, "Should find both reads"
    assert cropped_reads1[0]['seq'] == "AAAAGTCTAGCTAG", "Perfect match sequence incorrect"
    assert cropped_reads1[1]['seq'] == "AAAAGTCTAGCTAG", "Mismatch sequence incorrect"