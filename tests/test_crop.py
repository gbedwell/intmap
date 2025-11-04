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
import glob

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
        ltr3_alt = "AGAGAG"
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
        ttr = False
        min_ttr_len = 15
        flag_window_size = 12
        flag_min_diff = 4
        flag_search_limit = 0.5
        

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
        process_reads_parallel(
            chunk1, chunk2, patterns, params, is_zipped, out_nm, 
            processed_directory = f'{TEST_DATA_DIR}', chunk_num = 1, ttr = False
            )
        
    r1_tmp = glob.glob(f'{TEST_DATA_DIR}/out_nm_*_R1_tmp.fq.gz')[0]
    r2_tmp = glob.glob(f'{TEST_DATA_DIR}/out_nm_*_R2_tmp.fq.gz')[0]

    with gzip.open(r1_tmp) as f1, gzip.open(r2_tmp) as f2:
        cropped_reads1 = [line.decode() for line in f1.readlines()]
        cropped_reads2 = [line.decode() for line in f2.readlines()]
        
    for tmp_file in glob.glob(f'{TEST_DATA_DIR}/*_tmp.fq.gz'):
        os.remove(tmp_file)
    
    # Verify results
    num_reads = len(cropped_reads1) // 4
    assert num_reads == 2, "Should find both reads"

    # Check sequence lines specifically (every 2nd line)
    sequences = cropped_reads1[1::4]
    assert sequences[0] == "AAAAGTCTAGCTAG\n", "Perfect match sequence incorrect"
    assert sequences[1] == "AAAAGTCTAGCTAG\n", "Mismatch sequence incorrect"