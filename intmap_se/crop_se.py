import os 
import platform
import glob
import subprocess
import sys
import getopt
import regex
import gzip
import io
import math
from collections import defaultdict
import multiprocessing
from joblib import Parallel, delayed
from intmap_se.utils_se import *
import numpy as np

def check_crop_input(ltr3, linker3, ltr1_primer, ltr5, linker5,
                    contamination, ltr5_error_rate, linker5_error_rate,
                    ltr3_error_rate, linker3_error_rate):
    
    if ltr1_primer is not None:
        ltr1_check = ltr1_primer.replace(' ', '')
        ltr1_char = regex.sub('[ATGCN]', '', ltr1_check.upper())
        if ltr1_char != '':
            raise ValueError('Declared round 1 LTR primer may only contain A,T,G,C, and N nucleotides.')

        if len(ltr1_check) < 10:
            raise ValueError('Declared round 1 LTR primer must be at least 10 nucleotides long.')
        
    ltr5_check = ltr5.replace(' ', '')
    ltr5_char = regex.sub('[ATGC]', '', ltr5_check.upper())
    if ltr5_char != '':
        raise ValueError('Declared 5\' LTR sequence may only contain A,T,G, and C nucleotides.')

    if len(ltr5_check) < 5:
        raise ValueError('Declared 5\' LTR sequence must be at least 5 nucleotides long.')
    
    if ltr5_error_rate > 0.4:
        raise ValueError('5\' LTR sequence error rate cannot be > 0.4.')
    
    ltr3_check = ltr3.replace(' ', '')
    ltr3_char = regex.sub('[ATGC]', '', ltr3_check.upper())
    if ltr3_char != '':
        raise ValueError('Declared 3\' LTR sequence may only contain A,T,G, and C nucleotides.')

    if len(ltr3_check) < 5:
        raise ValueError('Declared 3\' LTR sequence must be at least 5 nucleotides long.')
    
    if ltr3_error_rate > 0.4:
        raise ValueError('3\' LTR sequence error rate cannot be > 0.4.')
        
    link5_check = linker5.replace(' ', '')
    link5_char = regex.sub('[ATGCN]', '', link5_check.upper())
    if link5_char != '':
        raise ValueError('Declared 5\' linker sequence may only contain A,T,G,C and N nucleotides.')

    if len(link5_check) < 5:
        raise ValueError('Declared 5\' linker sequence must be at least 5 nucleotides long.')
    
    if linker5_error_rate > 0.4:
        raise ValueError('5\' linker sequence error rate cannot be > 0.4.')
    
    link3_check = linker3.replace(' ', '')
    link3_char = regex.sub('[ATGCN]', '', link3_check.upper())
    if link3_char != '':
        raise ValueError('Declared 3\' linker sequence may only contain A,T,G,C, and N nucleotides.')

    if len(link3_check) < 5:
        raise ValueError('Declared 3\' linker sequence must be at least 5 nucleotides long.')
    
    if linker3_error_rate > 0.4:
        raise ValueError('3\' linker sequence error rate cannot be > 0.4.')
            
    if contamination is not None:
        for cc in contamination:
            cont = cc.replace(' ', '')
            cont_char = regex.sub('[ATGCN]', '', cont.upper())
            if cont_char != '':
                raise ValueError('Declared contamination(s) may only contain A,T,G,C, and N nucleotides.')
            
            if len(cont) < 10:
                raise ValueError('Declared contamination(s) must be at least 10 nucleotides long.')

def compile_patterns(ltr3, linker3, ltr5, linker5,
                    ltr3_error_rate, linker3_error_rate,
                    ltr5_error_rate, linker5_error_rate):
    ltr3_errors = math.floor(len(ltr3) * ltr3_error_rate)
    ltr5_errors = math.floor(len(ltr5) * ltr5_error_rate)
    linker3_errors = math.floor(len(linker3) * linker3_error_rate)
    linker5_errors = math.floor(len(linker5) * linker5_error_rate)
    
    ltr_rc = revcomp(ltr5 + ltr3)[:15]
    linker_rc = revcomp(linker5 + linker3)[:15]
    
    if len(ltr_rc) > 10:
        ltr_rc_errors = 2
    else:
        ltr_rc_errors = 1
        
    if len(linker_rc) > 10:
        linker_rc_errors = 2
    else:
        linker_rc_errors = 1
    
    return {
        'ltr': {
            'perfect': regex.compile(ltr5 + ltr3),
            'mismatch': regex.compile(f'({ltr5}){{s<={ltr5_errors}}}({ltr3}){{s<={ltr3_errors}}}'),
            'indel': regex.compile(f'({ltr5}){{e<={ltr5_errors}}}({ltr3}){{e<={ltr3_errors}}}')
        },
        'linker': {
            'perfect': regex.compile(linker5 + linker3),
            'mismatch': regex.compile(f'({linker5}){{s<={linker5_errors}}}({linker3}){{s<={linker3_errors}}}'),
            'indel': regex.compile(f'({linker5}){{e<={linker5_errors}}}({linker3}){{e<={linker3_errors}}}')
        },
        'ltr_rc': {
            'perfect': regex.compile(ltr_rc),
            'mismatch': regex.compile(f'({ltr_rc}){{s<={ltr_rc_errors}}}'),
            'indel': regex.compile(f'({ltr_rc}){{e<={ltr_rc_errors}}}')
        },
        'linker_rc': {
            'perfect': regex.compile(linker_rc),
            'mismatch': regex.compile(f'({linker_rc}){{s<={linker_rc_errors}}}'),
            'indel': regex.compile(f'({linker_rc}){{e<={linker_rc_errors}}}')
        },
        'ltr_suffix': {
            'perfect': str(ltr5 + ltr3)[-5:]
        },
        'linker_suffix': {
            'perfect': str(linker5 + linker3)[-5:]
        }
    }

def find_pattern_match(seq, patterns, pattern_type, no_error):
    # seq = seq.decode()
    
    match = patterns[pattern_type]['perfect'].search(seq)
    if match:
        return match
    
    if not no_error:
        match = patterns[pattern_type]['mismatch'].search(seq)
        if match:
            return match
        
        match = patterns[pattern_type]['indel'].search(seq)
        if match:
            return match
    
    return None

# # Check quality score type
# def detect_quality_offset(fastq_file, sample_size = 10000):
#     min_qual = float('inf')

#     with open_file(fastq_file, 'rt') as f:
#         for i, line in enumerate(f):
#             if i % 4 == 3:
#                 min_qual = min(min_qual, min(ord(c) for c in line.strip()))
#             if i >= sample_size * 4:
#                 break

#     return 64 if min_qual >= 64 else 33

# # Vectorized quality score conversion
# def qual_to_array(qual_str, offset):
#     return np.frombuffer(qual_str.encode(), dtype=np.int8) - offset

def init_params(args):    
    contam_patterns = None
    if args.c:
        contam_patterns = [regex.compile(f'({cc}){{e<={math.floor(len(cc) * 0.1)}}}') 
                            for cc in args.c]
        
    return {
        'file1': args.r1,
        'min_qual': 10 ** (args.min_qual / -10),
        'min_frag_len': args.min_frag_len,
        'ltr5_error': args.ltr5_error_rate,
        'linker5_error': args.linker5_error_rate,
        'ltr3_error': args.ltr3_error_rate,
        'linker3_error': args.linker3_error_rate,
        'min_len': args.min_frag_len,
        # 'qual_offset': detect_quality_offset(args.r1),
        'ltr_umi_len': args.ltr_umi_len,
        'ltr_umi_offset': args.ltr_umi_offset,
        'ltr_umi_pattern': args.ltr_umi_pattern,
        'contam_patterns': contam_patterns,
        'ltr5': args.ltr5,
        'ltr3': args.ltr3,
        'linker5': args.linker5,        
        'linker3': args.linker3,
        'nthr': args.nthr,
        'chunk_size': args.crop_chunk_size
    }
    
def extract_umis(seq, pattern_match, umi_len, umi_offset, 
                    umi_pattern = None):
    if umi_len > 0:
        start = pattern_match.start() - umi_offset - umi_len
        end = pattern_match.start() - umi_offset
        
        umi = seq[start:end]
        
        if umi_pattern:
            umi_regex_pattern = regex.sub(r'N+', lambda m: f"[ATGC]{{{len(m.group())}}}", umi_pattern)
            if not regex.search(umi_regex_pattern, umi):
                return None
        return umi
    return 'N'
    
def process_fastq_chunk(chunk_data, params, is_zipped):
    n_reads = len(chunk_data) // 4
    sequences = []
    qualities = []
    headers = []
    
    for i in range(n_reads):
        base_idx = i * 4
        head = chunk_data[base_idx].strip().split()[0]
        seq = chunk_data[base_idx + 1].strip()
        qual = chunk_data[base_idx + 3].strip()

        headers.append(head)
        sequences.append(seq)
        # qualities.append(qual_to_array(qual, params['qual_offset']))
        qualities.append(qual)
    
    return sequences, qualities, headers

def fastq_reader(filename, chunk_size):
    with open_file(filename, 'rt') as f:
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) == chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

def fastq_writer(file, reads):
    if isinstance(file, io.TextIOWrapper):
        for read in reads:
            file.write(f"{read['name']}\n{read['seq']}\n+\n{read['qual']}\n")
    else:
        with open_file(file, 'at') as f:
            for read in reads:
                f.write(f"{read['name']}\n{read['seq']}\n+\n{read['qual']}\n")

def process_reads_parallel(chunk1, patterns, params, is_zipped, 
                            out_nm, processed_directory, chunk_num):
    sequences1, qualities1, headers1 = process_fastq_chunk(chunk1, params, is_zipped)
    
    cropped_reads1 = []
    
    ltr_pattern = params['ltr5'] + params['ltr3']
    linker_pattern = params['linker5'] + params['linker3']
    
    no_error = False
    if (params['ltr5_error'] == 0 and params['linker5_error'] == 0 and
        params['ltr3_error'] == 0 and params['linker3_error'] == 0):
        no_error = True

    for i, (seq1, qual1, head1) in enumerate(zip(sequences1, qualities1, headers1)):
        seq1_str = seq1
        
        # if not patterns['ltr_suffix']['perfect'] in seq1_str:
        #     continue
        
        ltr_pattern_match = find_pattern_match(seq1_str, patterns, 'ltr', no_error)
        
        if ltr_pattern_match:
            start1 = ltr_pattern_match.end()
            linker_in_ltr = find_pattern_match(seq1_str, patterns, 'linker_rc', False)

            if linker_in_ltr:
                end1 = linker_in_ltr.start()
                cropped_seq1 = seq1_str[start1:end1]
                cropped_qual1 = qual1[start1:end1]
            else:
                cropped_seq1 = seq1_str[start1:]
                cropped_qual1 = qual1[start1:]
            
            ltr_umi = extract_umis(seq1_str, 
                                    ltr_pattern_match, 
                                    params['ltr_umi_len'], 
                                    params['ltr_umi_offset'],
                                    params['ltr_umi_pattern'])
            
            if not ltr_umi:
                continue
            
            toss = False
            if params['contam_patterns']:
                for pattern in params['contam_patterns']:
                    if regex.match(pattern, cropped_seq1):
                        toss = True
                        break

            if not toss and len(cropped_seq1) >= params['min_len']:
                ltr_found = ltr_pattern_match.group()
                    
                new_header1 = (f'{head1}\tCO:Z:1:{out_nm}\t'
                                    f'RX:Z:{ltr_umi}\t'
                                    f'OX:Z:{ltr_found}')
                
                cropped_reads1.append({
                    'name': new_header1,
                    'seq': cropped_seq1,
                    # 'qual': ''.join(chr(q + params['qual_offset']) for q in cropped_qual1)
                    'qual': cropped_qual1
                })
        else:
            continue
    
    tmp_file1 = os.path.join(processed_directory,    
                            f'{out_nm}_{os.getpid()}_{chunk_num}_R1_tmp.fq.gz')
    
    with gzip.open(tmp_file1, 'wt') as f1:
        fastq_writer(f1, cropped_reads1)

def concatenate_files(outfile, pattern):
    with gzip.open(outfile, 'wb') as out:
        for temp_file in sorted(glob.glob(pattern)):
            with gzip.open(temp_file, 'rb') as infile:
                shutil.copyfileobj(infile, out)

def concatenate_unix(outfile, pattern):
    cmd = ['cat'] + sorted(glob.glob(pattern))
    with open(outfile, 'wb') as out:
        subprocess.run(cmd, stdout=out)
                
def crop(file1, args, processed_directory, out_nm):
    out_file1 = os.path.join(processed_directory,
                            f'{out_nm}_R1_cropped.fq.gz')
    
    if os.path.exists(out_file1):
        os.remove(out_file1)
    
    params = init_params(args)
    
    is_zipped = zipped(params['file1'])
    
    chunks1 = list(fastq_reader(params['file1'], params['chunk_size']))
    
    patterns = compile_patterns(params['ltr3'], params['linker3'], params['ltr5'], params['linker5'],
                                params['ltr3_error'], params['linker3_error'],
                                params['ltr5_error'], params['linker5_error'])
    
    Parallel(n_jobs=params['nthr'])(
        delayed(process_reads_parallel)(
            chunk1, patterns, params, is_zipped, out_nm,
            processed_directory, chunk_num
        ) for chunk_num, chunk1 in enumerate(chunks1)
    )

    inputs = [
        (out_file1, os.path.join(processed_directory, f"{out_nm}_*_R1_tmp.fq.gz")),
        ]
    
    op_sys = platform.system()
    Parallel(n_jobs=params['nthr'])(
        delayed(concatenate_files if op_sys not in ('Darwin', 'Linux') else concatenate_unix)(*x)
        for x in inputs
    )
        
    for pattern in [
        os.path.join(processed_directory, f"{out_nm}_*_R1_tmp.fq.gz"),
        ]:
        for tmp_file in glob.glob(pattern):
            os.remove(tmp_file)