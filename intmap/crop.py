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
from joblib import Parallel, delayed
from itertools import tee, count
from intmap.utils import *
import numpy as np

def check_crop_input(ltr3, linker3, ltr1_primer, ltr5, linker5,
                    contamination, ltr5_error_rate, linker5_error_rate,
                    ltr3_error_rate, linker3_error_rate):
    
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

    if len(ltr5_check) < 10:
        raise ValueError('Declared 5\' LTR sequence must be at least 10 nucleotides long.')
    
    if ltr5_error_rate > 0.3:
        raise ValueError('5\' LTR sequence error rate cannot be > 0.3.')
    
    ltr3_check = ltr3.replace(' ', '')
    ltr3_char = regex.sub('[ATGC]', '', ltr3_check.upper())
    if ltr3_char != '':
        raise ValueError('Declared 3\' LTR sequence may only contain A,T,G, and C nucleotides.')

    if len(ltr3_check) < 5:
        raise ValueError('Declared 3\' LTR sequence must be at least 5 nucleotides long.')
    
    if ltr3_error_rate > 0.3:
        raise ValueError('3\' LTR sequence error rate cannot be > 0.3.')
        
    link5_check = linker5.replace(' ', '')
    link5_char = regex.sub('[ATGC]', '', link5_check.upper())
    if link5_char != '':
        raise ValueError('Declared 5\' linker sequence may only contain A,T,G, and C nucleotides.')

    if len(link5_check) < 10:
        raise ValueError('Declared 5\' linker sequence must be at least 10 nucleotides long.')
    
    if linker5_error_rate > 0.3:
        raise ValueError('5\' linker sequence error rate cannot be > 0.3.')
    
    link3_check = linker3.replace(' ', '')
    link3_char = regex.sub('[ATGC]', '', link3_check.upper())
    if link3_char != '':
        raise ValueError('Declared 3\' linker sequence may only contain A,T,G, and C nucleotides.')

    if len(link3_check) < 5:
        raise ValueError('Declared 3\' linker sequence must be at least 5 nucleotides long.')
    
    if linker3_error_rate > 0.3:
        raise ValueError('3\' linker sequence error rate cannot be > 0.3.')
            
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
    
    ltr_rc = revcomp(ltr5 + ltr3)[:12]
    linker_rc = revcomp(linker5 + linker3)[:12]
    rc_errors = 1
    
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
            'mismatch': regex.compile(f'({ltr_rc}){{s<={rc_errors}}}'),
            'indel': regex.compile(f'({ltr_rc}){{e<={rc_errors}}}')
        },
        'linker_rc': {
            'perfect': regex.compile(linker_rc),
            'mismatch': regex.compile(f'({linker_rc}){{s<={rc_errors}}}'),
            'indel': regex.compile(f'({linker_rc}){{e<={rc_errors}}}')
        }
    }

# Add error rates to this to skip looking for non-perfect matches if error = 0.
def find_pattern_match(seq, patterns, pattern_type):
    # seq = seq.decode()
    
    match = patterns[pattern_type]['perfect'].search(seq)
    if match:
        return match
    
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

# Cache frequently used values
def init_params(args):    
    contam_patterns = None
    if args.c:
        contam_patterns = [regex.compile(f'({cc}){{e<={math.floor(len(cc) * 0.1)}}}') 
                            for cc in args.c]
        
    return {
        'file1': args.r1,
        'file2': args.r2,
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
        'linker_umi_len': args.linker_umi_len,
        'linker_umi_offset': args.linker_umi_offset,
        'ltr_umi_pattern': args.ltr_umi_pattern,
        'linker_umi_pattern': args.linker_umi_pattern,
        'contam_patterns': contam_patterns,
        'ltr5': args.ltr5,
        'ltr3': args.ltr3,
        'linker5': args.linker5,        
        'linker3': args.linker3,
        'nthr': args.nthr,
        'chunk_size': args.crop_chunk_size
    }
    
def extract_umis(seq, pattern_match, umi_len, umi_offset, 
                    umi_pattern = None, position = None, no_error = None):
    if umi_len > 0:
        if no_error:
            start = position - umi_offset - umi_len
            end = position - umi_offset
        else:
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
            
def process_reads_parallel(chunk1, chunk2, patterns, params, is_zipped, 
                            out_nm, processed_directory, chunk_num):
    sequences1, qualities1, headers1 = process_fastq_chunk(chunk1, params, is_zipped)
    sequences2, qualities2, headers2 = process_fastq_chunk(chunk2, params, is_zipped)
    
    cropped_reads1 = []
    cropped_reads2 = []
    
    ltr_pattern = params['ltr5'] + params['ltr3']
    linker_pattern = params['linker5'] + params['linker3']
    
    no_error = False
    if (params['ltr5_error'] == 0 and params['linker5_error'] == 0 and
        params['ltr3_error'] == 0 and params['linker3_error'] == 0):

        no_error = True
        
        sequences1 = np.array(sequences1)
        sequences2 = np.array(sequences2)
        qualities1 = np.array(qualities1)
        qualities2 = np.array(qualities2)
        headers1 = np.array(headers1)
        headers2 = np.array(headers2)
        
        ltr_positions = np.char.find(sequences1, ltr_pattern)
        linker_positions = np.char.find(sequences2, linker_pattern)
        
        valid_mask = (ltr_positions >= 0) & (linker_positions >= 0)
        
        read1_pos = ltr_positions[valid_mask]
        read2_pos = linker_positions[valid_mask]
        
        sequences1 = list(sequences1[valid_mask])
        sequences2 = list(sequences2[valid_mask])
        qualities1 = list(qualities1[valid_mask])
        qualities2 = list(qualities2[valid_mask])
        headers1 = list(headers1[valid_mask])
        headers2 = list(headers2[valid_mask])
    else:
        ltr_positions = np.full(len(sequences1), -1)
        linker_positions = np.full(len(sequences2), -1)
    
    for i, ((seq1, qual1, head1, ltr_pos), 
            (seq2, qual2, head2, linker_pos)) in enumerate(zip(zip(sequences1, qualities1, headers1, ltr_positions), 
                                                            zip(sequences2, qualities2, headers2, linker_positions))):
        
        if head1 != head2:
            raise ValueError(f"Read pair mismatch: {head1} != {head2}")
            
        seq1_str = seq1
        seq2_str = seq2
        
        if no_error is False:
            ltr_pattern_match = find_pattern_match(seq1_str, patterns, 'ltr')
            linker_pattern_match = find_pattern_match(seq2_str, patterns, 'linker')
        else:
            ltr_pattern_match, linker_pattern_match = ltr_pattern, linker_pattern
        
        if ltr_pattern_match and linker_pattern_match:
            
            if no_error is False:
                start1 = ltr_pattern_match.end()
                start2 = linker_pattern_match.end()
            else:
                start1 = ltr_pos + len(ltr_pattern)
                start2 = linker_pos + len(linker_pattern)
            
            linker_in_ltr = find_pattern_match(seq1_str, patterns, 'linker_rc')
            ltr_in_linker = find_pattern_match(seq2_str, patterns, 'ltr_rc')

            if linker_in_ltr and ltr_in_linker:
                end1 = linker_in_ltr.start()
                end2 = ltr_in_linker.start()
                cropped_seq1 = seq1_str[start1:end1]
                cropped_seq2 = seq2_str[start2:end2]
                cropped_qual1 = qual1[start1:end1]
                cropped_qual2 = qual2[start2:end2]
            else:
                cropped_seq1 = seq1_str[start1:]
                cropped_seq2 = seq2_str[start2:]
                cropped_qual1 = qual1[start1:]
                cropped_qual2 = qual2[start2:]
            
            ltr_umi = extract_umis(seq1_str, 
                                    ltr_pattern_match, 
                                    params['ltr_umi_len'], 
                                    params['ltr_umi_offset'],
                                    params['ltr_umi_pattern'],
                                    position = ltr_pos if ltr_pos else None,
                                    no_error = no_error
                                    )
            linker_umi = extract_umis(seq2_str, 
                                    linker_pattern_match,
                                    params['linker_umi_len'],
                                    params['linker_umi_offset'],
                                    params['linker_umi_pattern'],
                                    position = linker_pos if linker_pos else None,
                                    no_error = no_error)
            
            if not ltr_umi or not linker_umi:
                continue
            
            toss = False
            if params['contam_patterns']:
                for pattern in params['contam_patterns']:
                    if regex.match(pattern, cropped_seq1):
                        toss = True
                        break

            if not toss and len(cropped_seq1) >= params['min_len'] and len(cropped_seq2) >= params['min_len']:
                ltr_found = (ltr_pattern_match if isinstance(ltr_pattern_match, str) 
                                else ltr_pattern_match.group())
                linker_found = (linker_pattern_match if isinstance(linker_pattern_match, str) 
                                else linker_pattern_match.group())
                    
                new_header1 = (f'{head1}\tCO:Z:1:{out_nm}\t'
                                    f'RX:Z:{ltr_umi}-{linker_umi}\t'
                                    f'OX:Z:{ltr_found}-{linker_found}')
                
                new_header2 = (f'{head1}\tCO:Z:2:{out_nm}\t'
                                    f'RX:Z:{ltr_umi}-{linker_umi}\t'
                                    f'OX:Z:{ltr_found}-{linker_found}')
                
                cropped_reads1.append({
                    'name': new_header1,
                    'seq': cropped_seq1,
                    # 'qual': ''.join(chr(q + params['qual_offset']) for q in cropped_qual1)
                    'qual': cropped_qual1
                })
                
                cropped_reads2.append({
                    'name': new_header2,
                    'seq': cropped_seq2,
                    # 'qual': ''.join(chr(q + params['qual_offset']) for q in cropped_qual2)
                    'qual': cropped_qual2
                })
    
    tmp_file1 = os.path.join(processed_directory,    
                            f'{out_nm}_{os.getpid()}_{chunk_num}_R1_tmp.fq.gz')
    tmp_file2 = os.path.join(processed_directory,
                            f'{out_nm}_{os.getpid()}_{chunk_num}_R2_tmp.fq.gz')
    
    with gzip.open(tmp_file1, 'wt') as f1, gzip.open(tmp_file2, 'wt') as f2:
        fastq_writer(f1, cropped_reads1)
        fastq_writer(f2, cropped_reads2)

def concatenate_files(outfile, pattern):
    with gzip.open(outfile, 'wb') as out:
        for temp_file in sorted(glob.glob(pattern)):
            with gzip.open(temp_file, 'rb') as infile:
                shutil.copyfileobj(infile, out)

def concatenate_unix(outfile, pattern):
    cmd = ['cat'] + sorted(glob.glob(pattern))
    with open(outfile, 'wb') as out:
        subprocess.run(cmd, stdout=out)
                
def crop(file1, file2, args, processed_directory, out_nm):
    out_file1 = os.path.join(processed_directory,
                            f'{out_nm}_R1_cropped.fq.gz')
    out_file2 = os.path.join(processed_directory,
                            f'{out_nm}_R2_cropped.fq.gz')
    
    for f in [out_file1, out_file2]:
        if os.path.exists(f):
            os.remove(f)
    
    params = init_params(args)
    
    is_zipped = zipped(params['file1'])
    
    chunks1 = list(fastq_reader(params['file1'], params['chunk_size']))
    chunks2 = list(fastq_reader(params['file2'], params['chunk_size']))
    
    patterns = compile_patterns(params['ltr3'], params['linker3'], params['ltr5'], params['linker5'],
                                params['ltr3_error'], params['linker3_error'],
                                params['ltr5_error'], params['linker5_error'])
    
    Parallel(n_jobs = params['nthr'])(
        delayed(process_reads_parallel)(
            chunk1, chunk2, patterns, params, is_zipped,
            out_nm, processed_directory, chunk_num
            )
        for chunk_num, (chunk1, chunk2) in enumerate(zip(chunks1, chunks2))
    )

    inputs = [
        (out_file1, os.path.join(processed_directory, f"{out_nm}_*_R1_tmp.fq.gz")),
        (out_file2, os.path.join(processed_directory, f"{out_nm}_*_R2_tmp.fq.gz"))
        ]
    
    op_sys = platform.system()
    if op_sys != 'Darwin' and op_sys != 'Linux':
        Parallel(n_jobs=2)(
            delayed(concatenate_files)(outfile=outfile, pattern=pattern)
            for outfile, pattern in inputs
        )
    else:
        Parallel(n_jobs=2)(
            delayed(concatenate_unix)(outfile=outfile, pattern=pattern)
            for outfile, pattern in inputs
        )
        
    for pattern in [
        os.path.join(processed_directory, f"{out_nm}_*_R1_tmp.fq.gz"),
        os.path.join(processed_directory, f"{out_nm}_*_R2_tmp.fq.gz")
        ]:
        for tmp_file in glob.glob(pattern):
            os.remove(tmp_file)