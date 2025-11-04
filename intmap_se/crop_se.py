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

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return ValueError('Sequences must be of equal length.')
    return sum(x1 != x2 for x1, x2 in zip(seq1, seq2))

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
            'perfect': regex.compile(revcomp(linker3) + revcomp(linker5)),
            'mismatch': regex.compile(f'({revcomp(linker3)}){{s<={linker3_errors}}}({revcomp(linker5)}){{s<={linker5_errors}}}'),
            'indel': regex.compile(f'({revcomp(linker3)}){{e<={linker3_errors}}}({revcomp(linker5)}){{e<={linker5_errors}}}')
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
    
def find_diagnostic_regions(ltr3, ltr3_alt, window, min_diff, search_limit):
    found_region = False
    for i in range(math.ceil(len(ltr3) * search_limit) - window + 1):
        ltr3_sub = ltr3[i:i + window]
        ltr3_alt_sub = ltr3_alt[i:i + window]
        distance = hamming_distance(ltr3_sub, ltr3_alt_sub)
        if distance >= min_diff:
            flag_info = (ltr3_sub, ltr3_alt_sub)
            found_region = True
            break
    if found_region:
        return flag_info
    else:
        return (ltr3[:window], ltr3_alt[:window])

def compile_patterns_ttr(ltr3, ltr3_alt, ltr5, min_ttr_len, ltr3_error_rate, 
                            ltr5_error_rate, linker3, linker5, linker3_error_rate,
                            linker5_error_rate, flag_window_size,
                            flag_min_diff, flag_search_limit):
    linker3_errors = math.floor(len(linker3) * linker3_error_rate)
    linker5_errors = math.floor(len(linker5) * linker5_error_rate)
    ltr5_errors = math.floor(len(ltr5) * ltr5_error_rate)
    
    linker_rc = revcomp(linker5 + linker3)[:15]
    if len(linker_rc) > 10:
        linker_rc_errors = 2
    else:
        linker_rc_errors = 1
    
    patterns = {
        'ltr3': [],
        'ltr3_alt': [],
        'ltr5': {
            'perfect': regex.compile(ltr5),
            'mismatch': regex.compile(f'({ltr5}){{s<={ltr5_errors}}}'),
            'indel': regex.compile(f'({ltr5}){{e<={ltr5_errors}}}')
            },
        'linker': {
            'perfect': regex.compile(revcomp(linker3) + revcomp(linker5)),
            'mismatch': regex.compile(f'({revcomp(linker3)}){{s<={linker3_errors}}}({revcomp(linker5)}){{s<={linker5_errors}}}'),
            'indel': regex.compile(f'({revcomp(linker3)}){{e<={linker3_errors}}}({revcomp(linker5)}){{e<={linker5_errors}}}')
            },
        'linker_rc': {
            'perfect': regex.compile(linker_rc),
            'mismatch': regex.compile(f'({linker_rc}){{s<={linker_rc_errors}}}'),
            'indel': regex.compile(f'({linker_rc}){{e<={linker_rc_errors}}}')
            }
        }
    
    if ltr3_alt is not None:
        flags = find_diagnostic_regions(
            ltr3 = ltr3, 
            ltr3_alt = ltr3_alt, 
            window = flag_window_size, 
            min_diff = flag_min_diff, 
            search_limit = flag_search_limit
            )
        
        if len(flags[0]) > 10:
            flag_errors = 2
        else:
            flag_errors = 1
            
        patterns.update(
            {
                'ltr3_flag': {
                    'perfect': regex.compile(flags[0]),
                    'mismatch': regex.compile(f'({flags[0]}){{s<={flag_errors}}}'),
                    'indel': regex.compile(f'({flags[0]}){{e<={flag_errors}}}')
                },
                'ltr3_flag_start': {
                    'perfect': regex.compile(ltr3[:flag_window_size]),
                    'mismatch': regex.compile(f'({ltr3[:flag_window_size]}){{s<={flag_errors}}}'),
                    'indel': regex.compile(f'({ltr3[:flag_window_size]}){{e<={flag_errors}}}')
                },
                'ltr3_alt_flag': {
                    'perfect': regex.compile(flags[1]),
                    'mismatch': regex.compile(f'({flags[1]}){{s<={flag_errors}}}'),
                    'indel': regex.compile(f'({flags[1]}){{e<={flag_errors}}}')
                },
                'ltr3_alt_flag_start': {
                    'perfect': regex.compile(ltr3_alt[:flag_window_size]),
                    'mismatch': regex.compile(f'({ltr3_alt[:flag_window_size]}){{s<={flag_errors}}}'),
                    'indel': regex.compile(f'({ltr3_alt[:flag_window_size]}){{e<={flag_errors}}}')
                }
            }
        )
    
    for length in range(len(ltr3), (min_ttr_len - 1), -1):
        ltr3_segment = ltr3[:length]
        ltr3_errors = math.floor(len(ltr3_segment) * ltr3_error_rate)
        full_ltr = ltr5 + ltr3_segment
        
        pattern_set = {
            'perfect': regex.compile(full_ltr),
            'mismatch': regex.compile(f'({ltr5}){{s<={ltr5_errors}}}({ltr3_segment}){{s<={ltr3_errors}}}'),
            'indel': regex.compile(f'({ltr5}){{e<={ltr5_errors}}}({ltr3_segment}){{e<={ltr3_errors}}}')
        }
        
        patterns['ltr3'].append(pattern_set)
        
    if ltr3_alt is not None:
        for length in range(len(ltr3_alt), (min_ttr_len - 1), -1):
            ltr3_alt_segment = ltr3_alt[:length]
            ltr3_alt_errors = math.floor(len(ltr3_alt_segment) * ltr3_error_rate)
            full_ltr_alt = ltr5 + ltr3_alt_segment
            
            alt_pattern_set = {
                'perfect': regex.compile(full_ltr_alt),
                'mismatch': regex.compile(f'({ltr5}){{s<={ltr5_errors}}}({ltr3_alt_segment}){{s<={ltr3_alt_errors}}}'),
                'indel': regex.compile(f'({ltr5}){{e<={ltr5_errors}}}({ltr3_alt_segment}){{e<={ltr3_alt_errors}}}')
            }
            
            patterns['ltr3_alt'].append(alt_pattern_set)
    
    return patterns
    
def compile_rc_ttr(ltr_match):
    ltr_rc = revcomp(ltr_match)[:15]
    
    if len(ltr_rc) > 10:
        ltr_rc_errors = 2
    else:
        ltr_rc_errors = 1
    
    return {
        'ltr_rc': {
            'perfect': regex.compile(ltr_rc),
            'mismatch': regex.compile(f'({ltr_rc}){{s<={ltr_rc_errors}}}'),
            'indel': regex.compile(f'({ltr_rc}){{e<={ltr_rc_errors}}}')
        }
    }

def find_pattern_match(seq, patterns, pattern_type, no_error):
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

def find_flag_match(patterns, pattern_key, start_pos, query_seq, no_error):
    flag_found = False
    
    flag_match = patterns[pattern_key]['perfect'].search(query_seq, pos = start_pos)
    
    if not flag_match:
        if not no_error:
            flag_match = patterns[pattern_key]['mismatch'].search(query_seq, pos = start_pos)
            if not flag_match:
                flag_match = patterns[pattern_key]['indel'].search(query_seq, pos = start_pos)
            if flag_match:
                flag_found = True
                
        if not flag_found:
            start_key = f"{pattern_key}_start"
            flag_match = patterns[start_key]['perfect'].search(query_seq, pos = start_pos)
            if not flag_match:
                if not no_error:
                    flag_match = patterns[start_key]['mismatch'].search(query_seq, pos = start_pos)
                    if not flag_match:
                        flag_match = patterns[start_key]['indel'].search(query_seq, pos = start_pos)
                    if flag_match:
                        flag_found = True
            else:
                flag_found = True
    else:
        flag_found = True
    
    return flag_found

def sum_term(sub_term, match_term):
    val1 = 1 if sub_term[0] == match_term[0] else 0
    val2 = 1 if sub_term[1] == match_term[1] else 0
    val3 = 2 if sub_term[2] == match_term[2] else 0
    
    return val1 + val2 + val3
    
def find_ttr(patterns, pattern_type, seq, no_error):
    match_found = False
    for subpat in patterns[pattern_type]:
        subseq = subpat['perfect'].pattern
        subseq_term = subseq[-3:]
        segment_match = subpat['perfect'].search(seq)
        if segment_match:
            match_found = True
            break
        elif not no_error:
            segment_match = subpat['mismatch'].search(seq)
            if segment_match and sum_term(subseq_term, segment_match.group()[-3:]) > 2:
                match_found = True
                break
            else:
                segment_match = None
            if not segment_match:
                segment_match = subpat['indel'].search(seq)
                if segment_match and sum_term(subseq_term, segment_match.group()[-3:]) > 2:
                    match_found = True
                    break
                else:
                    segment_match = None
        else:
            segment_match = None
            
    return segment_match

def find_pattern_match_ttr(seq, patterns, no_error, search_limit):
    ltr5_match = patterns['ltr5']['perfect'].search(seq)
    if not ltr5_match:
        if not no_error:
            ltr5_match = patterns['ltr5']['mismatch'].search(seq)
            if not ltr5_match:
                ltr5_match = patterns['ltr5']['indel'].search(seq)
    
    if not ltr5_match:
        return None
    
    ltr3_flag_found = False
    ltr3_alt_flag_found = False
    flag_query_start = ltr5_match.end()
    flag_query_end = ltr5_match.end() + math.ceil(len(patterns['ltr3'][0]['perfect'].pattern) * search_limit)
    flag_query_seq = seq[flag_query_start:flag_query_end]
    
    if 'ltr3_flag' in patterns:
        ltr3_flag_found = find_flag_match(
            patterns = patterns, 
            pattern_key = 'ltr3_flag', 
            start_pos = 0,
            query_seq = flag_query_seq, 
            no_error = no_error
            )
    else:
        ltr3_flag_found = True

    if 'ltr3_alt_flag' in patterns:
        ltr3_alt_flag_found = find_flag_match(
            patterns = patterns, 
            pattern_key = 'ltr3_alt_flag', 
            start_pos = 0,
            query_seq = flag_query_seq, 
            no_error = no_error
            )
            
    if not ltr3_flag_found and not ltr3_alt_flag_found:
        return None
    
    if ltr3_flag_found and not ltr3_alt_flag_found:
        segment_match = find_ttr(
            patterns = patterns,
            pattern_type = 'ltr3',
            seq = seq,
            no_error = no_error
            )
    elif not ltr3_flag_found and ltr3_alt_flag_found:
        segment_match = find_ttr(
                patterns = patterns,
                pattern_type = 'ltr3_alt',
                seq = seq,
                no_error = no_error
                )
    elif ltr3_flag_found and ltr3_alt_flag_found:     
        ltr3_match = find_ttr(
            patterns = patterns,
            pattern_type = 'ltr3',
            seq = seq,
            no_error = no_error
            )

        ltr3_alt_match = find_ttr(
            patterns = patterns,
            pattern_type = 'ltr3_alt',
            seq = seq,
            no_error = no_error
            )
        
        if ltr3_match and ltr3_alt_match:
            if len(ltr3_match.group()) >= len(ltr3_alt_match.group()):
                segment_match = ltr3_match
            else:
                segment_match = ltr3_alt_match
        elif ltr3_match:
            segment_match = ltr3_match
        elif ltr3_alt_match:
            segment_match = ltr3_alt_match
        else:
            segment_match = None
    
    return segment_match

# def detect_quality_offset(fastq_file, sample_size = 10000):
#     min_qual = float('inf')

#     with open_file(fastq_file, 'rt') as f:
#         for i, line in enumerate(f):
#             if i % 4 == 3:
#                 min_qual = min(min_qual, min(ord(c) for c in line.strip()))
#             if i >= sample_size * 4:
#                 break

#     return 64 if min_qual >= 64 else 33

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
        'linker_umi_len': args.linker_umi_len,
        'linker_umi_offset': args.linker_umi_offset,
        'linker_umi_pattern': args.linker_umi_pattern,
        'contam_patterns': contam_patterns,
        'ltr5': args.ltr5,
        'ltr3': args.ltr3,
        'linker5': args.linker5,        
        'linker3': args.linker3,
        'nthr': args.nthr,
        'chunk_size': args.crop_chunk_size,
        'ltr3_alt': args.ltr3_alt,
        'ttr': args.ttr,
        'min_ttr_len': args.min_ttr_len,
        'flag_window_size': args.flag_window_size,
        'flag_min_diff': args.flag_min_diff,
        'flag_search_limit': args.flag_search_limit
    }
    
def extract_umis(seq, pattern_match, umi_len, umi_offset,
                 umi_pattern = None, linker_end = False):

    if umi_len > 0:
        if not linker_end:
            start = pattern_match.start() - umi_offset - umi_len
            end = pattern_match.start() - umi_offset
        else:
            start = pattern_match.end() + umi_offset
            end = pattern_match.end() + umi_offset + umi_len

        umi = seq[start:end] if not linker_end else revcomp(seq[start:end])
        
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
                           out_nm, processed_directory, chunk_num, ttr, long_read):
    sequences1, qualities1, headers1 = process_fastq_chunk(chunk1, params, is_zipped)
    
    cropped_reads1 = []
    
    no_error = False
    if (params['ltr5_error'] == 0 and params['linker5_error'] == 0 and
        params['ltr3_error'] == 0 and params['linker3_error'] == 0):
        no_error = True

    for i, (seq1, qual1, head1) in enumerate(zip(sequences1, qualities1, headers1)):
        seq1_str = seq1
        
        if ttr:
            ltr_pattern_match = find_pattern_match_ttr(
                seq = seq1_str, 
                patterns = patterns, 
                no_error = no_error, 
                search_limit = params['flag_search_limit']
                )
            # if ltr_pattern_match:
            #     ttr_ltr_rc_pattern = compile_rc_ttr(ltr_pattern_match.group())
        else:
            ltr_pattern_match = find_pattern_match(
                seq = seq1_str, 
                patterns = patterns, 
                pattern_type = 'ltr', 
                no_error = no_error
                )
            
        if long_read:
            linker_pattern_match = find_pattern_match(
                    seq = seq1_str, 
                    patterns = patterns, 
                    pattern_type = 'linker', 
                    no_error = no_error
                    )
        
        if ltr_pattern_match and ((long_read and linker_pattern_match) or not long_read):
            start1 = ltr_pattern_match.end()
            
            if not long_read:
                linker_in_ltr = find_pattern_match(
                    seq = seq1_str, 
                    patterns = patterns, 
                    pattern_type = 'linker_rc', 
                    no_error = False
                    )

                if linker_in_ltr:
                    end1 = linker_in_ltr.start()
                    cropped_seq1 = seq1_str[start1:end1]
                    cropped_qual1 = qual1[start1:end1]
                else:
                    cropped_seq1 = seq1_str[start1:]
                    cropped_qual1 = qual1[start1:]
            else:
                end1 = linker_pattern_match.start()
                cropped_seq1 = seq1_str[start1:end1]
                cropped_qual1 = qual1[start1:end1]
                    
            ltr_umi = extract_umis(
                seq = seq1_str, 
                pattern_match = ltr_pattern_match, 
                umi_len = params['ltr_umi_len'], 
                umi_offset = params['ltr_umi_offset'],
                umi_pattern = params['ltr_umi_pattern'],
                linker_end = False
                )
            
            if long_read:
                linker_umi = extract_umis(
                    seq = seq1_str, 
                    pattern_match = linker_pattern_match, 
                    umi_len = params['linker_umi_len'], 
                    umi_offset = params['linker_umi_offset'],
                    umi_pattern = params['linker_umi_pattern'],
                    linker_end = True
                    )
            else:
                linker_umi = 'N'
            
            if not ltr_umi or (long_read and not linker_umi):
                continue
            
            toss = False
            if params['contam_patterns']:
                for pattern in params['contam_patterns']:
                    if regex.match(pattern, cropped_seq1):
                        toss = True
                        break

            if not toss and len(cropped_seq1) >= params['min_len']:
                ltr_found = ltr_pattern_match.group()
                linker_found = revcomp(linker_pattern_match.group()) if long_read else (revcomp(linker_in_ltr.group()) if linker_in_ltr else 'N')
                    
                new_header1 = (f'{head1}\tCO:Z:1:{out_nm}\t'
                              f'RX:Z:{ltr_umi}-{linker_umi}\t'
                              f'OX:Z:{ltr_found}-{linker_found}')
                
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
                
def crop(args, processed_directory, out_nm, long_read):
    out_file1 = os.path.join(processed_directory,
                            f'{out_nm}_R1_cropped.fq.gz')
    
    if os.path.exists(out_file1):
        os.remove(out_file1)
    
    params = init_params(args)
    
    is_zipped = zipped(params['file1'])
    
    chunks1 = list(fastq_reader(params['file1'], params['chunk_size']))
    
    if params['ttr']:
        patterns = compile_patterns_ttr(
            ltr3 = params['ltr3'],
            ltr3_alt = params['ltr3_alt'],
            ltr5 = params['ltr5'],
            min_ttr_len = params['min_ttr_len'],
            ltr3_error_rate = params['ltr3_error'],
            ltr5_error_rate =  params['ltr5_error'],
            linker3 = params['linker3'],
            linker5 = params['linker5'],
            linker3_error_rate = params['linker3_error'],
            linker5_error_rate = params['linker5_error'],
            flag_window_size = params['flag_window_size'], 
            flag_min_diff = params['flag_min_diff'], 
            flag_search_limit = params['flag_search_limit']
        )
    else:
        patterns = compile_patterns(
            ltr3 = params['ltr3'], 
            linker3 = params['linker3'], 
            ltr5 = params['ltr5'], 
            linker5 = params['linker5'],
            ltr3_error_rate = params['ltr3_error'], 
            linker3_error_rate = params['linker3_error'],
            ltr5_error_rate = params['ltr5_error'], 
            linker5_error_rate = params['linker5_error']
            )
    
    Parallel(n_jobs=params['nthr'])(
        delayed(process_reads_parallel)(
            chunk1 = chunk1, 
            patterns = patterns, 
            params = params, 
            is_zipped = is_zipped, 
            out_nm = out_nm,
            processed_directory = processed_directory, 
            chunk_num = chunk_num, 
            ttr = params['ttr'],
            long_read = long_read
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