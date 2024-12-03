#!/usr/local/bin/python

import os 
import sys
import getopt
import regex
import gzip
import io
import math
from collections import defaultdict
from joblib import Parallel, delayed
from itertools import tee, count
import shutil
import subprocess
import platform

def revcomp(seq):
    trans = str.maketrans("acgtumrwsykvhdbnACGTUMRWSYKVHDBN", "TGCAAKYWSRMBDHVNTGCAAKYWSRMBDHVN")
    seq_list = list(seq)
    seq_list.reverse()
    rev_seq = ''.join(seq_list)
    revcomp_seq = str.translate(rev_seq, trans)
    return revcomp_seq

def check_crop_input(ltr3, linker3, ltr1_primer, ltr5, linker5,
                    contamination, no_crop, ltr5_error_rate, 
                    linker5_error_rate):
    
    ltr1_check = ltr1_primer.replace(' ', '')
    ltr1_char = regex.sub('[ATGC]', '', ltr1_check.upper())
    if ltr1_char != '':
        raise ValueError('Declared round 1 LTR primer may only contain A,T,G, and C nucleotides.')

    if len(ltr1_check) < 10:
        raise ValueError('Declared round 1 LTR primer must be at least 10 nucleotides long.')
        
    ltr2_check = ltr5.replace(' ', '')
    ltr2_char = regex.sub('[ATGC]', '', ltr2_check.upper())
    if ltr2_char != '':
        raise ValueError('Declared round 2 LTR primer may only contain A,T,G, and C nucleotides.')

    if len(ltr2_check) < 8:
        raise ValueError('Declared round 2 LTR primer must be at least 8 nucleotides long.')
    
    if ltr5_error_rate > 0.3:
        raise ValueError('Round 2 LTR primer error rate cannot be > 0.3.')
        
    lp_check = linker5.replace(' ', '')
    lp_char = regex.sub('[ATGC]', '', lp_check.upper())
    if lp_char != '':
        raise ValueError('Declared linker primer may only contain A,T,G, and C nucleotides.')

    if len(lp_check) < 8:
        raise ValueError('Declared linker primer must be at least 8 nucleotides long.')
    
    if linker5_error_rate > 0.3:
        raise ValueError('Linker primer error rate cannot be > 0.3.')
            
    if contamination is not None:
        for cc in contamination:
            cont = cc.replace(' ', '')
            cont_char = regex.sub('[ATGCN]', '', cont.upper())
            if cont_char != '':
                raise ValueError('Declared contamination(s) may only contain A,T,G,C, and N nucleotides.')
            
            if len(cont) < 10:
                raise ValueError('Declared contamination(s) must be at least 10 nucleotides long.')
        
    linker_check = linker3.replace(' ', '')
    linker_char = regex.sub('[ATGC]', '', linker_check.upper())
    if linker_char != '':
        raise ValueError('The given linker sequence may only contain A, T, G, and C nucleotides.')
    
    if no_crop is False:
        if len(linker_check) < 5:
            raise ValueError('The given linker sequence must be >= 5 nucleotides long.')
        
    ltr_check = ltr3.replace(' ', '')
    test_ltr = regex.sub('[ATGC]', '', ltr_check.upper())
    if test_ltr != '':
        raise ValueError('The given LTR sequence may only contain A, T, G, and C nucleotides.')
    
    if no_crop is False:
        if len(ltr_check) < 5:
            raise ValueError('The given LTR sequence must be >= 5 nucleotides long.')

def zipped(file):
    if os.path.exists(file):
        with open(file, 'rb') as zip_test:
            return zip_test.read(2) == b'\x1f\x8b'
    else:
        return file.endswith('.gz')
    
def ensure_open(file):
    if isinstance(file, (str, bytes)):
        if zipped(file):
            return gzip.open(file, 'rt')
        else:
            return open(file, 'rt')
    elif isinstance(file, io.IOBase) and not file.closed:
        return file
    else:
        raise ValueError("Invalid file input or file is closed")

def concatenate_files(input_files, output_file):       
    buffer_size = 1024 * 1024 * 10
    with gzip.open(output_file, 'wb') as outfile:
        for file_name in input_files:
            with gzip.open(file_name, 'rb') as infile:
                while True:
                    buffer = infile.read(buffer_size)
                    if not buffer:
                        break
                    outfile.write(buffer)

def open_file(filename, mode):
    if zipped(filename):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
        
def concatenate_unix(input_files, output_file):
    cat_cmd = f"cat {' '.join(input_files)} > {output_file}"
    subprocess.call(cat_cmd, shell=True)

def chunk_fastq(file1_path, file2_path, chunk_size):
    with open_file(file1_path, 'rt') as file1, open_file(file2_path, 'rt') as file2:
        while True:
            chunk = []
            for _ in range(chunk_size):
                head1 = file1.readline().strip().split(' ')[0]
                seq1 = file1.readline().strip()
                opt1 = file1.readline().strip()
                qual1 = file1.readline().strip()

                head2 = file2.readline().strip().split(' ')[0]
                seq2 = file2.readline().strip()
                opt2 = file2.readline().strip()
                qual2 = file2.readline().strip()

                if not head1 or not head2:
                    break

                chunk.append((head1, seq1, opt1, qual1, head2, seq2, opt2, qual2))
            
            if not chunk:
                break
            
            yield chunk

def prepare_pattern(term, error_rate, error_type = "e"):
    term = term.replace(' ', '').upper()
    max_mismatches = math.floor(error_rate * len(term))
    return f'({term}){{{"s" if error_type == "s" else "e"}<={max_mismatches}}}'

def prepare_exact_pattern(term):
    term = term.replace(' ', '').upper()
    return rf'({term})'

def find_priority_matches(patterns, sequence):
    for pattern in patterns:
        match = pattern.search(sequence)
        if match:
            return list(pattern.finditer(sequence))[-1]
    return None

def crop_chunk(chunk, ltr_regex_list, linker_regex_list, ltr_rc_regex, linker_rc_regex,
                min_frag_len, no_crop, name, short_nm, contam_patterns, ltr_umi_len,
                ltr_umi_offset, linker_umi_len, linker_umi_offset,
                crop1_file, crop2_file):

    with gzip.open(crop1_file, 'wt') as crop1, gzip.open(crop2_file, 'wt') as crop2:
        for head1, seq1, opt1, qual1, head2, seq2, opt2, qual2 in chunk:
            if head1 != head2:
                raise Exception('ERROR: FASTQ read names do not match.')

            ltr_hit = find_priority_matches(ltr_regex_list, seq1)
            linker_hit = find_priority_matches(linker_regex_list, seq2)
            
            if not ltr_hit or not linker_hit:
                continue
            
            ltr_coords = ltr_hit.span()
            ltr_start, ltr_end = ltr_coords
            ltr_found = seq1[ltr_start:ltr_end]
            ltr_len = len(ltr_found)

            linker_coords = linker_hit.span()
            linker_start, linker_end = linker_coords
            linker_found = seq2[linker_start:linker_end]
            linker_len = len(linker_found)
            
            if not linker_hit:
                continue
            
            ltr_umi = seq1[(ltr_start - ltr_umi_offset - ltr_umi_len):(ltr_start - ltr_umi_offset)] if ltr_umi_len > 0 else 'N'
            linker_umi = seq2[(linker_start - linker_umi_offset - linker_umi_len):(linker_start - linker_umi_offset)] if linker_umi_len > 0 else 'N'
            
            linker_in_R1 = linker_rc_regex.search(seq1)
            if linker_in_R1:
                rlim_seq1 = linker_in_R1.span()[0]    
            else:
                rlim_seq1 = len(seq1)

            ltr_in_R2 = ltr_rc_regex.search(seq2)
            if ltr_in_R2:
                rlim_seq2 = ltr_in_R2.span()[0]
            else:
                rlim_seq2 = len(seq2)

            if no_crop:
                seq1_cropped = seq1[ltr_start:rlim_seq1]
                seq2_cropped = seq2[linker_start:rlim_seq2]
                contam_check = seq1_cropped[ltr_len:len(seq1_cropped)]
                qual1_cropped = qual1[ltr_start:rlim_seq1]
                qual2_cropped = qual2[linker_start:rlim_seq2]
            else:
                seq1_cropped = seq1[ltr_start + ltr_len:rlim_seq1]
                seq2_cropped = seq2[linker_start + linker_len:rlim_seq2]
                contam_check = seq1_cropped
                qual1_cropped = qual1[ltr_start + ltr_len:rlim_seq1]
                qual2_cropped = qual2[linker_start + linker_len:rlim_seq2]
                
            k = 0
            toss = False
            while k < len(contam_patterns) and toss != True:
                element = contam_patterns[k]
                if regex.match(element, contam_check):
                    toss = True 
                k += 1

            if ((len(seq1_cropped) < (min_frag_len + (ltr_len if no_crop is True else 0))) | 
                (len(seq2_cropped) < (min_frag_len + (linker_len if no_crop is True else 0))) |
                (toss == True)):
                continue

            # Write cropped reads to output
            crop1.write(f'{head1}\tCO:Z:1:{short_nm}\tRX:Z:{ltr_umi}-{linker_umi}\tOX:Z:{ltr_found}-{linker_found}\n')
            crop1.write(f'{seq1_cropped}\n{opt1}\n{qual1_cropped}\n')

            crop2.write(f'{head2}\tCO:Z:2:{short_nm}\tRX:Z:{ltr_umi}-{linker_umi}\tOX:Z:{ltr_found}-{linker_found}\n')
            crop2.write(f'{seq2_cropped}\n{opt2}\n{qual2_cropped}\n')

def crop(file1, file2, ltr3, virus, linker3, ltr1_primer, ltr5, linker5, 
        contamination, remove_internal_artifacts, U3, name, min_frag_len, no_crop, ltr3_error_rate, 
        linker3_error_rate, ltr5_error_rate, linker5_error_rate, ltr_umi_offset, ltr_umi_len, linker_umi_offset, 
        linker_umi_len, chunk_size, nthr, processed_dir):

    short_nm = name.rpartition('/')[-1] if '/' in name else name

    if remove_internal_artifacts:
        if virus is None:
            raise ValueError('Cannot define standard artifacts when the virus is not defined.')

        virus_artifacts = {
            'HIV1': ['GTCCCCCCTTTTCTT' if U3 else 'GTGGCGCCCGAA'],
            'HIV-1': ['GTCCCCCCTTTTCTT' if U3 else 'GTGGCGCCCGAA'],
            'MVV': ['GCTGGCGCCCAA'],
            'MLV': ['TTTGGGGGCTCG']
        }

        art = virus_artifacts.get(virus, None)
        if art is None:
            raise ValueError(
                "The stated virus is not a pre-defined virus type.\n"
                f"Pre-defined virus types are: {', '.join(virus_artifacts.keys())}"
            )
    else:
        art = []

    if contamination:
        art.extend(contamination)

    ltr3_pattern_exact = prepare_exact_pattern(ltr3)
    ltr5_pattern_exact = prepare_exact_pattern(ltr5)
    ltr_regex_exact = regex.compile(ltr5_pattern_exact + ltr3_pattern_exact)

    ltr3_pattern_sub = prepare_pattern(ltr3, ltr3_error_rate, error_type = "s")
    ltr5_pattern_sub = prepare_pattern(ltr5, ltr5_error_rate, error_type = "s")
    ltr_regex_sub = regex.compile(ltr5_pattern_sub + ltr3_pattern_sub)

    ltr3_pattern_indel = prepare_pattern(ltr3, ltr3_error_rate, error_type = "e")
    ltr5_pattern_indel = prepare_pattern(ltr5, ltr5_error_rate, error_type = "e")
    ltr_regex_indel = regex.compile(ltr5_pattern_indel + ltr3_pattern_indel)
    
    ltr_regex_list = [ltr_regex_exact, ltr_regex_sub, ltr_regex_indel]
    
    ltr_rc = revcomp(ltr5 + ltr3)[:11]
    ltr_rc_regex = regex.compile(f'({ltr_rc}){{e<=1}}')

    linker3_pattern_exact = prepare_exact_pattern(linker3)
    linker5_pattern_exact = prepare_exact_pattern(linker5)
    linker_regex_exact = regex.compile(linker5_pattern_exact + linker3_pattern_exact)

    linker3_pattern_sub = prepare_pattern(linker3, linker3_error_rate, error_type = "s")
    linker5_pattern_sub = prepare_pattern(linker5, linker5_error_rate, error_type = "s")
    linker_regex_sub = regex.compile(linker5_pattern_sub + linker3_pattern_sub)

    linker3_pattern_indel = prepare_pattern(linker3, linker3_error_rate, error_type = "e")
    linker5_pattern_indel = prepare_pattern(linker5, linker5_error_rate, error_type = "e")
    linker_regex_indel = regex.compile(linker5_pattern_indel + linker3_pattern_indel)
    
    linker_regex_list = [linker_regex_exact, linker_regex_sub, linker_regex_indel]

    linker_rc = revcomp(linker5 + linker3)[:11]
    linker_rc_regex = regex.compile(f'({linker_rc}){{e<=1}}')

    if art:
        contam_patterns = [regex.compile(f'({cc}){{e<={math.floor(len(cc) * 0.3)}}}') for cc in art]

    chunk_generator = chunk_fastq(file1, file2, chunk_size = chunk_size)
    chunk_gen, chunk_count = tee(chunk_generator)
    n_chunks = sum(1 for _ in chunk_count)
    
    tmp_dir = os.path.join(processed_dir, f'{name}_crop_tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        
    tmp_files_R1 = []
    tmp_files_R2 = []
    for i in range(n_chunks):
        tmp_files_R1.append(os.path.join(tmp_dir, f"{name}_tmp{str(i).zfill(3)}_R1_cropped.fq.gz"))
        tmp_files_R2.append(os.path.join(tmp_dir, f"{name}_tmp{str(i).zfill(3)}_R2_cropped.fq.gz"))

    Parallel(n_jobs = nthr)(
        delayed(crop_chunk)(
            chunk, ltr_regex_list, linker_regex_list, ltr_rc_regex, linker_rc_regex,
            min_frag_len, no_crop, name, short_nm, contam_patterns, ltr_umi_len,
            ltr_umi_offset, linker_umi_len, linker_umi_offset,
            tmp_files_R1[i], tmp_files_R2[i]
        )
        for i, chunk in enumerate(chunk_gen)
    )
    
    inputs_R1 = sorted(tmp_files_R1)
    inputs_R2 = sorted(tmp_files_R2)

    crop_out1 = os.path.join(processed_dir, f"{name}_R1_cropped.fq.gz")
    crop_out2 = os.path.join(processed_dir, f"{name}_R2_cropped.fq.gz")
    
    inputs = [inputs_R1, inputs_R2]
    outputs = [crop_out1, crop_out2]
    
    print('Concatenating cropped output...')
    
    op_sys = platform.system()
    if op_sys != 'Darwin' and op_sys != 'Linux':
        if nthr > 1:
            Parallel(n_jobs = 2)(
                delayed(concatenate_files)(input_files = inputs[i],
                                            output_file = outputs[i])
                for i in range(2))
        else:
            for i in range(2):
                concatenate_files(input_files = inputs[i],
                                output_file = outputs[i])
    else:
        if nthr > 1:
            Parallel(n_jobs = 2)(
                delayed(concatenate_unix)(input_files = inputs[i],
                                        output_file = outputs[i])
                for i in range(2))
        else:
            for i in range(2):
                concatenate_unix(input_files = inputs[i],
                                output_file = outputs[i])
        
    for f in inputs_R1:
        os.remove(f)
    for f in inputs_R2:
        os.remove(f)
        
    shutil.rmtree(tmp_dir)