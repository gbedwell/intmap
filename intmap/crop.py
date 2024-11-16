#!/usr/local/bin/python

import os 
import sys
import getopt
import regex
import gzip
import io
import math

def revcomp(seq):
    trans = str.maketrans("acgtumrwsykvhdbnACGTUMRWSYKVHDBN", "TGCAAKYWSRMBDHVNTGCAAKYWSRMBDHVN")
    seq_list = list(seq)
    seq_list.reverse()
    rev_seq = ''.join(seq_list)
    revcomp_seq = str.translate(rev_seq, trans)
    return revcomp_seq

def check_crop_input(ltr_term, linker_term, ltr1_primer, ltr2_primer, linker_primer,
                    contamination, no_crop, ltr2_primer_error_rate, 
                    linker_primer_error_rate):
    
    ltr1_check = ltr1_primer.replace(' ', '')
    ltr1_char = regex.sub('[ATGC]', '', ltr1_check.upper())
    if ltr1_char != '':
        raise ValueError('Declared round 1 LTR primer may only contain A,T,G, and C nucleotides.')

    if len(ltr1_check) < 10:
        raise ValueError('Declared round 1 LTR primer must be at least 10 nucleotides long.')
        
    ltr2_check = ltr2_primer.replace(' ', '')
    ltr2_char = regex.sub('[ATGC]', '', ltr2_check.upper())
    if ltr2_char != '':
        raise ValueError('Declared round 2 LTR primer may only contain A,T,G, and C nucleotides.')

    if len(ltr2_check) < 10:
        raise ValueError('Declared round 2 LTR primer must be at least 10 nucleotides long.')
    
    if ltr2_primer_error_rate > 0.3:
        raise ValueError('Round 2 LTR primer error rate cannot be > 0.3.')
        
    lp_check = linker_primer.replace(' ', '')
    lp_char = regex.sub('[ATGC]', '', lp_check.upper())
    if lp_char != '':
        raise ValueError('Declared linker primer may only contain A,T,G, and C nucleotides.')

    if len(lp_check) < 10:
        raise ValueError('Declared linker primer must be at least 10 nucleotides long.')
    
    if linker_primer_error_rate > 0.3:
        raise ValueError('Linker primer error rate cannot be > 0.3.')
            
    if contamination is not None:
        for cc in contamination:
            cont = cc.replace(' ', '')
            cont_char = regex.sub('[ATGCN]', '', cont.upper())
            if cont_char != '':
                raise ValueError('Declared contamination(s) may only contain A,T,G,C, and N nucleotides.')
            
            if len(cont) < 10:
                raise ValueError('Declared contamination(s) must be at least 10 nucleotides long.')
        
    linker_check = linker_term.replace(' ', '')
    linker_char = regex.sub('[ATGC]', '', linker_check.upper())
    if linker_char != '':
        raise ValueError('The given linker sequence may only contain A, T, G, and C nucleotides.')
    
    if no_crop is False:
        if len(linker_check) < 5:
            raise ValueError('The given linker sequence must be >= 5 nucleotides long.')
        
    ltr_check = ltr_term.replace(' ', '')
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

def crop(file1, file2, ltr_term, virus, linker_term, ltr1_primer, ltr2_primer, linker_primer, 
        contamination, remove_internal_artifacts, U3, name, min_frag_len, no_crop, ltr_error_rate, linker_error_rate, 
        ltr2_primer_error_rate, linker_primer_error_rate, ltr_umi_offset, ltr_umi_len, linker_umi_offset, 
        linker_umi_len, no_indels):
    
    file1 = ensure_open(file1)
    file2 = ensure_open(file2)
    
    if '/' in name:
        short_nm = name.rpartition('/')[-1]
    else:
        short_nm = name
    
    if remove_internal_artifacts:
        if virus is None:
            raise ValueError('Cannot define standard artifacts when the virus is not defined.')
        else:
            art = []
            if virus == 'HIV1' or virus == 'HIV-1':
                if U3 is False:
                    art.append('GTGGCGCCCGAA')
                else:
                    art.append('GTCCCCCCTTTTCTT')
            elif virus == 'MVV':
                art.append('GCTGGCGCCCAA')
            elif virus == 'MLV':
                art.append('TTTGGGGGCTCG')
            else:
                raise ValueError('The stated virus is not a pre-defined virus-type.',
                    '\n',
                    'Pre-defined virus types are: HIV1, HIV-1, MVV, and MLV'
                    '\n')
    
    try: 
        art
    except NameError:
        art = None
    
    if art is not None and contamination is not None:
        art.extend(contamination)
        
    linker_term = linker_term.replace(' ', '').upper()
    ltr_term = ltr_term.replace(' ', '').upper()
    ltr2_primer = ltr2_primer.replace(' ', '').upper()
    linker_primer = linker_primer.replace(' ', '').upper()
    
    ltr_term_mm = math.floor(ltr_error_rate * len(ltr_term))
    ltr_term_pattern = '({}){{{}<={}}}'.format(ltr_term, 's' if no_indels is True else 'e', ltr_term_mm)
    ltr2_primer_mm = math.floor(ltr2_primer_error_rate * len(ltr2_primer))
    ltr2_primer_pattern = '({}){{{}<={}}}'.format(ltr2_primer, 's' if no_indels is True else 'e', ltr2_primer_mm)
    ltr_pattern = ltr2_primer_pattern + ltr_term_pattern
    ltr_regex = regex.compile(ltr_pattern)
        
    linker_term_mm = math.floor(linker_error_rate * len(linker_term))
    linker_term_pattern = '({}){{{}<={}}}'.format(linker_term, 's' if no_indels is True else 'e', linker_term_mm)
    linker_primer_mm = math.floor(linker_primer_error_rate * len(linker_primer))
    linker_primer_pattern = '({}){{{}<={}}}'.format(linker_primer, 's' if no_indels is True else 'e', linker_primer_mm)
    linker_pattern = linker_primer_pattern + linker_term_pattern
    linker_regex = regex.compile(linker_pattern)
    
    if art is not None:
        contam = []
        contam_count = []
        for cc in art:
            contam_error = math.floor(len(cc) * 0.3)
            contam_pattern = '({}){{e<={}}}'.format(cc, contam_error)
            contam.append(regex.compile(contam_pattern))
            contam_count.append(0)
        
    crop1 = gzip.open(name + '_R1_cropped.fq.gz', 'wt')
    crop2 = gzip.open(name + '_R2_cropped.fq.gz', 'wt')
    
    nl = 0
    keep = False
    ltr_pos = -1
    ltr_rc_pos = -1
    linker_pos = -1
    linker_rc_pos = -1
    
    for line_r1, line_r2 in zip(file1, file2):
        lf1 = line_r1.strip()
        lf2 = line_r2.strip()
        if nl % 4 == 0:                
            if lf1.split(' ')[0] == lf2.split(' ')[0]:
                seq_name = lf1.split(' ')[0]
            else:
                raise Exception('ERROR: FASTQ read names do not match.')
                
        elif nl % 4 == 1:
            seq1 = lf1.strip('\n')
            seq2 = lf2.strip('\n')
            
            ltr_hit = ltr_regex.search(seq1)
            if ltr_hit:
                ltr_coords = list(ltr_regex.finditer(seq1))[-1].span()
                ltr_pos = ltr_coords[0]
                ltr_len = ltr_coords[1] - ltr_coords[0]
                ltr_found = seq1[ltr_coords[0]:ltr_coords[1]]
            else:
                continue

            linker_hit = linker_regex.search(seq2)
            if linker_hit:
                linker_coords = list(linker_regex.finditer(seq2))[-1].span()
                linker_pos = linker_coords[0]
                linker_len = linker_coords[1] - linker_coords[0]
                linker_found = seq2[linker_coords[0]:linker_coords[1]]
            else:
                continue
                        
            ltr_rc = revcomp(ltr_term)[0:10]
            ltr_found_rc = revcomp(ltr_found)[0:10]
            ltr_rc_pattern = '({}){{s<=1}}|({}){{s<=1}}'.format(ltr_rc, ltr_found_rc)
            ltr_rc_regex = regex.compile(ltr_rc_pattern)
            
            linker_rc = revcomp(linker_term)[0:10]
            linker_found_rc = revcomp(linker_found)[0:10]
            linker_rc_pattern = '({}){{s<=1}}|({}){{s<=1}}'.format(linker_rc, linker_found_rc)
            linker_rc_regex = regex.compile(linker_rc_pattern)
                
            if ltr_umi_len > 0:
                ltr_umi_end = ltr_pos - ltr_umi_offset
                ltr_umi = seq2[(ltr_umi_end - ltr_umi_len):ltr_umi_end]
            else:
                ltr_umi = 'N'
                
            if linker_umi_len > 0:
                linker_umi_end = linker_pos - linker_umi_offset
                linker_umi = seq2[(linker_umi_end - linker_umi_len):linker_umi_end]
            else:
                linker_umi = 'N'
                
            if ((ltr_pos > 0) & (linker_pos > 0)):
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
                    
                if no_crop is True:
                    seq1_cropped = seq1[ltr_pos:rlim_seq1]
                    seq2_cropped = seq2[linker_pos:rlim_seq2]
                    contam_check = seq1_cropped[ltr_len:len(seq1_cropped)]
                else:
                    seq1_cropped = seq1[ltr_pos + ltr_len:rlim_seq1]
                    seq2_cropped = seq2[linker_pos + linker_len:rlim_seq2]
                    contam_check = seq1_cropped
                try: 
                    contam
                except NameError:
                    contam = []
                    
                k = 0
                toss = False
                
                while k < len(contam) and toss != True:
                    element = contam[k]
                    if regex.match(element, contam_check):
                        toss = True 
                        contam_count[k] += 1 
                    k += 1

                if ((len(seq1_cropped) >= (min_frag_len + (ltr_len if no_crop is True else 0))) & 
                    (len(seq2_cropped) >= (min_frag_len + (linker_len if no_crop is True else 0))) & 
                    (toss == False)):
                    keep = True
                    crop1.write(seq_name + '\tCO:Z:1:' + short_nm + '\tRX:Z:' + 
                                ltr_umi + '-' + linker_umi + '\tOX:Z:' +
                                ltr_found + '-' + linker_found + '\n')
                    crop1.write(seq1_cropped + '\n')      
                    crop2.write(seq_name + '\tCO:Z:2:' + short_nm + '\tRX:Z:' + 
                                ltr_umi + '-' + linker_umi + '\tOX:Z:' +
                                ltr_found + '-' + linker_found + '\n')
                    crop2.write(seq2_cropped + '\n')
            
                else:
                    keep = False
                    ltr_pos = -1
                    linker_pos = -1
                    ltr_rc_pos = -1
                    linker_rc_pos = -1
            else:
                keep = False
                ltr_pos = -1
                linker_pos = -1
                ltr_rc_pos = -1
                linker_rc_pos = -1
            
        elif ((nl % 4 == 2) & (keep == True)):
            crop1.write('+' + '\n')
            crop2.write('+' + '\n')
            
        elif ((nl % 4 == 3) & (keep == True)):
            if no_crop is True:
                qual1 = lf1.strip('\n')[ltr_pos:rlim_seq1]
                qual2 = lf2.strip('\n')[linker_pos:rlim_seq2]
            else:
                qual1 = lf1.strip('\n')[ltr_pos + ltr_len:rlim_seq1]
                qual2 = lf2.strip('\n')[linker_pos + linker_len:rlim_seq2]
            crop1.write(qual1 + '\n')
            crop2.write(qual2 + '\n')
            
        nl += 1
            
    crop1.close()
    crop2.close()
    file1.close()
    file2.close()