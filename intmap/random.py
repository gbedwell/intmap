import math
import regex
import pysam
from intmap import crop
import os
import subprocess

def align_random(bt2_path, sam_path, bt2_idx_dir, bt2_idx_name, 
                file1, file2, name, min_frag_len, max_frag_len, nthr):
    
    bt2_idx = '{}/{}'.format(bt2_idx_dir, bt2_idx_name)
        
    final_out = name + '_aln.bam'
    index_out = name + '_aln.bai'

    bt2_cmd = ('{} --very-sensitive -x {} -1 {} -2 {} -p {} -I {} -X {} -f '
            '--sam-append-comment --no-mixed --no-discordant --no-unal | ').format(
            bt2_path, bt2_idx, file1, file2, nthr, min_frag_len, max_frag_len)
    
    bt2_to_bam = '{} view -@ {} -h -b | '.format(sam_path, nthr)
    sort_cmd = '{} sort -@ {} -o {} -'.format(sam_path, nthr, final_out)
    
    bt2_cmd = bt2_cmd + bt2_to_bam + sort_cmd
    
    idx_cmd = '{} index -b -@ {} -o {} {}'.format(sam_path, nthr, index_out, final_out)

    print('Aligning reads with bowtie2...')
    subprocess.call(bt2_cmd, shell=True)
    
    print('Indexing output bam file...')
    subprocess.call(idx_cmd, shell=True)

def process_random_bam(out_bam):
    bamfile = pysam.AlignmentFile(out_bam, 'rb')
    reads_dict = {}
    read_pairs = []

    for read in bamfile.fetch():
        if read.is_proper_pair:
            if read.query_name not in reads_dict:
                reads_dict[read.query_name] = read
            else:
                previous_read = reads_dict.pop(read.query_name)
                
                if read.is_read1:
                    read1, read2 = read, previous_read
                else:
                    read1, read2 = previous_read, read
                    
                if read1.reference_name == read2.reference_name:
                    if ((read1.is_forward and read2.is_reverse) or 
                        (read1.is_reverse and read2.is_forward)):
                        read_pairs.append((read1, read2))

    bamfile.close()
    return read_pairs

def process_random_read(read, max_frag_len, min_mapq, no_mm):
    
    cigar = read.cigarstring
    cig_split = regex.findall('\\d+|\\D', cigar)
    if cig_split[1] != 'M' and len(cig_split) != 2:
        return None
    
    # Make sure AS and XS tags make sense
    as_tag = read.get_tag('AS')
    if read.has_tag('XS'):
        xs_tag = read.get_tag('XS')
        if xs_tag > as_tag:
            return None
    else:
        xs_tag = None
    
    mapq = read.mapping_quality
    
    multimapping = False
    if (as_tag == xs_tag and (mapq == 1 or mapq == 0)):
        multimapping = True
    
    # Option to remove multimapping reads
    if no_mm and multimapping:
        return None
        
    # Filter by MAPQ
    # Only applied to uniquely mapping reads
    if not multimapping and mapq < min_mapq:
        return None
    
    tlen = read.template_length
    if abs(tlen) > max_frag_len:
        return None

    edit_dist = read.get_tag('NM')
    if edit_dist > 0:
        return None
    else:
        strand = '-' if tlen < 0 else '+'
        return {
            'reference_name': read.reference_name,
            'reference_start': read.reference_start,
            'reference_end': read.reference_end,
            'query_name': read.query_name,
            'strand': strand,
            'tlen': read.tlen,
            'multimapping': multimapping
        }
    
def process_random_pair(read1, read2, max_frag_len, min_mapq, no_mm):
    read1_info = process_random_read(read = read1, 
                                    max_frag_len = max_frag_len,
                                    min_mapq = min_mapq,
                                    no_mm = no_mm)
    
    read2_info = process_random_read(read = read2, 
                                    max_frag_len = max_frag_len,
                                    min_mapq = min_mapq,
                                    no_mm = no_mm)
        
    if read1_info and read2_info:
        return read1_info, read2_info
    else:
        return None, None