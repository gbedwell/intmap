import math
import regex
import pysam
from intmap import crop
import os
import subprocess

def align_random(bwa_path, sam_path, bwa_idx_dir, bwa_idx_name, 
                file1, file2, name, min_frag_len, max_frag_len, nthr):
    
    bwa_idx = f'{bwa_idx_dir}/{bwa_idx_name}'
    final_out = name + '_aln.bam'
    index_out = name + '_aln.bai'
    
    seed_len = 20
    
    bwa_cmd1 = (f'{bwa_path} mem -t {nthr} -k {seed_len} -C -L 1000,1000 {bwa_idx} {file1} {file2} | ')
    bwa_cmd2 = f'{sam_path} view -@ {nthr} -h -f 2 -F 2316 -b | '
    bwa_cmd3 = f'{sam_path} sort -@ {nthr} -o {final_out} -'
    
    align_cmd = bwa_cmd1 + bwa_cmd2 + bwa_cmd3
    
    idx_cmd = f'{sam_path} index -b -@ {nthr} {final_out} -o {index_out}'
    
    print('Performing random alignment...')
    subprocess.call(align_cmd, shell=True)
    
    print('\n')
    print('Indexing output bam file...')
    subprocess.call(idx_cmd, shell=True)

def process_random_bam(out_bam):
    bamfile = pysam.AlignmentFile(out_bam, 'rb')
    reads_dict = {}
    read_pairs = []

    n_reads = 0
    n_pairs = 0
    for read in bamfile.fetch():
        n_reads += 1
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
                        n_pairs += 1

    print(f'Number of aligned reads: {n_reads}')
    
    if n_reads > 0:
        print(f'Number of proper pairs: {n_pairs} ({(((n_pairs * 2) / n_reads) * 100):.2f}%)')
        bamfile.close()
        return read_pairs
    else:
        print(f'\nNo reads aligned to the target genome. Exiting.')
        bamfile.close()
        sys.exit()

def process_random_read(read, aln_mismatch_rate, aln_indel_rate, max_frag_len,
                        min_mapq, no_mm, min_qual, match_after):
    
    # Check CIGAR string starts with match/mismatch
    cigar = read.cigarstring
    cig_split = regex.findall('\\d+|\\D', cigar)
    if cig_split[1] != 'M':
        return None
    if 'S' in cig_split:
        return None
    
    # Get MD tag
    md = read.get_tag('MD')
    
    # Check alignment scores
    as_tag = read.get_tag('AS')
    if read.has_tag('XS'):
        xs_tag = read.get_tag('XS')
        if xs_tag > as_tag:
            return None
    else:
        xs_tag = None
    
    mapq = read.mapping_quality
    
    # Define multimapping status
    multimapping = False
    if (xs_tag is not None and abs(as_tag - xs_tag) <= 10 and mapq <= 0):
        multimapping = True
    
    if no_mm and multimapping:
        return None
        
    if not multimapping and mapq < min_mapq:
        return None
    
    # Check base qualities
    q_array = np.array(read.query_qualities)
    mean_qual = np.mean(10 ** (q_array / -10))
    if mean_qual > min_qual:
        return None
    
    # Check fragment length
    tlen = read.template_length
    if abs(tlen) > max_frag_len:
        return None
    
    # Check mismatch and indel rates
    aln_len = read.query_alignment_length
    n_mm = math.floor(aln_len * aln_mismatch_rate)
    n_indel = math.floor(aln_len * aln_indel_rate)
    edit_dist = read.get_tag('NM')
    
    if edit_dist <= (n_mm + n_indel):
        strand = '-' if tlen < 0 else '+'
        return {
            'reference_name': read.reference_name,
            'reference_start': read.reference_start,
            'reference_end': read.reference_end,
            'query_name': read.query_name,
            'strand': strand,
            'tlen': read.tlen,
            'multimapping': multimapping,
            'mean_quality': mean_qual
        }
    return None
    
def process_random_pair(read1, read2, aln_mismatch_rate, aln_indel_rate, max_frag_len, 
                        min_mapq, no_mm, min_qual, match_after):
    read1_info = process_random_read(read = read1, 
                                    aln_mismatch_rate = aln_mismatch_rate,
                                    aln_indel_rate = aln_indel_rate,
                                    max_frag_len = max_frag_len,
                                    min_mapq = min_mapq,
                                    no_mm = no_mm,
                                    min_qual = min_qual,
                                    match_after = match_after)
    
    read2_info = process_random_read(read = read2,
                                    aln_mismatch_rate = aln_mismatch_rate,
                                    aln_indel_rate = aln_indel_rate,
                                    max_frag_len = max_frag_len,
                                    min_mapq = min_mapq,
                                    no_mm = no_mm,
                                    min_qual = min_qual,
                                    match_after = match_after)
        
    if read1_info and read2_info:
        return read1_info, read2_info
    else:
        return None, None