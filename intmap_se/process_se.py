import math
import regex
import pysam
import sys
import numpy as np
from intmap_se.utils_se import *
from collections import defaultdict

def process_bam(out_bam):
    bamfile = pysam.AlignmentFile(out_bam, 'rb')
    reads_set = set()
    reads = []

    n_reads = 0
    for read in bamfile.fetch():
        if read.query_name not in reads_set:
            reads_set.add(read.query_name)
            reads.append(read)
            n_reads += 1

    print(f'Number of aligned reads: {n_reads}', flush = True)
    
    if n_reads > 0:
        bamfile.close()
        return reads
    else:
        print(f'\nNo reads aligned to the target genome. Exiting.', flush = True)
        bamfile.close()
        sys.exit()

def process_read(read, aln_mismatch_rate, aln_indel_rate, min_mapq, 
                 no_mm, min_qual, match_after):
    
    cigar = read.cigarstring
    cig_split = regex.findall('\\d+|\\D', cigar)
    if cig_split[1] != 'M':
        return None
    if 'S' in cig_split: # Should never happen with end-to-end
        return None
        
    md = read.get_tag('MD')
    read_seq = read.seq
    ns = read_seq.count('N')

    first_match = regex.match(r'^([ACGT]+)?(\d+)', md) or 'N' in read_seq[:match_after]
    if first_match:
        start_mm = len(first_match.group(1) or '')
        start_n = read_seq[:match_after].count('N')
        if (start_mm + start_n) > match_after:
            return None
    
    as_tag = read.get_tag('AS')
    if read.has_tag('XS'):
        xs_tag = read.get_tag('XS')
        if xs_tag > as_tag:
            return None
    else:
        xs_tag = None
    
    mapq = read.mapping_quality
    
    multimapping = False
    if (xs_tag is not None and as_tag == xs_tag and mapq <= 1):
        multimapping = True
    
    if no_mm and multimapping:
        return None
  
    if not multimapping and mapq < min_mapq:
        return None

    q_array = np.array(read.query_qualities)
    mean_qual = np.mean(10 ** (q_array / -10))
    if mean_qual > min_qual:
        return None

    aln_len = read.query_alignment_length
    n_mm = math.floor(aln_len * aln_mismatch_rate)
    n_indel = math.floor(aln_len * aln_indel_rate)
    edit_dist = read.get_tag('NM')
    
    if edit_dist <= (n_mm + n_indel):
        if read.has_tag('XM'):
            # N's are counted in the XM tag
            mismatch = read.get_tag('XM')
        else:
            # N's are not counted in the MD tag
            mismatch = len(regex.findall('[ATCG]', md))
            mismatch += ns
        
        if mismatch <= n_mm:
            indel = edit_dist - mismatch
            if indel <= n_indel:

                strand = '-' if read.is_reverse else '+'
                                
                ltr_umi, linker_umi = read.get_tag('RX').split('-')

                if strand == '-':
                    seq = revcomp(read.query_sequence)
                else:
                    seq = read.query_sequence

                matches = read.get_tag('OX').split('-')
                ltr_match = matches[0]
                linker_match = matches[1]
                    
                read_name = read.query_name
                        
                return{
                    'reference_name': read.reference_name,
                    'reference_start': read.reference_start,
                    'reference_end': read.reference_end,
                    'query_name': read_name,
                    'strand': strand,
                    'ltr_umi': ltr_umi,
                    'linker_umi': linker_umi,
                    'duplicate_count': 1,
                    'sequence': seq,
                    'multimapping': multimapping,
                    'map_quality': mapq,
                    'mean_quality': mean_qual,
                    'ltr_match': ltr_match,
                    'linker_match': linker_match
                }
            else:
                return None
        else:
            return None
    else:
        return None
      
def process_bam_lr(out_bam):
    bamfile = pysam.AlignmentFile(out_bam, 'rb')
    reads_dict = defaultdict(list)
    # reads_set = set()
    reads = []

    for read in bamfile.fetch():
        aln_score = read.get_tag('AS')
        reads_dict[read.query_name].append((aln_score, read))
        
    # n_reads = 0
    # for read in bamfile.fetch():
    #     if read.query_name not in reads_set:
    #         reads_set.add(read.query_name)
    #         reads.append(read)
    #         n_reads += 1

    n_reads = 0
    for _, alignments in reads_dict.items():
        if len(alignments) == 1:
            reads.append(alignments[0][1])
            n_reads += 1
        else:
            primary_alignment = next((aln for score, aln in alignments if (not aln.is_secondary or aln.get_tag('tp') == 'P')), None)
            primary_as_score = primary_alignment.get_tag('AS')
            # minimap2 will assign multimapping MAPQ scores to reads with very close alignment scores.
            # intmap considers these 'equal'
            equal_score_alignments = [aln for score, aln in alignments if abs(primary_as_score - score) < 25]
            multi_best = len(equal_score_alignments) > 1
            if multi_best:
                primary_alignment.set_tag('XS', primary_as_score)
            
            reads.append(primary_alignment)
            n_reads += 1
        
    print(f'Number of aligned reads: {n_reads}', flush=True)
        
    if n_reads > 0:
        bamfile.close()
        return reads
    else:
        print(f'\nNo reads aligned to the target genome. Exiting.', flush=True)
        bamfile.close()
        sys.exit()
        
def process_read_lr(read, aln_mismatch_rate, aln_indel_rate, min_mapq, min_frag_len,
                    max_frag_len, no_mm, min_qual, match_after):
    
    cigartuples = read.cigartuples
    
    n_clipped = 0
    if cigartuples[0][0] == 4:  # 4 corresponds to 'S' (soft-clipping)
        n_clipped = cigartuples[0][1]
    
    if cigartuples[0][0] != 0:  # 0 corresponds to 'M'
      return None
    
    n_clipped_end = 0
    if cigartuples[-1][0] == 4:
        n_clipped_end = cigartuples[-1][1]
        
    md = read.get_tag('MD')
    read_seq = read.seq
    ns = read_seq.count('N')

    md_prefix = ''
    for char in md:
        if char.isalpha():
            md_prefix += char
        elif char.isdigit():
            break
    first_match = md_prefix or 'N' in read_seq[:match_after] or n_clipped > 0
    if first_match:
        start_mm = len(md_prefix)
        start_n = read_seq[:match_after].count('N')
        total_start_mm = start_mm + n_clipped + start_n
        if total_start_mm > match_after:
            return None
    
    # multimapping = True if read.has_tag('MA') else False
    as_tag = read.get_tag('AS')
    if read.has_tag('XS'):
        xs_tag = read.get_tag('XS')
        if xs_tag > as_tag:
            return None
    else:
        xs_tag = None
    
    mapq = read.mapping_quality
    
    multimapping = False
    if (xs_tag is not None and as_tag == xs_tag and mapq <= 1):
        multimapping = True
    
    if no_mm and multimapping:
        return None

    if not multimapping and mapq < min_mapq:
        return None

    q_array = np.array(read.query_qualities)
    mean_qual = np.mean(10 ** (q_array / -10))
    if mean_qual > min_qual:
        return None
    
    tlen = len(read_seq)
    if abs(tlen) > max_frag_len or abs(tlen) < min_frag_len:
        return None

    aln_len = read.query_alignment_length
    # Modify aln_len to include soft-clipped bases
    # Important for treating soft-clipped bases as mismatches
    aln_len += n_clipped
    aln_len += n_clipped_end
    n_mm = math.floor(aln_len * aln_mismatch_rate)
    n_indel = math.floor(aln_len * aln_indel_rate)
    edit_dist = read.get_tag('NM')
    
    if edit_dist <= (n_mm + n_indel):
        if read.has_tag('XM'):
            mismatch = read.get_tag('XM')
            mismatch += n_clipped
            mismatch += n_clipped_end
        else:
            mismatch = sum(1 for char in md if char in 'ATCG')
            mismatch += ns
            mismatch += n_clipped
            mismatch += n_clipped_end
        
        if mismatch <= n_mm:
            indel = edit_dist - mismatch
            if indel <= n_indel:

                strand = '-' if read.is_reverse else '+'
                                
                ltr_umi, linker_umi = read.get_tag('RX').split('-')

                if strand == '-':
                    seq = revcomp(read.query_sequence)
                else:
                    seq = read.query_sequence

                matches = read.get_tag('OX').split('-')
                ltr_match = matches[0]
                linker_match = matches[1]
                    
                read_name = read.query_name
                        
                return{
                    'reference_name': read.reference_name,
                    'reference_start': read.reference_start - n_clipped,
                    'reference_end': read.reference_end + n_clipped_end,
                    'query_name': read_name,
                    'strand': strand,
                    'ltr_umi': ltr_umi,
                    'linker_umi': linker_umi,
                    'duplicate_count': 1,
                    'sequence': seq,
                    'multimapping': multimapping,
                    'map_quality': mapq,
                    'mean_quality': mean_qual,
                    'ltr_match': ltr_match,
                    'linker_match': linker_match
                }
            else:
                return None
        else:
            return None
    else:
        return None
