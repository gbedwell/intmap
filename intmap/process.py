import math
import pysam
import sys
import numpy as np
from intmap.utils import *
from collections import defaultdict

def process_bam(out_bam, single_end, long_read, remove_chr):
    bamfile = pysam.AlignmentFile(out_bam, 'rb')
    reads = []
    n_reads = 0
    n_pairs = 0 if not single_end else None
    if remove_chr:
        rm_reads = set()
    
    if single_end and not long_read:
        reads_set = set()
    else:
        reads_dict = defaultdict(list)


    for read in bamfile.fetch():
        if remove_chr:
            if read.reference_name in remove_chr:
                if not single_end:
                    if read.query_name in reads_dict:
                        reads_dict.pop(read.query_name)
                    else:
                        rm_reads.add(read.query_name)
                continue
        if long_read:
            aln_score = read.get_tag('AS')
            reads_dict[read.query_name].append((aln_score, read))
        elif single_end and not long_read:
            if read.query_name not in reads_set:
                reads_set.add(read.query_name)
                reads.append(read)
                n_reads += 1
        else:
            if remove_chr:
                if read.query_name in rm_reads:
                    continue
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
                            reads.append((read1, read2))
                            n_pairs += 1

    if long_read:
        n_reads = len(reads_dict)
        for _, alignments in reads_dict.items():
            if len(alignments) == 1:
                reads.append(alignments[0][1])
            else:
                primary_alignment = next(
                    (aln for _, aln in alignments if (not aln.is_secondary or aln.get_tag('tp') == 'P')),
                    None
                )
                if primary_alignment:
                    primary_as_score = primary_alignment.get_tag('AS')
                    # minimap2 will assign multimapping MAPQ scores to reads with very close alignment scores.
                    # intmap considers these 'equal'
                    equal_score_alignments = [
                        aln for score, aln in alignments if abs(primary_as_score - score) < 25
                    ]
                    multi_best = len(equal_score_alignments) > 1
                    if multi_best:
                        primary_alignment.set_tag('XS', primary_as_score)
                    reads.append(primary_alignment)

    print(f'Number of aligned reads: {n_reads}', flush = True)
    
    if n_reads > 0:
        if not single_end and not long_read:
            print(f'Number of proper pairs: {n_pairs} ({(((n_pairs * 2) / n_reads) * 100):.2f}%)', flush = True)
        bamfile.close()
        return reads
    else:
        print(f'\nNo reads aligned to the target genome. Exiting.', flush = True)
        bamfile.close()
        sys.exit()

def process_read(read, aln_mismatch_rate, aln_indel_rate, max_frag_len,
                 min_frag_len, min_mapq, no_mm, min_qual, match_after,
                 single_end, long_read):
    cigartuples = read.cigartuples
    n_clipped = 0
    n_clipped_end = 0
    
    if cigartuples[0][0] == 4:  # 4 corresponds to 'S' (soft-clipping)
        n_clipped = cigartuples[0][1]
    elif cigartuples[0][0] != 0:  # 0 corresponds to 'M'
        return None

    if cigartuples[-1][0] == 4:
        n_clipped_end = cigartuples[-1][1]
        
    if not long_read and (n_clipped > 0 or n_clipped_end > 0):
        return None
        
    md = read.get_tag('MD')
    read_seq = read.seq
    ns = read_seq.count('N')

    if read.is_read1:
        md_prefix = ''
        for char in md:
            if char.isalpha():
                md_prefix += char
            elif char.isdigit() and char != '0':
                break
            elif char == '^':
                continue
        first_match = len(md_prefix) > 0 or 'N' in read_seq[:match_after] or n_clipped > 0
        if first_match:
            start_mm = len(md_prefix)
            start_n = read_seq[:match_after].count('N')
            total_start_mm = start_mm + n_clipped + start_n
            if total_start_mm > match_after:
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

    tlen = read.query_length if single_end else read.template_length
    if abs(tlen) > max_frag_len or abs(tlen) < min_frag_len:
        return None

    aln_len = read.query_alignment_length
    # Add soft-clipped bases to aln_len -- important for treating soft-clipped bases as mismatches.
    aln_len += n_clipped
    aln_len += n_clipped_end
    n_mm = math.floor(aln_len * aln_mismatch_rate)
    n_indel = math.floor(aln_len * aln_indel_rate)
    edit_dist = read.get_tag('NM')
    
    if edit_dist <= (n_mm + n_indel):
        if read.has_tag('XM'):
            # N's are counted in the XM tag
            mismatch = read.get_tag('XM')
            mismatch += n_clipped
            mismatch += n_clipped_end
        else:
            # N's are not counted in the MD tag
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
                    'reference_start': read.reference_start,
                    'reference_end': read.reference_end,
                    'query_name': read_name,
                    'strand': strand,
                    'tlen': read.tlen,
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

def process_read_pair(read1, read2, aln_mismatch_rate, aln_indel_rate, max_frag_len, 
                      min_frag_len, min_mapq, no_mm, min_qual, match_after, single_end = False,
                      long_read = False):
    read1_info = process_read(
      read = read1, aln_mismatch_rate = aln_mismatch_rate, aln_indel_rate = aln_indel_rate,
      min_frag_len = min_frag_len, max_frag_len = max_frag_len, min_mapq = min_mapq,
      no_mm = no_mm, min_qual = min_qual, match_after = match_after, single_end = single_end,
      long_read = long_read
    )
    
    read2_info = process_read(
      read = read2, aln_mismatch_rate = aln_mismatch_rate, aln_indel_rate = aln_indel_rate,
      min_frag_len = min_frag_len, max_frag_len = max_frag_len, min_mapq = min_mapq,
      no_mm = no_mm, min_qual = min_qual, match_after = match_after, single_end = single_end,
      long_read = long_read
    )
        
    if read1_info and read2_info:
        return read1_info, read2_info
    else:
        return None, None