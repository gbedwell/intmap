import math
import regex
import pysam
import sys
import numpy as np
from intmap_se.utils_se import *

def process_bam(out_bam):
    bamfile = pysam.AlignmentFile(out_bam, 'rb')
    reads_dict = {}
    reads = []

    n_reads = 0
    for read in bamfile.fetch():
        if read.query_name not in reads_dict:
            reads_dict[read.query_name] = read
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

# Perform QC on each read.
# Check mismatch rate, indel rate, fragment length, clip rate, MAPQ, and base quality.abs
# Return None if read does not pass QC.
# Also defines multimapping vs. unique reads, fragment UMIs, sequenced portion of fragments, etc.
def process_read(read, aln_mismatch_rate, aln_indel_rate, max_frag_len,
                min_frag_len, min_mapq, U3, no_mm, min_qual, match_after):
    
    # Ensure that each read starts on a match/mismatch (M)
    cigar = read.cigarstring
    cig_split = regex.findall('\\d+|\\D', cigar)
    if cig_split[1] != 'M':
        return None
    if 'S' in cig_split: # Should never happen with end-to-end
        return None
        
    # Get MD tag. Used in a couple of places downstream.
    md = read.get_tag('MD')
    
    # Check that the 5' end of read 1 matches the reference after
    # match_after aligned positions.
    first_match = regex.match(r'^([ACGT]+)?(\d+)', md)
    if first_match:
        start_mm = len(first_match.group(1) or '')
        if start_mm > match_after:
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
    if (xs_tag is not None and as_tag == xs_tag and mapq <= 1):
        multimapping = True
    
    # Option to remove multimapping reads
    if no_mm and multimapping:
        return None
        
    # Filter by MAPQ
    # Only applied to uniquely mapping reads
    if not multimapping and mapq < min_mapq:
        return None
    
    # Quantify base quality scores
    # Work in probabilities, not Phred scores
    q_array = np.array(read.query_qualities)
    mean_qual = np.mean(10 ** (q_array / -10))
    if mean_qual > min_qual:
        return None
    
    # Filter by mismatch and indel rates
    # Quantified relative to alignment length
    aln_len = read.query_alignment_length
    n_mm = math.floor(aln_len * aln_mismatch_rate)
    n_indel = math.floor(aln_len * aln_indel_rate)
    edit_dist = read.get_tag('NM')
    
    if edit_dist <= (n_mm + n_indel):
        if read.has_tag('XM'):
            mismatch = read.get_tag('XM')
        else:
            mismatch = len(regex.findall('[ATCG]', md))
        
        if mismatch <= n_mm:
            indel = edit_dist - mismatch
            if indel <= n_indel:
                
                # All filtering is done at this point
                # Fill-in read information
                strand = '-' if read.is_reverse else '+'
                                
                ltr_umi = read.get_tag('RX')

                # Define sequences relative to sequenced strand
                # Since BAM files report everything relative to the forward strand,
                # sequences must be converted back to proper strand.
                if strand == '-':
                    seq = revcomp(read.query_sequence)
                else:
                    seq = read.query_sequence

                # Report LTR/linker matches as-sequenced.
                # This allows easy assessment of proper LTR/linker matches based on
                # the expected sequence.
                matches = read.get_tag('OX')
                ltr_match = matches
                    
                read_name = read.query_name
                        
                return{
                    'reference_name': read.reference_name,
                    'reference_start': read.reference_start,
                    'reference_end': read.reference_end,
                    'query_name': read_name,
                    'strand': strand,
                    'ltr_umi': ltr_umi,
                    'duplicate_count': 1,
                    'sequence': seq,
                    'multimapping': multimapping,
                    'mean_quality': mean_qual,
                    'ltr_match': ltr_match
                }
            else:
                return None
        else:
            return None
    else:
        return None