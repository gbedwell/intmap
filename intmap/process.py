import math
import regex
import pysam
from intmap import crop

# Read in BAM file and check R1 and R2 for
# read name, location, and orientation.
def process_bam(out_bam):
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
                        
    print('Number of aligned reads:', n_reads)
    print('Number of proper pairs:', n_pairs)
    print(f'Proper pair percentage: {(((n_pairs * 2) / n_reads) * 100):.2f}%')

    bamfile.close()
    return read_pairs

# Perform QC on each read.
# Check mismatch rate, indel rate, fragment length, clip rate, MAPQ, and base quality.abs
# Return None if read does not pass QC.
# Also defines multimapping vs. unique reads, fragment UMIs, sequenced portion of fragments, etc.
def process_read(read, aln_mismatch_rate, aln_indel_rate, max_frag_len, 
                clip_rate, min_mapq, U3, no_mm, min_qual, match_after):
    
    # Ensure that each read starts on a match/mismatch (M)
    # (i.e., not an insertion, deletion, or is not soft-clipped)
    # Apply to both ends -- important for mapping both integration site and fragmentation breakpoint.
    cigar = read.cigarstring
    cig_split = regex.findall('\\d+|\\D', cigar)
    if read.is_forward:
        if cig_split[1] != 'M':
            return None
    elif read.is_reverse:
        if cig_split[-1] != 'M':
            return None
    
    # Check that the 5' end of read 1 matches the reference within
    # match_after aligned positions.
    if read.is_read1:
        md = read.get_tag('MD')
        md_split = regex.findall('[1-9]\\d*|[ACGT]', md)
        if len(md_split) > 1 and not md_split[0].isdigit():
            start_mm = 0
            for entry in md_split[0:]:
                if entry in 'ACGT':
                    start_mm += 1
                else:
                    break
            if start_mm > match_after:
                return None
    
    # Calculate the soft-clipped ratio
    cig_tups = read.cigartuples
    cig_counts = {i: 0 for i in range(10)}
    for identifier, count in cig_tups:
        cig_counts[identifier] += count
    count_m = cig_counts[0]
    count_i = cig_counts[1]
    count_d = cig_counts[2]
    count_s = cig_counts[4]
    
    clip_ratio = count_s / (count_m + count_i - count_d + count_s)
    if clip_ratio > clip_rate:
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
    
    # Quantify base quality scores
    # Work in probabilities, not Phred scores
    q_array = read.query_qualities
    qs = []
    for n in q_array:
        qs.append(10 ** (n / -10))
    # mean_qual = -10 * math.log10(sum(qs) / len(qs))
    mean_qual = sum(qs) / len(qs)
    if mean_qual > min_qual:
        return None
    
    # Filter by fragment length
    tlen = read.template_length
    if abs(tlen) > max_frag_len:
        return None
    
    # Filter by mismatch and indel rates
    # Quantified relative to alignment length
    aln_len = read.query_alignment_length
    n_mm = math.floor(aln_len * aln_mismatch_rate)
    n_indel = math.floor(aln_len * aln_indel_rate)
    edit_dist = read.get_tag('NM')
    
    if edit_dist <= (n_mm + n_indel):
        mismatch = read.get_tag('XM')
        
        if mismatch <= n_mm:
            indel = edit_dist - mismatch
            if indel <= n_indel:
                
                # All filtering is done at this point
                # Fill-in read information
                strand = '-' if tlen < 0 else '+'
                
                bcs = read.get_tag('RX').split('-')
                
                if U3:
                    # Define LTR UMI relative to provirus orientation
                    # and linker UMI relative to the top strand.
                    # Because cropping operates on the sequenced strand,
                    # UMIs from reverse strand must be corrected for orientation.
                    ltr_umi = crop.revcomp(bcs[0])
                    linker_umi = bcs[1]
                else:
                    ltr_umi = bcs[0]
                    linker_umi = crop.revcomp(bcs[1])

                # Define sequences relative to sequenced strand
                # Since BAM files report everything relative to the forward strand,
                # sequences must be converted back to proper strand.
                if strand == '-':
                    seq = crop.revcomp(read.query_sequence)
                else:
                    seq = read.query_sequence

                # Report LTR/linker matches as-sequenced.
                # This allows easy assessment of proper LTR/linker matches based on
                # the expected sequence.
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

# Perform QC on each read pair
def process_read_pair(read1, read2, aln_mismatch_rate, aln_indel_rate, max_frag_len, 
                        clip_rate, min_mapq, U3, no_mm, min_qual, match_after):
    read1_info = process_read(read = read1, 
                            aln_mismatch_rate = aln_mismatch_rate, 
                            aln_indel_rate = aln_indel_rate,
                            max_frag_len = max_frag_len,
                            clip_rate = clip_rate,
                            min_mapq = min_mapq,
                            U3 = U3,
                            no_mm = no_mm,
                            min_qual = min_qual,
                            match_after = match_after)
    
    read2_info = process_read(read=read2, 
                            aln_mismatch_rate = aln_mismatch_rate, 
                            aln_indel_rate = aln_indel_rate,
                            max_frag_len = max_frag_len,
                            clip_rate = clip_rate,
                            min_mapq = min_mapq,
                            U3 = U3,
                            no_mm = no_mm,
                            min_qual = min_qual,
                            match_after = match_after)
        
    if read1_info and read2_info:
        return read1_info, read2_info
    else:
        return None, None