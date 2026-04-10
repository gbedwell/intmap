# intmap
intmap is a python package for mapping the locations of genomic integration from second- and third-generation sequencing data.

## Current version
The latest version of intmap is v2.0.0. Versions of intmap <= v1.0.1 are considered obsolete and should not be used. It is recommended that users of older versions upgrade to intmap >= v2.0.0. See the CHANGELOG for more information.

## Installation
The easiest way to install intmap is by cloning this repository, navigating to the generated 'intmap' directory, and running `bash install.sh`. This bash script will:

1. Create an intmap conda environment,
2. install all (version-controlled) intmap dependencies, and
3. install the intmap package within the intmap conda environment.

`install.sh` requires an initialized conda instance. If conda is not available on your computer, you can easily install conda by following the recommendations [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). `install.sh` can additionally install new versions of intmap seamlessly from either the base environment or an activated intmap environment.

## Package Structure
There are 3 distinct intmap modules: intmap_demux, intmap, and intmap_multi. I will explain these in more detail below.

### intmap_demux
intmap_demux is written to demultiplex multiplexed FASTQ files. The software can demultiplex reads based on in-line indexes (barcodes) on either or both read ends, indexes stored in FASTQ headers, or indexes stored in separate FASTQ files. For in-line indexes, intmap_demux uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/#) to look for target sequences and allocate them to the appropriate file. For other index formats, intmap uses fuzzy regex matching for sample assignment. Sample-index pairs are defined in the input FASTA files. FASTA headers should denote each sample-of-interest and must be unique. The index sequence paired with each sample should be given as the sequence of each FASTA entry. All index sequences must also be unique.

NOTE: integrant-end sequences are universally referred to as 'ltr' in intmap syntax. However, these sequences do not have to be retrovirus/retrotransposon-specific. Any integrant-specific sequence is valid.

#### Examples
Given multiplexed paired-end FASTQ files file_R1.fq.gz and file_R2.fq.gz with in-line barcodes on both reads, intmap_demux can be run like so, assuming R1 contains integrant-end reads and R2 linker-end reads:

```
intmap_demux \
    -ltr_reads file_R1.fq.gz \
    -linker_reads file_R2.fq.gz \
    -bc_ltr ltr_bc.fa \
    -bc_linker linker_bc.fa
```

If indexes are stored in a separate file, say, file_I1.fq.gz:

```
intmap_demux \
    -ltr_reads file_R1.fq.gz \
    -linker_reads file_R2.fq.gz \
    -i1 file_I1.fq.gz \
    -bc_ltr ltr_bc.fa
```

A description of all arguments accepted by intmap_demux are provided below.

### intmap

Once the multiplexed FASTQ files have been demultiplexed, integration sites can be mapped with intmap. To reiterate, <b>intmap expects input files to be demultiplexed</b>. Below is a minimal<sup>†</sup> command to run intmap. In addition to this core set of arguments, there are over 70 additional user-adjustable parameters that can help fine-tune the analysis to match most experimental circumstances and analytical considerations.

```
intmap \
    -ltr_reads /path/to/demux/sample1_R1.fq.gz \
    -linker_reads /path/to/demux/sample1_R2.fq.gz \
    -ltr TGTGGAAAATCTCTAGCA \
    -linker TCAGGCATCGACGTCTCAG \
    -nm sample1 \
    -g_idx_dir /path/to/host/genome/index/directory \
    -g_idx_name host_index_name \
    -v_idx_dir /path/to/virus/genome/index/directory \
    -v_idx_name virus_index_name
```

<small>† Not quite -- the linker sequence isn't strictly required, but should you be analyzing paired-end data without a linker sequence if one's available?</small>

Broadly speaking, intmap automates 4 distinct processes. These are: removal of non-host sequences from the ends of the provided reads (a.k.a, "cropping"), alignment of the cropped reads to the host genome, post-alignment read processing, and coordinate retrieval. The cropping step looks for the given LTR/integrant sequence in `-ltr_reads` and the given linker sequence in `-linker_reads` (when present). For paired-end and long-read data, a linker sequence should, in the vast majority of cases, be provided. Both the LTR/integrant and linker sequences should be given in the 5' to 3' orientation. Reverse complementation is handled internally, when necessary. For fixed-end integrants, the LTR/integrant and linker sequences are sequences ≥ 10 bp that define the expected termini of the known, non-genomic sequences used for priming during library preparation. For integrants with truncated terminal repeats (e.g., AAV integrants), a full-length representation of the ITR (or equivalent) should be given (e.g., D-A'-C-C'-B-B'-A). An alternate ITR/integrant orientation can be concurrently specified with the `-ltr2` argument. When sequence matches are found, the matched sequences are removed from the respective reads and the cropped reads are saved to new FASTQ files (*_cropped.fq.gz). Cropped reads are then aligned using either [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [minimap2](https://github.com/lh3/minimap2) for short-read and long-read data, respectively. For all data-types, alignment is a two-step process. Reads are first aligned to the provided virus genome/vector sequence to filter out artifactual reads. The remaining unaligned reads are then aligned to the provided host genome. Aligned reads are subsequently QC'd for alignment quality, sequence quality, fragment length, read orientation, mismatch/indel rates, etc. The QC parameters are highly tunable (see the argument descriptions below). Mapped fragments are categorized as either uniquely- or multi-mapping during QC.

Following QC, the remaining reads are deduplicated. For uniquely-mapping fragments deduplication proceeds by collapsing fragments with identical alignment positions and UMIs (when present) into a single representative read. The representative read is chosen based on MAPQ and average Phred scores. Next, "fuzzy" deduplication groups fragments whose start/end coordinates are within `-len_diff` distance of the start/end coordinates of a nearby fragment. Within each fuzzy grouping, fragments are further deduplicated by UMI (when present) using an approach similar to the "directional" algorithm implemented in [UMI-tools](https://github.com/CGATOxford/UMI-tools). Multi-mapping reads are deduplicated by sequence similarity. For paired-end data, read pairs are concatenated if they are < `-seq_sim` * 100% similar. Otherwise, the `-ltr_read` read is used. MinHash fingerprints are then generated for each sequence. These MinHash fingerprints are then binarized and quickly compared using Faiss's binary index search. The Hamming distances between fingerprints are used to quickly estimate the Jaccard indexes according to $J \approx 1 - \frac{2d_H}{L}$, where $d_H$ is the Hamming distance and $L$ is the fingerprint length. Multi-mapping fragment sequences are then clustered according the specified `-hash_similarity` threshold. By default, this threshold is internally defined as $\frac{s^t}{2 - s^t} - 0.05$, where $s$ is `-seq_sim` and $t$ is `-token_size`. Coarse clusters are then sub-clustered into groups with sequence similarities >= `-seq_sim`, as evaluated using normalized Levenshtein distance. These clusters are then deduplicated by UMI analogously to groups of uniquely-mapping fragments.

When the `-reassign_mm` flag is set, intmap will attempt to reassign multi-mapping fragments according to sequence similarity with uniquely-mapping fragments or other multi-mapping fragments. This functionality improves, say, clonal abundance estimates when a clonal population is partially or completely composed of multi-mapping fragments.

After deduplication, intmap makes a final pass through the mapped fragments. The primary purpose of the final pass is to clean up the mapped sites. To do this, intmap collapses ISs within `-cluster_win` distance of an "abundant" IS (defined with the `-min_count` parameter). If two abundant ISs are nearby each other, the less abundant IS will be collapsed into the more abundant IS if the more abundant IS is > `-abundant_fc` more abundant. Low-confidence ISs are also removed in this step. Low-confidence ISs are defined as ISs sequenced < `-low_confidence_fc`-fold of the number of times the average IS was sequenced in the dataset.

When sequencing from the U3-end of integrants with a target site duplication (TSD), the `--U3` flag should be set and the `-tsd` argument should be defined as the expected size of the TSD. For HIV-1, the expected TSD width is 5 nucleotides. Setting these parameters is required to properly map the first position of the TSD. 

A final note to add is that several intmap arguments can be effectively turned off by setting them to certain values. For example, fuzzy deduplication can be turned off by setting `-len_diff` = 0 (for uniquely-mapping reads) and `-seq_sim` = 1 (for multi-mapping reads). Similarly, other parameters such as `-frag_ratio`, `-min_count`, and `-abundant_fc` can be rendered effectively irrelevant for a given analysis by setting their values to a high number (like 1000).

### intmap_multi
Experiments involving IS mapping usually involve more than one sample. intmap_multi automates the analysis of multiple datasets in a single command. The most common call to intmap_multi will look something like this:

```
intmap_multi \
    -s setup_file.txt \
    -n target_params > output_file.txt
```

In this example, `-s` is the setup file that tells the program what analysis parameters to use for each dataset and `-n` is the base name of an internally generated JSON file (the .json extension is added automatically). In this example, the software would output a file named 'target_params.json' that holds the information for each individual analysis being performed. The structure of the setup file was chosen to maximize flexibility. There are three requirements for the setup file. These are: line 1 must define the number of analyses being performed, line 2 must define what arguments are being passed to intmap, and line 3 on must define each of the arguments given in line 2 for each sample being analyzed. If a parameter is the same for all analyses (say, the same `-ltr` parameter for every analysis), that parameter only has to be defined once. If a parameter changes for any of the analyses, that parameter must be defined for each sample (i.e., the number of times indicated on line 1). The example setup file given below should help clarify these points.

```
3
ltr_reads, linker_reads, ltr, linker, nm, nthr, g_idx_dir, g_idx_name, v_idx_dir, v_idx_name
demux/sample1_R1.fq.gz, demux/sample2_R1.fq.gz, demux/sample3_R1.fq.gz
demux/sample1_R2.fq.gz, demux/sample2_R2.fq.gz, demux/sample3_R2.fq.gz
TGTGGAAAATCTCTAGCA
CAGGCATCGACGTCTCAG
sample1, sample2, sample3
4
/path/to/host/genome/index/directory
host_index_name
/path/to/virus/genome/index/directory
virus_index_name
```

In this example, 3 analyses will be performed (line 1 is 3) and only the required arguments are given (compare the arguments listed on line 2 to the minimal intmap example above). For each analysis, the parameters for `-ltr`, `-linker`, `-nthr`, `-g_idx_dir`, `-g_idx_name`, `-v_idx_dir`, and `-v_idx_name` are the same, so those values are only listed once. The arguments for `-ltr_reads`, `-linker_reads`, and `nm`, on the other hand, change across analyses. As such, those parameters are each defined individually (one parameter for each analysis).

For Boolean arguments/flags, intmap_multi expects those to be explicitly enumerated as True/False when given as an argument on line 2. For example, if I wanted to analyze a dataset with and without the inclusion of multi-mapping reads, the setup file would look something like:

```
2
ltr_reads, linker_reads, ltr, linker, nm, nthr, g_idx_dir, g_idx_name, v_idx_dir, v_idx_name, no_mm
demux/sample1_R1.fq.gz
demux/sample1_R2.fq.gz
TGTGGAAAATCTCTAGCA
CAGGCATCGACGTCTCAG
sample1, sample1_no_mm
4
/path/to/host/genome/index/directory
host_index_name
/path/to/virus/genome/index/directory
virus_index_name
False, True
```

Finally, intmap accepts 3 arguments that can be defined with more than a single value. These are `-contam`, `-remove_chr`, and `-annotations`. To pass multiple values to these arguments in the setup file, separate the intra-analysis values with semi-colons and the inter-analysis values with commas. As an example, assume that we are analyzing two IS datasets, that each dataset has two distinct contaminating sequences, and that the contaminating sequence differ between datasets. The setup file would look something like:

```
2
ltr_reads, linker_reads, ltr, linker, nm, nthr, g_idx_dir, g_idx_name, v_idx_dir, v_idx_name, contam
demux/sample1_R1.fq.gz, demux/sample2_R1.fq.gz
demux/sample1_R2.fq.gz, demux/sample2_R2.fq.gz
TGTGGAAAATCTCTAGCA
CAGGCATCGACGTCTCAG
sample1, sample2
4
/path/to/host/genome/index/directory
host_index_name
/path/to/virus/genome/index/directory
virus_index_name
ATTGCATCAA; TTTGCGGAGC, GCGATATAGC; ATGGACTACT
```

### Arguments
#### intmap
1. `-ltr_reads`: The path to integrant-end reads. Expects the defined integrant sequence to be near the 5' end of the read.
2. `-linker_reads`: The path to linker-end or non-integrant-end reads (i.e., the mate to ltr_reads). Expects the defined linker sequence near the 5' end of the read. Only used with paired-end data.
3. `-ltr`: The LTR/integrant sequence to match. Must be given in 5' to 3' orientation. Must be at least 10 nt long.
4. `-ltr2`: An alternative orientation of the LTR/integrant search sequence. Used for e.g., distinguishing between flip and flop orientations of AAV ITRs.
5. `-linker`: The linker sequence to match. Required for paired-end and long-read data. Must be given in the 5' to 3' orientation, even when sequenced in the 3' to 5' direction on the 3' end of ltr_reads (e.g., in long-read data). Must be at least 10 nt long.
6. `-ltr_error_rate`: The acceptable error rate in the LTR sequence. Defaults to 0.1. Cannot exceed 0.4.
7. `-linker_error_rate`: The acceptable error rate in the linker sequence. Defaults to 0.1. Cannot exceed 0.4.
8. `-contam`: Possible contaminating sequences immediately downstream of the defined integrant search sequence.
9. `-ltr_primer`: The 3' end of the integrant primer used for PCR amplification. Used to control for mispriming during fragment amplification. Must be at least 10 bp.
10. `-U3` or `--U3`: Whether or not the integrant-host junctions are sequenced from the U3 end of the provirus (or equivalent). Used for defining the 1st bp of the target site duplication, defining integrant orientation, etc.
11. `-nthr`: The number of threads to use for processing. Defaults to 1.
12. `-nm`: The name of the dataset. Used for naming output files.
13. `-ttr`, `--ttr`, `-truncated_terminal_repeat`, or `--truncated_terminal_repeat`: Tells the software to look truncation of the defined integrant search sequence. When set/True, the first min_ttr_len bases of the search sequence are assumed to be obligatory (an anchor sequence).
14. `-min_ttr_len`: The minimum allowed length of the truncated terminal repeat. Only used when ttr is True. Defaults to 10.
15. `-cap_ttr_error` or `--cap_ttr_error`: When set/True, the number of mismatches allowed for truncated integrant ends is capped at min_ttr_len * ltr_error_rate.
16. `-leftmost` or `--leftmost`: When set/True, cutadapt looks for the leftmost occurrence of the defined LTR/integrant search sequence. When not set/False (default), the rightmost occurrence is found. Does not affect linker search sequence position.
17. `-mixed` or `--mixed`: When set/True, the cropping step accommodates situations where library preparation can result in target LTR-ends on either terminus of a given fragment. Cannot be used with ttr or with single-end short-read data.
18. `-min_frag_len`: The minimum allowed fragment length. Defaults to 25.
19. `-max_frag_len`: The maximum allowed fragment length. Defaults to 2000.
20. `-g_idx_dir`: The path to the Bowtie2 or minimap2 genome index directory. Note that this should not contain the index prefix.
21. `-g_idx_name`: The genome index prefix (e.g., hs1).
22. `-v_idx_dir`: The path to the virus genome index directory. Note that this should not contain the index prefix.
23. `-v_idx_name`: The virus genome index prefix.
24. `-aln_mismatch_rate`: The allowable mismatch rate in read alignment. N's and soft-clipped bases are counted as mismatches. Defaults to 0.10.
25. `-aln_indel_rate`: The allowable indel rate in read alignment. Defaults to 0.02.
26. `-min_mapq`: The minimum acceptable MAPQ score. Does not apply to multimapping reads. Defaults to 10.
27. `-match_after`: The number of mismatched bases allowed at the start of the integrant-end read alignment. Defaults to 2.
28. `-tsd`: The size of the target site duplication. Defaults to 1. Especially important when sequencing from the U3-end of retroviruses. Must be >= 1.
29. `-disorg` or `--disorg`: When set/True, files are not moved to type-specific subdirectories.
30. `-ltr_umi_offset`: The number of basepairs between the first basepair of the provided LTR/integrant sequence and the last basepair of the UMI. Defaults to 0 (i.e., the UMI immediately precedes the given LTR sequence).
31. `-ltr_umi_len`: The length of the LTR/integrant-end UMI. Defaults to 0 (no UMI present).
32. `-ltr_umi_pattern`: Any pattern that defines the LTR/integrant-end UMI. For example, 4 random bases, followed by 4 known bases, followed by 4 random bases should be given as NNNNATGCNNNN, where ATGC is replaced with the known bases in the UMI sequence. Reads whose UMI does not match the given pattern are removed.
33. `-linker_umi_offset`: The number of basepairs between the first basepair of the provided linker sequence to the last basepair of the UMI. Defaults to 0 (i.e., the UMI immediately precedes the given linker sequence).
34. `-linker_umi_len`: The length of the linker-end UMI. Defaults to 0 (no UMI present).
35. `-linker_umi_pattern`: Any pattern that defines the linker UMI. For example, 4 random bases, followed by 4 known bases, followed by 4 random bases should be given as NNNNATGCNNNN, where ATGC is replaced with the known bases in the UMI sequence. Reads whose UMI does not match the given pattern are removed.
36. `-no_mm` or `--no_mm`: When set/True, multimapping reads are discarded.
37. `-min_qual`: The minimum acceptable average base quality score across a read. Used to filter out reads with low-quality base calls. Defaults to 10.
38. `-len_diff`: The maximum allowed length difference on either end of uniquely mapping fragments for fragments to be grouped together during duplicate removal. Alternatively, the maximum absolute length difference for multimapping fragments. Defaults to 5.
39. `-cluster_win`: The size of the window used to cluster nearby sites in the final pass. Defaults to 5.
40. `-frag_ratio`: The ratio of fragment counts used for UMI comparisons. A match is called if the higher number of counts > (frag_ratio \* lower number). Defaults to 2.
41. `-genome_fasta`: The file path to the genome fasta file. The file must be either bgzip compressed or uncompressed. Used when ltr_primer is given to look for mispriming.
42. `-genome_fasta_index`: The file path to the genome fasta index file. This is not strictly necessary if the genome fasta index file is named \<genome_fasta\>.fai for an uncompressed genome_fasta file, or \<genome_fasta\>.gzi for a bgzip compressed genome_fasta file.
43. `-check_mispriming`: When set/True, look for mispriming around each mapped IS. Aligment parameters are currently set internally.
44. `-num_perm`: The number of permutations in the minhash fingerprint. Defaults to 128. Must be divisible by 8 and >= 32.
45. `-seq_sim`: Sequence similarity threshold. Used for validating multimapping read groups and, in certain instances, assigning multimapping read positions. Defaults to 0.95. Must be > 0 and <= 1.
46. `-token_size`: The token size for MinHash generation. Defaults to 4 for short-read data and 8 for long-read data. A value of -1 tells the program to use default values.
47. `-hash_similarity`: Defines the similarity threshold between MinHashes. Used for coarse clustering of multimapping fragments. Must be >= 0 and <= 1. By default, internally defined as seq_sim^token_size / (2 - seq_sim^token_size) - 0.05.
48. `-b_bits`: Defines the number of bits with which to represent each MinHash. Must be between 1 and 32. Defaults to 32.
49. `-minimizer_win`: Defines the minimizer window size. By default, internally defined as 2 \* token_size for long-read sequences and token_size for short-read sequences. Must be > 0 and <= min_frag_len - token_size + 1. When -1, all k-mers are included in the MinHash.
50. `-min_count`: The minimum number of mapped sites required to define a particular integration site as abundant. Used during the final pass to clean up mapped ISs. Defaults to 1.
51. `-abundant_fc`: The fold-change threshold between nearby abundant sites for them to be considered the same. If two nearby abundant sites have abundance fold-change values < abundant_fc, they are considered different. Defaults to 2.
52. `-low_confidence_fc`: The fold-change threshold relative to the average number of observations per fragment for a mapped site to be considered low-confidence. If a site has an observed value < low_confidence_fc-fold from the average number of observations per fragment, the site is removed. Defaults to 100.
53. `-reassign_mm` or `--reassign_mm`: When set/True, attempts to reassign multimapping reads based on concordance with uniquely mapping reads or similar multimapping reads.
54. `-mm_k`: The size of the k-mers used to reassign multimapping reads. Defaults to 15. Must be <= min_frag_len - len_diff.
55. `-mm_group_threshold`: Defines the lower multimapping group size threshold for whole-group reassignment to a probable location. Expressed as a fraction of the total number of multimapping reads. Defaults to 0.002.
56. `-um_reassign_diff`: The fuzzy similarity difference enforced when choosing the best uniquely-mapping fragment position for multi-mapping fragment reassignment. Used in conjection with um_reassign_fc. Must be >= 0 and < 0.25. Defaults to 0.05.
57. `-um_reassign_fc`: Defines the count fold-change required between two uniquely-mapping fragments with similar sequence similarity scores to a multi-mapping fragment in order to re-define the best match. Used in conjunction with um_reassign_diff. Must be >= 1. Defaults to 2.
58. `-no_dedup` or `--no_dedup`: When set/True, the pipeline stops after alignment and QC.
59. `-check_consensus` or `--check_consensus`: When set/True, reads from input files are sampled to identify common variants of the given LTR and linker sequences. Results are written to stdout.
60. `-consensus_sample_size`: The number of reads to sample for consensus checking. Defaults to 10000.
61. `-consensus_threshold`: The fractional threshold for sequence significance. Used during consensus checking. If the frequency of any sequence is >= this value, that sequence will be reported. Defaults to 0.01.
62. `-consensus_error_rate`: The error rate for consensus searching. Any sequence with <= len(LTR/linker sequence) \* consensus_error_rate errors will be reported as a variant provided sufficient representation (see consensus_threshold). Defaults to 0.3.
63. `-extract_consensus` or `--extract_consensus`: When set/True, a consensus sequence of length consensus_length is determined from provided FASTQ files. Useful for determining integrant/linker sequences when they are unknown/ambiguous.
64. `-consensus_length`: The length of the consensus sequence to generate. Used with extract_consensus. Defaults to 100.
65. `-consensus_base_threshold`: The threshold required to call a base the consensus. Positions with base frequency values below this are returned as N. Defaults to 0.65.
66. `-consensus_cpt_penalty`: The penalty term for changepoint estimation. Used in defining the estimated consensus sequence. Higher values increase sensitivity, but also potentially increase noise. Defaults to 1.5.
67. `-consensus_peak_threshold`: The fractional threshold relative to the maximum segment mean for a changepoint segment to be considered part of the consensus. Defaults to 0.9.
68. `-umi_dedup_fp`, `--umi_dedup_fp`: When set/True, checks for conflicting UMIs during the final pass. These would be similar UMIs on fragments mapping to different genomic positions.
69. `-mle` or `--mle`: When set/True, site abundance is estimated using maximum likelihood estimation. This cannot be set when either linker_umi_len or ltr_umi_len > 0 or when single-end short-read sequencing is being analyzed. Requires random genome fragmentation.
70. `-long_read`, `--long_read`, `-lr`, or `--lr`: When set/True, reads are processed as long-reads (ONT/PB).
71. `-lr_type`: The type of long-read sequencing. Fills in the -x value in minimap2. Choices: map-ont, map-pb, map-hifi, lr:hq. Defaults to map-ont.
72. `-annotations`: File paths to genomic annotation files (comma-separated). When not None, ISs overlapping each feature-set and the distance from each IS to the nearest feature in each feature-set are determined. Annotation files are assumed to be in BED format. Nearest distance is determined independently of feature strandedness.
73. `-write_peaks` or `--write_peaks`: When set/True, mapped sites are collapsed into peaks and peak ranges our output.
74. `-peak_win`: Maximum distance (bp) between sites to merge into a single peak. Only used when write_peaks is True. Defaults to 25.
75. `-peak_alpha`: Significance threshold for retaining single-strand IS positions in peaks. Expressed as the probability of only observing mapped sites on a given strand. Sites where this probability is >= peak_alpha are deemed high quality. Only used when write_peaks is True. Defaults to 0.05.
76. `-remove_chr`: The names of chromosomes to remove during mapping.
77. `-bt2_path`: The path to the Bowtie2 executable. Defaults to 'bowtie2'. OK to ignore for standard installations.
78. `-sam_path`: The path to the samtools executable. Defaults to 'samtools'. OK to ignore for standard installations.
79. `-cut_path`: The path to the cutadapt executable. Defaults to 'cutadapt'. OK to ignore for standard installations.
80. `-seqtk_path`: The path to the seqtk executable. Defaults to 'seqtk'. OK to ignore for standard installations.
81. `-bed_path`: The path to the bedtools executable. Defaults to 'bedtools'. OK to ignore for standard installations.
82. `-mm2_path`: The path to the minimap2 executable. Defaults to 'minimap2'. OK to ignore for standard installations.

#### intmap_demux
1. `-ltr_reads`: Integrant-end FASTQ file.
2. `-linker_reads`: Linker-end/non-integrant-end FASTQ file.
3. `-bc_ltr`: The FASTA file containing sample-index pairs for ltr_read in-line barcodes or I1 reads.
4. `-bc_linker`: The FASTA file containing sample-index pairs for linker_read in-line barcodes or I2 reads.
5. `-i1`: The I1 FASTQ file, if applicable.
6. `-i2`: The I2 FASTQ file, if applicable.
7. `-nthr`: The number of threads to use for demultiplexing. Defaults to 1.
8. `-error_rate`: The acceptable error rate in index matching. Defaults to 0.1.
9. `-keep_all` or `--keep_all`: Tells the program to keep all combinatorial matches. Only applies to in-line barcodes.
10. `-cutadapt_path`: The path to the cutadapt executable. Only potentially need for non-standard installations.
11. `-ltr_only` or `--ltr_only`: Tells the program to only use ltr_reads for demultiplexing.
12. `-linker_only` or `--linker_only`: Tells the program to only use linker_reads for demultiplexing.
13. `-header_index` or `--header_index`: Tells the program to look for reads in FASTQ headers.
14. `-three_prime` or `--three_prime`: Tells the program to only look for barcodes on the 3' end of ltr_reads.
15. `-no_index` or `--no_index`: Tells cutadapt not to build barcode indexes for faster searching.
16. `-file_match`: A two-column TSV or CSV file used for renaming files.

#### intmap_multi
1. `-s` or `-setup_file`: The setup file. A text file holding the arguments and parameters for each intended run.
2. `-n` or `-json_name`: The base name of the JSON file created by the software. This JSON file is fed back into the software for analysis. The .json extension is appended automatically.
3. `-a` or `-args_json`: An arguments JSON typically created internally. Only given if an arguments JSON storing the parameters of interest already exists.
4. `-arg_sets_only` or `--arg_sets_only`: Tells the program to stop after creating the arguments JSON.

### Helper scripts
#### intwrap
A bash script -- intwrap.sh -- is provided to wrap demultiplexing with intmap_demux and serialized analyses with intmap_multi into a single command. The script is copied to `~/CONDA_PREFIX/intmap/bin/` when intmap is installed via the install script, making intwrap.sh callable from the command line. The script takes as arguments the multiplexed sequencing data, barcode/index FASTA file(s), sample index files (if applciable), the intmap_multi setup file, the name of the argument JSON file (the actual file is created internally), and the name of the output file to store the analysis information. Other arguments are also accepted. Run `intwrap -h` for more information. Below is an example of using intwrap for paired-end data:

```
intwrap \
    -r1 multiplexed_R1.fq.gz \
    -r2 multiplexed_R2.fq.gz \
    -bc_r1 ltr_bc.fa \
    -bc_r2 linker_bc.fa \
    -setup example_setup.txt \
    -json example_args \
    -out example_out.txt
```

#### prepim
prepim.sh is another helper script useful for creating setup files for intmap_multi. When many samples are being analyzed, for example, it can be cumbersome to manually list every file name, every sample name, every linker sequence (if they are different), etc. It can also be cumbersome to create the requisite FASTA files for demultiplexing. prepim is intended to alleviate these annoyances. prepim takes as input a tab-delimited text file with at least three columns: column 1 stores sample names, column 2 stores the sequenced integrant-end sequence, and column 3 stores the sequenced linker-end sequence. By default, prepim expects the integrant and linker sequences to be everything that follows the appended sequencing adapters (for in-line barcodes). If in-line barcodes are not used in the experiment, the `--no_bc` flag can be set. If other index sequences are used, the `-index_col1` and `-index_col2` arguments should be set to values > 0, defining the column numbers of the respective index sequences. If the `--no_bc` flag is set and `-index_col1/2` arguments are not given, no FASTA files are written. Besides the FASTA files, prepim writes to STDOUT simple comma-separated vectors of sample names, sample read file names, integrant-end sequences, and linker-end sequences. Integrant- and linker-end sequences are derived from the right-most (3') ends of the provided sequences. For single-end data (i.e., any data without paired reads), the `--se` flag should be set. The prepim output makes it possible to easily copy/paste the string of values required in the intmap_multi setup file. Like intwrap, prepim is copied to `~/CONDA_PREFIX/intmap/bin/` upon intmap installation with the provided install script, making it callable from the command line. A basic prepim example is as follows:

```
prepim input.tsv
```

## Licensing
intmap is licensed under the GPL-3 License (see LICENSE). It additionally relies on several third-party tools, which are governed by their own licenses. See `THIRD_PARTY_LICENSES.md` for details.