# Changelog

## v1.0.1 (2025-11-17)
+ Updated fragment choosing logic to choose fragments based on MAPQ and average Phred score.

## v2.0.0 (2026-04-10)
+ Changed r1/r2 arguments in the intmap_demux and intmap modules to ltr_reads and linker_reads, respectively. This is intended to clarify which read file should be specified with each argument.
+ Changed read file output formats from *_R1/R2.fq.gz to *_Rltr/Rlink.fq.gz in both intmap_demux and the cropping step of intmap.
+ Refactored the intmap module to handle single-end, paired-end, and long-read data. Accordingly, eliminated the intmap_se module.
  + Single-end data is automatically inferred by specification of ltr_reads but not linker_reads.
  + Long-read data should be specified by setting the `long_read` flag.
+ Modified the cropping routine to use cutadapt instead of fuzzy regex matching for identifying/removing integrant and linker sequences.
  + LTR and linker sequences are now given in a single sequence instead of split across 5' and 3' chunks.
+ Added the `leftmost` flag to tell cutadapt to find the leftmost occurrence of the given LTR/integrant-end search sequence during cropping.
+ Added the `low_confidence_fc` parameter.
  + Expressed as a fold-change from the average number of times individual fragments are <i>sequenced</i> in the given dataset.
  + Fragments sequenced less-than `low_confidence_fc`-fold of the average are deemed low-confidence and removed.
+ Changed MinHash generation from the purely Pythonic datasketch to the faster mmh3.
+ Added the option to change MinHash b-bit value.
  + Must be between 1 and 32.
  + Still defaults to the full 32-bit MinHash values.
+ Updated token_size defaults to be 4 for short-read data and 8 for long-read data.
+ Changed default hash_similarity value to be seq_sim^token_size / (2 - seq_sim^token_size) - 0.05.
+ Incorporated a minimizer for MinHash generation. Default win_sizes are token_size for short-read data and 2 * token_size for long-read data.
+ Changed `count_fc` argument name to `abundant_fc`.
+ Corrected the Jaccard index calculation for multimapping fragments to $J \approx 1 - \frac{2d_H}{L}$.
+ Modified `ranged_groupby()`:
  + Replaced the Pythonic grouping approach with a vectorized NumPy approach.
  + Added a post-grouping graph-based merging step to look for possible missed connections.
+ Optimized UMI deduplication algorithm.
  + Replaced union-find with a two-phase approach that resolves connected components through iterative BFS from the highest count UMI.
  + Reduced memory overhead by searching for 1-edit variants only among observed UMIs instead of all possible 1-edit variants.
+ Updated mm_coarse_cluster() to choose the coarse group representative based on both length and highest Jaccard index.
+ Updated consensus generating functions to return:
  + A raw consensus based on the defined fractional cutoff.
  + An estimated consensus based on change point estimation of the fractional abundance values across the entire consensus length.
+ Updated consensus checking functions:
  + cutadapt is now used instead of regex matching to find matching consensus sequences (similar to cropping).
  + Overrepresented sequences immediately adjacent to the consensus sequence are automatically returned.
+ All output coordinate files are now output in BEDn+ format (https://github.com/samtools/hts-specs/blob/master/BEDv1.pdf).
+ Added the `annotations` option to define genomic features with which to compare IS targeting.
  + Fractional IS overlap with each genomic feature-set is printed in the intmap output.
  + Files denoting overlapping features and distance to the nearest feature are output when annotation files are given.
+ Added the `remove_chr` argument to remove reads/pairs aligning to unwanted chromosomes (e.g., 'chrM').
+ Replaced fuzzy regex mispriming analysis with Smith-Waterman alignment.
+ Added the `cap_ttr_error` option to cap truncated integrant end fragment errors to min_ttr_len * ltr_error_rate errors, regardless of fragment length.
+ Added the `write_peaks` flag for calling hotspots of AAV integration into Cas-mediated DSBs.
  + Quality ISs are identified by the presence of integration in both orientations.
  + Sites lacking integrants in both orientations can be deemed quality sites if it is statistically probable to only see integration in one orientation.
  + Sites within `peak_win` distance of each other are rolled up into a single peak. 
  + Peak boundaries and the peak apex (the position within the peak with the most ISs) are returned.
+ Added the `-um_reassign_diff` and `-um_reassign_fc` arguments that can modify the behavior of which uniquely-mapping fragment is deemed the best match for a multi-mapping fragment. The old hard-coded behavior is matched with `um_reassign_diff` = 0 and `um_reassign_fc` = 1.
+ Removed `ltr_cufp` and `linker_cufp` flags in favor of a combined `umi_dedup_fp` flag.
+ Updated unit tests.

## v2.0.1 (2026-04-27)
+ Fixed a bug that could cause cropping to fail for single-end data.
+ Expanded `--write_peaks` functionality to non-AAV-like systems.
  + Includes a `--stranded_peaks` flag to call peaks by strand.
+ Expanded `--mixed` functionality to truncated terminal repeats.
+ Added more arguments to the initial print statement.
+ Added a dynamic seq_sim adjustment for k-mer frequency comparisons.
  + Replaces the static and universal seq_sim - 0.05 adjustment.
  + Now calculated as seq_sim - max(min(0.2, k / (min(len_i, len_j) - k + 1)), 0.05)