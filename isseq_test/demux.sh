intmap-demux \
    -r1 data/single_cell_clones_R1.fastq.gz \
    -r2 data/single_cell_clones_R2.fastq.gz \
    -bc_r1 scc_ltr.fa \
    -bc_r2 scc_linker.fa

intmap-demux \
    -r1 data/serial_dil_250_500_set1_R1.fastq.gz \
    -r2 data/serial_dil_250_500_set1_R2.fastq.gz \
    -bc_r1 dil1_ltr.fa \
    -bc_r2 dil1_linker.fa \
    -file_match dil1_file_match.txt

intmap-demux \
    -r1 data/serial_dil_250_set2_R1.fastq.gz \
    -r2 data/serial_dil_250_set2_R2.fastq.gz \
    -bc_r1 dil2_ltr.fa \
    -bc_r2 dil2_linker.fa \
    -file_match dil2_file_match.txt

intmap \
    -r1 demux/CL.6_R1.fq.gz \
    -r2 demux/CL.6_R2.fq.gz \
    -ltr5 AGTCAGTGTGGA \
    -ltr3 AAATCTCTAGCA \
    -linker5 AGATCTGGAA \
    -linker3 TGAACTGGCC \
    -ltr1_primer TTGAGTGCTTCA \
    -c TGGCGCCCGAAC \
    -nm CL.6 \
    -nthr 4 \
    -bt2_idx_dir /Users/gbedwell/Documents/github/T2T_genome/indexes/bowtie2 \
    -bt2_idx_name hs1 \
    -v_idx_dir /Users/gbedwell/Desktop/isseq_test/pRRLSIN \
    -v_idx_name pRRLSIN \
    -genome_fasta /Users/gbedwell/Documents/github/T2T_genome/UCSC/hs1.fa.bgz \
    -linker_umi_offset 0 \
    -linker_umi_len 18 \
    -no_mm


















intmap-demux \
    -r1 data/serial_dil_250_set2_R1.fastq.gz \
    -r2 data/serial_dil_250_set2_R2.fastq.gz \
    -bc_r1 scc_ltr.fa \
    -bc_r2 scc_linker.fa


