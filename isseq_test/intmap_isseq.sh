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

# This will create isseq_params.json and then run it.
intmap-multi \
    -setup_file isseq_setup.txt \
    -json_name isseq_params

# If isseq_params.json already exists, run:
intmap-multi \
    -args_json isseq_params.json