# Append LTR-end barcodes to the LTR-end FASTQ file
# Outputs Undetermined_S0_L001_R2_001_bc_append.fastq.gz
intmap_demux \
    -r1 Undetermined_S0_L001_R2_001.fastq.gz \
    -r2 Undetermined_S0_L001_I1_001.fastq.gz \
    --append_index

# Generate clone-specific FASTQ files
# Uses linker-end barcodes on R1
intmap_demux \
    -r1 Undetermined_S0_L001_R1_001.fastq.gz \
    -r2 Undetermined_S0_L001_R2_001_bc_append.fastq.gz \
    -bc_r1 linker_bc.fa \
    --r1_only

# Generate replicate-specific FASTQ files
# Uses LTR-end barcodes on R2
# LTR-end barcodes for each clone are stored in separate files
for f in demux/*_R1.fq.gz; do
    filename=$(basename "$f" _R1.fq.gz)
    intmap_demux \
        -r1 "demux/${filename}_R1.fq.gz" \
        -r2 "demux/${filename}_R2.fq.gz" \
        -bc_r2 "${filename}_ltr_bc.fa" \
        --r2_only
done

# Run intmap
# The setup file includes all of the parameters needed for each sample
# Note that the INSPIIRED workflow sequences the LTR-end in R2 and the linker-end in R1
# intmap expects the LTR-end to be R1
# intmap can be "tricked" by simply giving actual R2 reads as R1 and vice versa
intmap_multi \
    -s inspiired_setup.txt \
    -n inspiired_params > inspiired_out.txt

# If inspiired_params.json already exists, run:
intmap-multi \
    -s inspiired_params.json > inspiired_out.txt