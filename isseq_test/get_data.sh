#$ -S /bin/bash
#$ -N get_data
#$ -pe pvm 6
#$ -l m_mem_free=24G
#$ -wd /home/user/isseq_test
#$ -j n
#$ -m bea

# sra environment created as follows:
# conda create -n sra sra-tools python>=3.12 
source "/path/to/miniforge3/etc/profile.d/conda.sh"
conda activate sra

python get_fastq.py \
    -a accessions.txt \
    -o data \
    -nthr 30

export -f process_file

process_file() {
    f=$1
    filename=$(basename "$f" .fastq.gz)
    zcat "$f" | awk '{if (NR % 4 == 3) {$0 = "+"} print}' | gzip > "data/${filename}_clean.fastq.gz"
    rm $f
}

export -f process_file

find data -name "*.fastq.gz" | parallel -j 6 process_file

# modified from https://unix.stackexchange.com/a/646449
while read -r accession sample; do
    for name in data/"$accession".*; do
        [ -e "$name" ] || continue
        experiment=data/"$sample".fastq.gz
        mv -- "$name" "$experiment"
    done
done < file_match.txt

conda deactivate
