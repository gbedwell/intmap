# intmap

intmap is a python package for mapping the locations of retroviral genomic integration from paired-end Illumina sequencing data.

## Installation
The easiest way to install intmap is by running <code>bash install.sh</code>. This bash script will:

1. Create an intmap conda environment,
2. install all (version-controlled) intmap dependencies, and
3. install the intmap package.

<code>install.sh</code> requires an initialized conda instance. If conda is not available on your computer, you can easily install conda by following the recommendations [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). <code>install.sh</code> can additionally install new versions of intmap seamlessly from either the base environment or an activated intmap environment.

## Package Structure
There are 3 distinct intmap modules. These are: intmap-demux, intmap, and intmap-multi. I will explain each of these in more detail below.

### intmap-demux
Currently, intmap-demux is written to demultiplex multiplexed FASTQ files based on unique LTR- and linker-end sequences. Using sequenced LTR- and linker-end sequences helps to catch instances of index hopping that can occur during Illumina sequencing. Behind the scenes, intmap-demux uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/#) to look for target sequences and allocate them to the appropriate file.

Given multiplexed paired-end FASTQ files file_R2.fq.gz and file_R2.fq.gz, intmap-demux can be run like so:

```
intmap-demux \
    -r1 file_R1.fq.gz \
    -r2 file_R2.fq.gz \
    -bc_r1 ltr_bc.fa \
    -bc_r2 linker_bc.fa
```

The fasta files ltr_bc.fa and linker_bc.fa contain the matched LTR- and linker-end sequences that the program looks for. Note that all of the fasta headers must be unique. By default, intmap-demux allows for a substitution rate of 0.1. This value can be changed with the <code>-e</code> argument. When cutadapt is finished running, intmap-demux will place the proper R1 and R2 files into a subdirectory called 'demux' within the current working directory. These files can then be used directly with intmap. Note that multiplexed FASTQ files do not have to be demultiplexed with intmap-demux to be compatible with intmap. Any demultiplexing procedure will suffice, provided the LTR- and linker-end sequences immediately adjacent to the provirus-host genome junction are retained after demultiplexing. Usually, this is not a problem, as demultiplexing typically only looks for, and optionally removes, short sequences at the very 5' ends of reads.

intmap-demux can be used for demultiplexing multiplexed reads with barcode structures that deviate from the simplest case of unique, sample-specific barcodes on LTR- and/or linker-ends. See the inspiired_test and isseq_test directories in this repository for examples.

### intmap

Once the original FASTQ files have been demultiplexed (however you prefer to do that), integration sites can be mapped with intmap. Below is a (near-)minimal command to run intmap. I say <b>near</b>-minimal because one of the given arguments, <code>-nthr</code>, is not <i>strictly</i> required, but given the hardware available on modern computers and HPC clusters, I can't really think of a reason to use only one core. All of the other arguments are required. These include: file paths, sample-specific sequences, virus-specific sequences, genome information, etc. There are a myriad of other arguments available to the user that I will not go over here. See the help documentation for more information.

```
intmap \
    -r1 /path/to/demux/sample1_R1.fq.gz \
    -r2 /path/to/demux/sample1_R2.fq.gz \
    -ltr5 ccttttagtcag \
    -ltr3 TGTGGAAAATCTCTAGCA \
    -linker5 TACATCTCAG \
    -linker3 GCATCGACGTCTCAG \
    -ltr1_primer gagatccctc \
    -nm sample1 \
    -nthr 4 \
    -bt2_idx_dir /path/to/host/genome/index/directory \
    -bt2_idx_name host_index_name \
    -v_idx_dir /path/to/virus/genome/index/directory \
    -v_idx_name virus_index_name \
    -genome_fasta /path/to/host_genome.fa.bgz
```

Broadly speaking, intmap automates 4 distinct processes. These are: removal of non-host sequences from the ends of both R1 and R2 (referred to in the program as "cropping"), alignment of the cropped reads to the host genome, post-alignment read processing, and coordinate retrieval. The cropping step looks for the concatenated ltr5 and ltr3 sequences in R1 and the concatenated linker5 and linker3 sequences in R2. The identified sequences are removed from the respective reads and the reads are output to new *_cropped.fq.gz files. These new files are then given to [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment. Alignment is a two-step process. The reads are first aligned to the provided virus genome/vector sequence to filter out potentially unwanted reads. The remaining unaligned reads are then aligned to the provided host genome using the <code>--very_sensitive</code> end-to-end alignment mode. Aligned reads are then QC'd for sequence quality, read orientation, the number of mismatches/indels, etc. Reads passing QC are then deduplicated and the final set of fragments are checked for possible mispriming artifacts. Finally, the coordinates of the final mapped fragments and integration sites are returned.

### intmap-multi
In most situations, integration site mapping experiments are designed to include more than one sample. To streamline the mapping of integration sites from multiple samples, users can use intmap-multi to serially analyze a given set of samples in a single command. The parameters for the respective samples can be input into intmap-multi interactively, but it is more likely that setup file will be given to the program. The setup file should be structured in the following way:

```
3
r1, r2, ltr5, ltr3, linker5, linker3, ltr1_primer, nm, nthr, bt2_idx_dir, bt2_idx_name, v_idx_dir, v_idx_name, genome_fasta
demux/sample1_R1.fq.gz, demux/sample2_R1.fq.gz, demux/sample3_R1.fq.gz
demux/sample1_R2.fq.gz, demux/sample2_R2.fq.gz, demux/sample3_R2.fq.gz
AGTCAGTGTGGA
AAATCTCTAGCA
AGATCTGGAA
TGAACTGGCC
TTGAGTGCTTCA
TGGCGCCCGAAC
sample1, sample2, sample3
4
/path/to/host/genome/index/directory
host_index_name
/path/to/virus/genome/index/directory
virus_index_name
/path/to/host_genome.fa.bgz
```

The first line of the setup file is the number of samples being analyzed. The second line is a comma-delimited set of arguments that are given to intmap for each sample run. All of the subsequent lines must either be a single entry that the program will interpret as the same for each sample or a comma-delimited set of values equal in length to the number of samples defined in line 1. Boolean arguments must be explicitly enumerated as True or False, if used. This structure is intended to provide flexibility in the argument values used for each call to intmap.

Once the setup file is constructed, intmap-multi can be run as follows:

```
intmap-multi \
    -setup_file setup_file.txt \
    -json_name target_params
```

This command will create the file target_params.json and use the entries therein to run intmap for each enumerated sample. If the output JSON file already exists, intmap-multi will not overwrite it. The parameter JSON can be used in the intmap-multi command directly:

```
intmap-multi \
    -args_json target_params.json
```