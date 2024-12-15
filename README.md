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
Currently, intmap-demux is written to demultiplex multiplexed FASTQ files based on unique LTR- and linker-end sequences. Eventually, I hope to include functionality that allows for the serial demultiplexing of paired-end files that have been demultiplexed with e.g., Illumina's i5 and i7 indices. Using sequenced LTR- and linker-end sequences helps to catch instances of index hopping that can occur during Illumina sequencing. Behind the scenes, intmap-demux uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/#) to look for target sequences and allocate them to the appropriate file.

Given multiplexed paired-end FASTQ files file_R2.fq.gz and file_R2.fq.gz, intmap-demux can be run like so:

```
intmap-demux \
    -r1 file_R1.fq.gz \
    -r2 file_R2.fq.gz \
    -bc_ltr ltr_bc.fa \
    -bc_linker linker_bc.fa \
    -nthr 2
```

The fasta files ltr_bc.fa and linker_bc.fa contain the matched LTR- and linker-end sequences that the program looks for. Note that all of the fasta headers must be unique. By default, intmap-demux allows for a substitution rate of 0.1. This value can be changed with the <code>-e</code> argument. When cutadapt is finished running, intmap-demux will place the proper R1 and R2 files into a subdirectory called 'demux' within the current working directory. These files can then be used directly with intmap. Note that multiplexed FASTQ files do not have to be demultiplexed with intmap-demux to be compatible with intmap. Any demultiplexing procedure will suffice, provided the LTR- and linker-end sequences immediately adjacent to the provirus-host genome junction are retained after demultiplexing. Usually, this is not a problem, as demultiplexing typically only looks for, and optionally removes, short sequences at or near the very 5' ends of reads.

### intmap

Once the original FASTQ files have been demultiplexed (however you prefer to do that), integration sites can be mapped with intmap. Below is a (near-)minimal command to run intmap. I say <b>near</b>-minimal because one of the given arguments, <code>-nthr</code>, is not <i>strictly</i> required, but given the hardware available on modern laptops, desktops, and HPC clusters, I can't really think of a reason to use only one core. The other arguments are required. These include: file paths, sample-specific sequences, virus-specific sequences, genome information, etc. There are a myriad of other arguments available to the user that I will not go over here. See the help documentation for more information.

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

Broadly speaking, intmap automates 3 distinct processes. These are: removal of non-host sequences from the ends of both R1 and R2, alignment of the cropped reads to the host genome, and post-alignment processing.