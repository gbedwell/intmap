# intmap

intmap is a python package for mapping the locations of genomic integration from second- and third-generation sequencing data generated via ligation-mediated or linear amplification-mediated PCR approaches.

## Installation
The easiest way to install intmap is by cloning this repository, navigating to the generated 'intmap' directory, and running <code>bash install.sh</code>. This bash script will:

1. Create an intmap conda environment,
2. install all (version-controlled) intmap dependencies, and
3. install the intmap package.

<code>install.sh</code> requires an initialized conda instance. If conda is not available on your computer, you can easily install conda by following the recommendations [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). <code>install.sh</code> can additionally install new versions of intmap seamlessly from either the base environment or an activated intmap environment.

## Package Structure
There are 4 distinct intmap modules: intmap_demux, intmap, intmap_se, and intmap_multi. I will explain these in more detail below.

### intmap_demux
intmap_demux is written to demultiplex multiplexed FASTQ files based on unique integrant- and/or linker-end sequences. Behind the scenes, intmap_demux uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/#) to look for target sequences and allocate them to the appropriate file. 

NOTE: integrant-end sequences are universally referred to as 'ltr' in intmap syntax. However, these sequences do not have to be retrovirus/retrotransposon-specific. Any integrant-specific sequence is valid.

Given multiplexed paired-end FASTQ files file_R2.fq.gz and file_R2.fq.gz, intmap_demux can be run like so:

```
intmap_demux \
    -r1 file_R1.fq.gz \
    -r2 file_R2.fq.gz \
    -bc_r1 ltr_bc.fa \
    -bc_r2 linker_bc.fa
```

This command would look for sequences on both R1 and R2. Options are avaible to restrict the search to one of the two reads. The FASTA files ltr_bc.fa and linker_bc.fa contain the matched integrant- and linker-end sequences that the program looks for. Note that all FASTA headers must be unique, as they are used to name the respective demultiplexed files. By default, intmap_demux allows for a substitution rate of 0.1. This value can be changed with the <code>-e</code> argument. When cutadapt is finished running, intmap_demux will place the proper R1 and R2 files into a subdirectory called 'demux' within the current working directory. These files can then be used directly with intmap. Note that multiplexed FASTQ files do not have to be demultiplexed with intmap_demux to be compatible with intmap. Any demultiplexing procedure will suffice, provided the integrant- and linker-end sequences immediately adjacent to the integrant-host genome junction are retained after demultiplexing.

### intmap and intmap_se

Once the multiplexed FASTQ files have been demultiplexed, integration sites can be mapped with intmap or intmap_se. The difference between these two modules is that intmap expects paired-end data and intmap_se expects single-end data (both short- and long-read data are considered single-end in intmap). Beyond that, their functionality is identical. Below is a minimal command to run intmap. In addition to this core set of arguments, there are over 60 additional user-adjustable parameters that can help fine-tune the analysis to match most experimental circumstances and analytical considerations.

```
intmap \
    -r1 /path/to/demux/sample1_R1.fq.gz \
    -r2 /path/to/demux/sample1_R2.fq.gz \
    -ltr5 TGTGGAAA \
    -ltr3 ATCTCTAGCA \
    -linker5 TCAGGCATC \
    -linker3 GACGTCTCAG \
    -nm sample1 \
    -bt2_idx_dir /path/to/host/genome/index/directory \
    -bt2_idx_name host_index_name \
    -v_idx_dir /path/to/virus/genome/index/directory \
    -v_idx_name virus_index_name
```

Broadly speaking, intmap automates 4 distinct processes. These are: removal of non-host sequences from the ends of the provided reads (referred to in the program as "cropping"), alignment of the cropped reads to the host genome, post-alignment read processing, and coordinate retrieval. The cropping step looks for the concatenated ltr5 and ltr3 sequences in R1 and the concatenated linker5 and linker3 sequences in R2, when present. The identified sequences are removed from the respective reads and the reads are output to new *_cropped.fq.gz files. These new files are then given to [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment. Alignment is a two-step process. The reads are first aligned to the provided virus genome/vector sequence to filter out artifactual reads. The remaining unaligned reads are then aligned to the provided host genome. Aligned reads are then QC'd for sequence quality, read orientation, the number of mismatches/indels, etc. Reads passing QC are then deduplicated. Finally, the coordinates of the final mapped fragments and integration sites are returned.

Integrant- and linker-end search sequences are given to intmap in two chunks. The reason for this is more granular error control in certain parts of the search sequences. Following is an example of what I mean by "two chunks". Suppose that the integrant terminus is 'TGTGGAAAATCTCTAGCA'. I could split this sequence into 2 chunks in a number of ways, but intmap requires that each chunk be at least 5 nt long. Arbitrarily, I could decide that the first chunk will be 'TGTGGAAA' and the second will be 'ATCTCTAGCA'. Since 'TGTGGAAA' comes from the 5' end of the target sequence, it is given as ltr5. 'ATCTCTAGCA' is given as ltr3, as it comes from the 3' end. intmap concatenates these sequences when it looks for the target sequence, but to allow for different error rates in the two chunks, intmap builds regex patterns of the two chunks separately before concatenation. The target linker sequence is defined in a similar way.

### intmap_multi
In most situations, integration site mapping experiments are designed to include more than one sample. To streamline the mapping of integration sites from multiple samples, users can use intmap_multi to serially analyze a given set of samples in a single command. The parameters for the respective samples can be input into intmap_multi interactively, but it is more likely that a setup file will be given to the program. The setup file should be structured in the following way:

```
3
r1, r2, ltr5, ltr3, linker5, linker3, nm, nthr, bt2_idx_dir, bt2_idx_name, v_idx_dir, v_idx_name
demux/sample1_R1.fq.gz, demux/sample2_R1.fq.gz, demux/sample3_R1.fq.gz
demux/sample1_R2.fq.gz, demux/sample2_R2.fq.gz, demux/sample3_R2.fq.gz
TGTGGAAA \
ATCTCTAGCA \
TCAGGCATC \
GACGTCTCAG \
sample1, sample2, sample3
4
/path/to/host/genome/index/directory
host_index_name
/path/to/virus/genome/index/directory
virus_index_name
```

This file is intended to be highly flexible. There are 4 requirements:

1. The first line must be an integer indicating the number of samples being analyzed.
2. The second line must be a comma-delimited set of arguments to be given to intmap for each sample being analyzed.
3. Every line after the second must contain either a single value or a comma-delimited set of values with a total number of values equal to the number given in line 1.
4. No line can be blank. Boolean values should be given as True or False, the program will correctly parse them.

This structure allows the user flexibility in deciding which arguments to pass to intmap (beyond the default values), and to change values between analyses. Once the setup file is constructed, intmap_multi can be run as follows:

```
intmap_multi \
    -s setup_file.txt \
    -n target_params
```

This command will create the file target_params.json and use the entries therein to run intmap for each enumerated sample. If the output JSON file already exists, intmap_multi will not overwrite it (the program will stop). The parameter JSON can also be used with intmap_multi directly:

```
intmap_multi \
    -a target_params.json
```

For mapping single-end data with intmap_multi, the <code>--single_end</code> flag must be set when intmap_multi is called.

### Wrapper script
A bash script -- wrapim.sh -- is provided to wrap demultiplexing with intmap_demux and serialized analyses with intmap_multi into a single command. That script is copied to '~/CONDA_PREFIX/intmap/bin/' when intmap is installed via the install script. The script takes as arguments the multiplexed sequencing data, the barcode FASTA file(s), the intmap_multi setup file, the intmap_multi JSON name, and the name of the output file to store the analysis information. Other arguments are also accepted. Run <code>wrapim -h</code> for more information. Below is an example of using wrapim for paired-end data:

```
wrapim \
    -r1 multiplexed_R1.fq.gz \
    -r2 multiplexed_R2.fq.gz \
    -bc_r1 ltr_bc.fa \
    -bc_r2 linker_bc.fa \
    -setup example_setup.txt \
    -json example_params \
    -out example_out.txt
```

## Licensing
This software is licensed under the GPL-3 License (see LICENSE). It additionally relies on several third-party tools, which are governed by their own licenses. See <code>THIRD_PARTY_LICENSES.md</code> for details.