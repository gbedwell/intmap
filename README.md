# intmap

intmap is a python package for mapping the location of retroviral genomic integration from paired-end Illumina sequencing data.

## Installation
The easiest way to install intmap is by running <code>bash install.sh</code>. This bash script will:

1. Create an intmap conda environment,
2. install all (version-controlled) intmap dependencies, and
3. install the intmap package.

<code>install.sh</code> requires an initialized conda instance. If conda is not available on your computer, you can easily install conda by following the recommendations [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). If you're unfamiliar with conda, be sure to read the messages returned during creation of the conda environment, they're generally pretty informative. <code>install.sh</code> has been successfully executed using conda 24.5.8. I make no promises about older version of conda.

The intmap package itself can be installed by entering the intmap directory and running <code>pip install -e .</code>. This is useful for installing updated versions of the package, etc.

## Package Structure
There are 4 distinct intmap modules. These are: intmap-demux, intmap, intmap-multimap, and intmap-random. I will explain each of these in more detail below.

### intmap-demux
intmap-demux is written for demultiplexing multiplexed Illumina data. Under the hood, it relies on cutadapt to stratify reads into the appropriate bin. By default, intmap-demux expects to be given:

1. multiplexed R1 and R2 fastq files (<code>-r1</code> and <code>-r2</code>, respectively),
2. a fasta file containing the appropriate R1 barcode sequences (<code>-bc_ltr</code>), and
3. a fasta file containing the appropriate terminal R2 linker sequences (<code>-bc_linker</code>).

You will also probably want to set <code>-nthr</code> > 1 to enable parallelization, but this is not strictly required. Other options are also available to set the allowable error rate (mismatches only), demultiplex based only on the LTR-end, etc.

### intmap

