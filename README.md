# Analysis pipeline for BLISS data

<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Setup](#setup)
* [Usage](#usage)
* [Misc](#misc)

<!-- vim-markdown-toc -->

Description
===========

A fairly generic pipeline for the detection of double-strand DNA breaks from
[BLISS](https://www.nature.com/articles/s41596-020-0397-2) data.

Input is single or paired-end fastq files from BLISS libraries and reference
genome.

The main output is aligned and de-duplicated bam files, peak files, bigwig
files and a few diagnostic plots.

Setup
=====

If not already done, install `conda`, `mamba`, and configure for
[bioconda](https://bioconda.github.io/). Check the proper documentation, but
something like this should do on Linux:

```
# Install conda
curl -O -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh # Follow interactive instructions, default settings probably ok.

# Install mamba
conda install mamba -n base -c conda-forge

# Configure bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Create an environment dedicated to this project. Change the project name
`bliss` to something more appropriate if you wish:

```
conda create --yes -n bliss
conda activate bliss
mamba install -n bliss --yes --freeze-installed --file requirements.txt
```

Usage
=====

Edit the [sample_sheet.tsv](sample_sheet.tsv) with your samples, genomes and
fastq files. Edit [genomes.tsv](genomes.tsv) to indicate where to find the
reference genome.

The following command runs the pipleine on the test data in dry-run mode (`-n` option):

```
snakemake -p -n -j 10 -d tmp --nt \
   -C ss=$PWD/test/sample_sheet.tsv \
      genomes=$PWD/test/genomes.tsv \
      min_depth=0
```

The real analysis may be somethng like this. Remove `-n` for actual execution and edit `-d` to point to the output directory:

```
snakemake --rerun-trigger mtime -p -n -j 10 -d output/20211214_ross_bliss \
   -C ss=$PWD/sample_sheet.tsv \
      genomes=$PWD/genomes.tsv \
      min_depth=10
```

Misc
====

To compile this markdown to pdf:

```
pandoc -V colorlinks=true -V geometry:margin=1in README.md -o README.pdf
```
