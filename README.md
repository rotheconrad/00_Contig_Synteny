# Workflow to explore shared genes and gene synteny for a group of metagenome assembled contigs in fasta format.

This workflow produces a gene synteny or gene neighborhood style plot placing arrows in the position and orientation of predicted gene regions on each contig in the set. Contigs are ordered along the y-axis by the user defined fasta file order, or, optionally, they can be ordered by contig length with the -y flag. Contigs are aligned along the x-axis by the first gene cluster that is shared between all contigs in the set, or, optionaly, they can be aligned by a user defined gene with the -x flag. The plot is output in png format and the legend is written to a separate png file.

## Step 00: Required tools :: Python 3.6+, Blast+, Prodigal and CD-HIT.

### Python 3.6+ for running the Python scripts in this repo.

Information for installing and running Python can be found [here](https://www.python.org/). I recommend installing [mini conda](https://docs.conda.io/en/latest/miniconda.html) first and then creating an environment for Python 3.6+ and other tools for the project at hand.

*All Python scripts in this repo were written for Python 3.6+. If you get a syntax error the first time you run a script, please first check your Python version.*

### Blast+

Installation Details for Blast+ can be found [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

Or Blast+ can be easily installed using a [conda environment](https://docs.conda.io/en/latest/miniconda.html):

```bash
conda create -n blastplus
conda activate blastplus
conda install bioconda::blast=2.7.1 conda-forge::gnutls conda-forge::nettle
```

### Prodigal for protein coding gene prediction.
 
Information and installation instructions for Prodigal can be found [here](https://github.com/hyattpd/Prodigal). The publication is [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/).

Prodigal can also be installed with a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```bash
conda create -n prodigal
conda activate prodigal
conda install -c bioconda prodigal
```

### CD-HIT to cluster predicted gene sequence.

Information and installation for CD-HIT can be found [here](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide). The publication can be found [here](https://academic.oup.com/bioinformatics/article/22/13/1658/194225).

CD-HIT can also be installed with a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```bash
conda create -n cdhit
conda activate cdhit
conda install -c bioconda cd-hit
```

## Step 01: Collect contigs to test in fasta format

The recommended strategy here is to use blast or some other sequence search tool to find contigs assembled from one or metagenomic sample that have a region which matches a gene of interest. So first you need to have a gene of interest. Then you metagenome assembles. Then you make a Blast+ (or other search tool) database from your contigs. Then you search your gene against the contigs. Then you filter results to find good matches you are satisfied with. Then you collect the contig sequences and arrange them in fasta format.

*For this pipeline to work, all contigs provided must share at least 1 gene*

## Step 02: Predict gene regions on selected contigs with Prodigal.

## Step 03: Cluster predicted genes with CD-HIT (or CD-HIT-EST).

Depends on level to cluster to. Nucleotide or Amino Acid. 90% or 40% etc..

## Step 04: Build the plot.