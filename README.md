# Workflow to explore shared genes and gene synteny for a group of metagenome assembled contigs (or draft genome assembly contigs) in fasta format.

** WORK IN PROGRESS - THIS PROTOCOL IS NOT COMPLETE **

This workflow produces a gene synteny or gene neighborhood style plot placing arrows in the position and orientation of predicted gene regions on each contig in the set. Contigs are ordered along the y-axis by the user defined fasta file order, or, optionally, they can be ordered by contig length with the -y flag. Contigs are aligned along the x-axis by the first gene cluster that is shared between all contigs in the set, or, optionaly, they can be aligned by a user defined gene cluster with the -x flag. The plot is output in png format and the legend is written to a separate png file.

The idea is to identify a set of contigs containing a gene of interest using Blast+ (or other sequence alignment tool), predict CDS regions on the contigs (with Prodigal), cluster the predicted genes (with CD-HIT), and then plot the results for visual interpretation (with a Python script from this GitHub Repo). In the case that the sequence search identifies a certain gene on all contigs, but the gene prediction and clustering tools do not find the same gene cluster to be on all contigs, the contig(s) missing the alleged shared gene cluster will be centered on the x-axis for visual inspection compared to the other contigs. The CD-HIT clstr file can be modified to merge or adjust gene clusters manually if needed.

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

The recommended strategy here is to use Blast+ or some other sequence search tool to identify contigs assembled from one or more metagenomic samples or draft genomes that all share sequence similarity to a gene(s) of interest. So first you need to have a gene(s) of interest (this part is not covered in the protocol). You also need metagenome or draft genome assemblies (also not covered in this protocol). Once you have these, make a Blast+ (or other search tool) database with the assemblies and search your gene(s) of interest against the contigs. I typically work with the tabular blast output format and customize the options to output the query and sequence lengths in columns 13 and 14. Using the tabular blast output with added sequence lengths I can filter the results for matches with a good (alignment length / query sequence length) ratio and also select contigs that are longer. Longer contigs will have more information about gene synteny and neighborhoods.

#### Example Command: Tabular Blast+ with added sequence lengths

```bash
blastn -task 'megablast' -evalue 0.01 -max_target_seqs 20 -num_threads 2 -db {pathto_db} -query {input_fasta} -out {outfile_name} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
```

If you've used the above sequence search method, you can filter your results with 01a_BlastTab_BestHit_Filter.py Python script included in this GitHub Repo. You can decide the minimum percent sequence similarity (Default = 70%), the (alignment length / query sequence length) ratio * 100 (Percent match length: Default = 50), and minimum contig size (Default = 1000 base pairs).

#### Example Command: Filter tabular Blast+ output with added sequence length

```bash
# To Display the program description and parameter options
python 01a_BlastTab_BestHit_Filter.py -h

# Example execution:
python 01a_BlastTab_BestHit_Filter.py -i tabular_blast_output.blast

# Loop through a list of tabular blast files:
for file in *blast; do echo $file; python 01a_BlastTab_BestHit_Filter.py -i $file; done
```

If you want to visualize the distribution of sequence similarity for genes of interests in the metagenome or draft genome assembly, the 01b_ContigFishing_Summary_Plot.py script included in this GitHub Repo provides a good image.

Assuming a select set of genes are used as blast query sequences, and assembled metagenome or draft genome contigs are used as the database or subject sequences, this script plots the percent IDs of the query matches along the x-axis for each query along the y-axis to show the distribution of similar sequences that are present in an assembly.

It is assumed the blast table has been filtered for quality matches.

Optionally, an annotation metadata file can be given to match the query sequence names to a short or long gene name and to specify a custom  color. The metadata file format is a comma separated file with a header and three columns. You can specify colors and some annotations by leaving the blanks in the second column. You can specify annotations and use default colors by leaving all color fields blank.

Example:

    query name, annotation, color
    name_one, geneX, #006837
    name_two, , #0868ac
    name_three, geneY, #e6550d

#### Example Command: Plot distribution of sequence similarity to each gene of interest in the metagenome or draft genome assembly.

```bash
# To Display the program description and parameter options
python 01b_ContigFishing_Summary_Plot.py -h

# Example execution:
python 01b_ContigFishing_Summary_Plot.py -b tabular_blast_output.fltrdBstHts.blst
```

Once you have a currated list of contigs that you want to explore, you need to retrieve the sequences and place them in a fasta formatted file. This can be accomplished with ####### Python script included in this GitHub Repo.

#### Example Command: Retrieving selected contigs in fasta format

```bash
# To Display the program description and parameter options
python 01c -h

# Example execution:
python 01c
```


## Step 02: Predict gene regions on selected contigs with Prodigal.

Now that you have a fasta file containing all the contigs that potentially contian a copy of your gene(s) of interest, predict all the CDS regions on the contig set. I like to have Prodigal output a GFF file, a fasta file of nucleotide sequences (.fna), and a fasta file of translated amino acid sequences (.faa) so I have them if I need want them at a later point.

#### Example Command: Predict CDS regions with Prodigal

```bash
prodigal -i contig_set.fasta -o my_genes -a my_proteins.faa -d my_genes.fna
```

## Step 03: Cluster predicted genes with CD-HIT (or CD-HIT-EST).

CD-HIT has many options. You can cluster your predicted gene sequences using the nucleotide sequence (with CD-HIT-EST) for higher resolution with less divergent gene sequences, or you can use CD-HIT with amino acid sequences for more divergent gene sequences. For now we will use amino acid sequence and cluster using 40% sequence identity with an alignment over at least 70% of the shorter sequence. It is good to think about and test the best parameters for your project and questions.

```bash
cd-hit -i my_proteins.faa -o my_proteins.cdhit.faa -c 0.4 -n 2 -G 0 -g 1 -aS 0.7 -d 0
```

This will output two files. The my_proteins.cdhit.faa file contains representative (longest) protein sequence for each gene cluster. These sequences can be used for Step 05 to add functional annotation information to the gene cluster. The my_proteins.cdhit.faa.clstr contains the actual gene cluster information. We will use the clstr file to build the plot in Step 04.

## Step 04: Build the plot.

The ###### script included in this GitHub Repo will take the contig_set.fasta, the my_proteins.faa and the my_proteins.cdhit.faa.clstr file and the build the plot.



```bash
python 
```

## Step 05: Customize and Add Annotations







