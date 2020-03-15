#!/usr/bin/env python

'''Build a contig synteny plot.

This script needs 3 files:
    - Fasta file of the contig sequences
    - Prodigal file of nucleotide gene sequences in fasta format
    - CD-HIT-EST clstr file result from Prodigal gene predictions

This script writes a gene synteny plot as a .png to output file.

Optional annotation file format is comma separated in the order:
Cluster_Number, Gene_Name, Hex_Color

The user can provide a gene annotation with the Gene_Name position or 
specifiy the color for the gene cluster with the Hex_Color position.
To specify the color without annotation information just repeat the
Cluster_Number for the Gene_Name position. Ex:

Cluster_0, NifH, #31a354
Cluster_1, NifG, #756bb1
Cluster_2, IS1 Transposase, #de2d26
Cluster_3, Cluster_3, #3182bd
Cluster_4, ABC Transporter, #e6550d
Cluster_5, Zinc finger BED domain containing protein, ##636363

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 5th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
from collections import defaultdict
from collections import OrderedDict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

def read_fasta(fp):
    '''Parses fasta file format and returns name, seq strings'''

    # initialize name, seq objects
    name, seq = None, []
    # iterate through lines in file
    for line in fp:
        # strip new line characters
        line = line.rstrip()
        # check if line is read name
        if line.startswith(">"):
            # return prior name, seq object at each new name
            if name: yield (name, ''.join(seq))
            # define new name, empty seq object for next sequence
            name, seq = line, []
        # else it is sequence so add to current seq object
        else:
            seq.append(line)
    # this returns the last name, seq pair of the file.
    if name: yield (name, ''.join(seq))


def parse_contig_file(infile):
    '''Parses contig fasta file. Returns dict sorted by contig length:
    {contig_name: length}'''
    
    # initialize ordered dictionary to preserve user provided order
    contigs = OrderedDict()

    # open file
    with open(infile, 'r') as file:
        # iterate through fasta file
        for name, seq in read_fasta(file):
            # populate dictionary with key: value = contig_name: length
            contigs[name[1:]] = len(seq)

    # sort contigs by length - longest contig to smallest
    cntg_order = [ 
        k for k, v in sorted(
                        contigs.items(),
                        key=lambda item: item[1],
                        reverse=True
                            )
                                ]
    # use to center contigs on longest and to set figure/plot width.

    return contigs, cntg_order


def parse_gene_file(infile):
    '''Parses gene fasta file. Returns dict
       {gene_name: [start, stop, strand, partial]}'''
    
    # initialize dictionary
    genes = {}

    # open file
    with open(infile, 'r') as file:
        # iterate through fasta file
        for name, seq in read_fasta(file):
            # select values from name string
            X = name.split(' ')
            gene_name = X[0][1:]
            start = int(X[2])
            stop = int(X[4])
            strand = X[6]
            partial = X[8].split('=')[2].split(';')[0]
            genes[gene_name] = [start, stop, strand, partial]

    return genes


def parse_cluster_file(infile):
    '''Parses cd-hit clstr file. Returns dict
       {cluster_number: [gene names in cluster]}'''
    
    # initialize dictionary
    clusters = {}

    # open file
    with open(infile, 'r') as file:
        # iterate through cdhit clstr file
        for line in file:
            # look for start of cluster, set cluster number
            if line.startswith('>'):
                cluster_number = f"Cluster_{line.rstrip().split(' ')[1]}"
                clusters[cluster_number] = []
            # otherwise add gene name to current cluster
            else:
                gene_name = line.split('>')[1].split('...')[0]
                clusters[cluster_number].append(gene_name)

    # which gene cluster is represented on the most contigs?
    v = max(clusters, key=lambda key: len(clusters[key]))
    clstr_order = [ 
                k for k, v in sorted(
                                clusters.items(),
                                key=lambda item: len(item[1]),
                                reverse=True
                                    )
                                        ]
    # use to center the synteny plot around the most represented gene cluster.

    return clusters, clstr_order


def order_by_contig_length(axs, contigs, cntg_order, xalign, lw, fs):
    '''Plot contigs ordered by length'''

    # track start and end position to get leftmost and rightmost points of plot
    leftmost = []
    rightmost = []
    # initialize dictionary to store contig stop start sites on plot
    contig_pos = {}
    # set number of congtigs and length of longest contigs
    total_contigs = len(cntg_order)
    longest_contig = cntg_order[0]
    longest_contig_length = contigs[longest_contig]
    ypos = 0

    # iterate through remaining contigs and plot lines for each.
    for i, cntg in enumerate(cntg_order):
        # set current contig length
        current_contig_length = contigs[cntg]
        # if contig has gene from most representative gene cluster
        # center contig on the gene
        if cntg in xalign:
            xstart = -xalign[cntg]
            xend = current_contig_length - xalign[cntg]
        # if it doesn't, center contig at it's center.
        else:
            # round current contig length to even value.
            if current_contig_length % 2 != 0: current_contig_length += 1
            # divide it in half to center contig on zero.
            half_contig = int(current_contig_length / 2)
            # define start and stop of centered contig
            xstart = -half_contig
            xend = half_contig
        # append xstart to leftmost to find minium x value later
        leftmost.append(xstart)
        # append xend to rightmost to find maximum x value later
        rightmost.append(xend)
        # define xvalues to plot contig
        xvalues = range(xstart, xend)
        # record contig position 
        contig_pos[cntg] = [xstart, xend, ypos, i]
        # define yvalues to plot contig
        yvalues = [ypos] * current_contig_length
        # plot contig
        axs[i].plot(xvalues, yvalues, linestyle='-', lw=lw, color='#252525')

    # Write contig names at left edge of plot
    xlabpos = min(leftmost) - (longest_contig_length*.01)
    # Get max length of combined contigs for x-axis
    xmaxlen = max(rightmost) - min(leftmost)

    # Write contig names at left edge of plot
    for i, cntg in enumerate(cntg_order):
        axs[i].text(
            xlabpos, ypos, cntg,
            fontsize=fs, color='#252525',
            horizontalalignment='right',
            verticalalignment='center'
            )

    return axs, contig_pos, xmaxlen


def order_by_user(axs, contigs, cntg_order, xalign, lw, fs):
    '''Plot contigs in order user provided'''

    # track start and end position to get leftmost and rightmost points of plot
    leftmost = []
    rightmost = []
    # initialize dictionary to store contig stop start sites on plot
    contig_pos = {}
    # set number of congest and length of longest contigs
    total_contigs = len(cntg_order)
    longest_contig = cntg_order[0]
    longest_contig_length = contigs[longest_contig]
    #yscaler = longest_contig_length*.02
    ypos = 0

    # iterate through contigs and plot lines for each.
    for i, (cntg, length) in enumerate(contigs.items()):
        current_contig_length = length
        xstart = -xalign[cntg]
        xend = current_contig_length - xalign[cntg]
        # append xstart to leftmost to find minium x value later
        leftmost.append(xstart)
        # append xend to rightmost to find maximum x value later
        rightmost.append(xend)
        # define xvalues to plot contig
        xvalues = range(xstart, xend)
        # record contig position
        contig_pos[cntg] = [xstart, xend, ypos, i]
        # define yvalues to plot contig
        yvalues = [ypos] * current_contig_length
        # plot contig
        axs[i].plot(xvalues, yvalues, linestyle='-', lw=lw, color='#252525')

    # Write contig names at left edge of plot
    xlabpos = min(leftmost) - (longest_contig_length*.01)
    # Get max length of combined contigs for x-axis
    xmaxlen = max(rightmost) - min(leftmost)

    # Write contig names at left edge of plot
    for i, (cntg, length) in enumerate(contigs.items()):
        axs[i].text(
            xlabpos, ypos, cntg,
            fontsize=fs, color='#252525',
            horizontalalignment='right',
            verticalalignment='center'
            )

    return axs, contig_pos, xmaxlen


def parse_annotation_file(infile):
    '''Parses the user provided annotation file. This file should be
       three comma separated files with columns of:
       Cluster Number, Gene Annotation, Color'''

    # initialize dictionary to store annotations
    annotations = {}

    # open file
    with open(infile, 'r') as file:
        # iterate through file
        for line in file:
            # split each line by comma
            X = line.rstrip().split(', ')
            cluster = X[0]
            anno = X[1]
            color = X[2]
            annotations[cluster] = [anno, color]

    return annotations


def find_first_core_clust(clusters, cntg_order):
    '''Finds the 1st cluster that contains a gene from all contigs'''
    # initialize variable to find the first core cluster.
    first_core_clust = None
    # initialize dict to keep track of how many genomes are represented
    # by each cluster.
    cntgs_per_cluster = {}
    # read through the cluster dictionary and calculate the number of
    # genome represented in each gene cluster
    # clusters dict is of {cluster_number: [list of genes]}
    for clust, gene_list in clusters.items():
        # for each cluster store a list of genomes in that cluster
        cntg_list = {}
        # for each gene in the cluster
        for g in gene_list:
            # get the contig name from the gene name
            c = '_'.join(g.split('_')[:-1])
            # add contig name to cntg_list dict as key with blank value
            # this builds a dereplicated list
            cntg_list[c] = ''
        # count the number of contigs with a gene in the cluster.
        cntg_count = len(cntg_list)
        # Add cluster as key and count of represented contigs.
        cntgs_per_cluster[clust] = cntg_count
        # check if cntg_list is the same length as cntg_order
        # cntg_order contains a full list of the contigs
        # if true the gene cluster contains at least one gene from
        # each contig in the set. This is the first_core_clust.
        # Set the cluster_number and use this to align contigs on
        # the x-axis (default behavior if user does not specify)
        if len(cntg_order) == len(cntg_list):
            first_core_clust = clust

    # if a cluster was not found that contains at least 1 gene for all
    # contigs, then set first_core_clust equal to a gene cluster 
    # containing a gene from the maximum number of contigs in the set.
    if not first_core_clust:
        print(
            '\nGene Cluster Error:\n'
            'No Gene Cluster found that is represented on all contigs.\n'
            'The gene cluster containing genes from the most contigs '
            'will be used for\nx-axis alignment. Contigs without a gene '
            'in this gene cluster will be centered\non the x-axis out of '
            'alignment with the contigs containing a gene from this\n'
            'most representative gene cluster.\n\n'
            )
        # Get a cluster with maximum number of contigs represented.
        max_rep_clust = max(
            cntgs_per_cluster, key=lambda k: cntgs_per_cluster[k]
            )
        # asign this cluster as the first_core_clust to use for x-axis alignment
        first_core_clust = max_rep_clust


        # print out contigs in the cluster
        # print out contigs not in the cluster
        # 1st select contig names from gene names in max_rep_clust
        maxrep_cntg_list = [
            '_'.join(g.split('_')[:-1]) for g in clusters[first_core_clust]
            ]
        # initialize lists for contigs in or out of cluster
        cntg_in = []
        cntg_out = []
        # iterate through all contigs with cntg order
        for c in cntg_order:
            # if c in cluster
            if c in maxrep_cntg_list: cntg_in.append(c)
            # if c not in cluster
            else: cntg_out.append(c)

        print('Contigs with gene from the most represented gene cluster:')
        for c in cntg_in: print(c)
        print('\n\nContigs without this gene:')
        for c in cntg_out: print(c)
        print('\n\n')

    return first_core_clust

def contig_synteny_plot(
        contigs, cntg_order, genes, clusters, clstr_order, outfile,
        annotations, sorty, alignx, width, height
                    ):
    '''Builds a contig synteny plot and writes .png to file'''
    
    # set default colors to use for genes
    colors = [
            '#08306b', '#67000d', '#00441b', '#3f007d', '#7f2704', '#1a1a1a',
            '#08519c', '#a50f15', '#006d2c', '#54278f', '#a63603', '#252525',
            '#2171b5', '#cb181d', '#238b45', '#6a51a3', '#d94801', '#525252',
            '#4292c6', '#ef3b2c', '#41ab5d', '#807dba', '#f16913', '#737373',
            '#6baed6', '#fb6a4a', '#74c476', '#9e9ac8', '#fd8d3c', '#969696',
            '#9ecae1', '#fc9272', '#a1d99b', '#bcbddc', '#fdae6b', '#bdbdbd',
            '#c6dbef', '#fcbba1', '#c7e9c0', '#dadaeb', '#fdd0a2', '#d9d9d9',
            '#deebf7', '#fee0d2', '#e5f5e0', '#efedf5', '#fee6ce', '#f0f0f0',
            '#f7fbff', '#fff5f0', '#f7fcf5', '#fcfbfd', '#fff5eb', '#ffffff'
            ]

    # set font size
    fs = 18
    # set line width
    lw = 2

    # initialize the plot objects
    total_contigs = len(cntg_order)
    total_clusters = len(clusters)
    longest_contig_length = contigs[cntg_order[0]]

    # set default figure width and height if not specified
    if not height: height = total_contigs
    if not width: width = 30

    # build separate figure for the legend
    fig_legend, ax_legend = plt.subplots(figsize=(width,total_clusters*1.5), dpi=300)
    ax_legend.set_axis_off()

    # Build a figure with height and ncols equal to number of contigs.
    fig, axis_array = plt.subplots(
                                nrows=total_contigs, ncols=1,
                                squeeze=False, # make 2D for nrows and ncols
                                figsize=(width,height), dpi=300,
                                sharex=True, #sharey=True
                                )
    axs = axis_array.flat # flatten the axis array to loop over it.

    # Set Plot and axis titles
    #plt.suptitle('Contig Synteny Plot', fontsize=28, y=1.02)

    # Removes Spines ticks and labels. Set axis equal length to preserve scale.
    for ax in axs:
        ax.axis('equal')
        ax.set_axis_off()

    #####################################################################
    ################## CONSIDER MOVING TO FUNCTION ######################
    #####################################################################
    # Fuction would look like this:
    #xalign = parse_gene_positions(genes, clusters, alignx, clstr_order)

    # initialize dictionary to store x alignment info for each contig
    xalign = {}
    # Select genes from cluster to align and parse position information.
    if alignx:
        selected = clusters[alignx]
        for gene in selected:
            # select contig name from gene name
            contig = '_'.join(gene.split('_')[:-1])
            # add contig name to xalign with gene position to align to.
            xalign[contig] = genes[gene][0]
    # Set default alignx cluster to largest cluster if user did not specify.
    else:
        first_core_clust = find_first_core_clust(clusters, cntg_order)
        selected = clusters[first_core_clust]
        #selected = clusters[clstr_order[0]]
        for gene in selected:
            # select contig name from gene name
            contig = '_'.join(gene.split('_')[:-1])
            # add contig name to xalign with gene position to align to.
            xalign[contig] = genes[gene][0]
    #####################################################################
    #####################################################################


    # Lay down contig order and alignment based on user preference
    if sorty:
        axs, contig_pos, xmaxlen = order_by_contig_length(
                                        axs, contigs, cntg_order, xalign, lw, fs
                                        )
    else:
        axs, contig_pos, xmaxlen = order_by_user(
                                        axs, contigs, cntg_order, xalign, lw, fs
                                        )

    #####################################################################
    #####################################################################

    # Initialize list for legend elements
    legend_elements = []

    # Add legend key elements
    complete = mpatches.Patch(
                    label='Complete Gene',
                    facecolor='#FFFFFF',
                    edgecolor='#000000',
                    linestyle='-',
                    lw=lw
                    )
    legend_elements.append(complete)

    partial_left = mpatches.Patch(
                    label='Left Side Partial',
                    facecolor='#FFFFFF',
                    edgecolor='#000000',
                    linestyle='--',
                    hatch='\\',
                    lw=lw
                    )
    legend_elements.append(partial_left)

    partial_right = mpatches.Patch(
                    label='Right Side Partial',
                    facecolor='#FFFFFF',
                    edgecolor='#000000',
                    linestyle='--',
                    hatch='/',
                    lw=lw
                    )
    legend_elements.append(partial_right)

    partial_both = mpatches.Patch(
                    label='Both Side Partial',
                    facecolor='#FFFFFF',
                    edgecolor='#000000',
                    linestyle='--',
                    hatch='x',
                    lw=lw
                    )
    legend_elements.append(partial_both)

    #####################################################################
    ################## CONSIDER MOVING TO FUNCTION ######################
    #####################################################################

    # scale arrow height based on x-axis max length
    arrowhieght = xmaxlen * 0.025
    arrowhalf = arrowhieght / 2
    arrownudge = arrowhieght * 0.01
    print('X-axis Max Length:', xmaxlen)
    print('Arrow Height:', arrowhieght)

    # Draw arrows on contigs for each gene cluster
    # Read through clusters dictionary of cluster_number: gene_list
    for i, (clust, gene_list) in enumerate(clusters.items()):
        # set the cluster number
        cluster_number = clust.split('_')[1]

        # choose the color and legend label
        if annotations:
            cluster_label = annotations[clust][0]
            cluster_color = annotations[clust][1]
        else:
            cluster_label = clust
            if i >= len(colors): i = i - len(colors)
            cluster_color = colors[i]

        # Add gene cluster to the legend
        legend_element = mpatches.Patch(
                                label=cluster_label,
                                facecolor=cluster_color,
                                )
        legend_elements.append(legend_element)

        # for each gene in the cluster parse values and plot arrow
        for gene in gene_list:
            # define contig name from gene name
            contig = '_'.join(gene.split('_')[:-1])
            X = contig_pos[contig]
            contig_start = X[0]
            contig_end = X[1]
            contig_y = X[2]
            contig_ax = X[3]
            # retrieve gene info from genes dictionary
            X = genes[gene]
            gene_start = X[0]
            gene_stop = X[1]
            strand = X[2]
            partial = X[3]
            # Choose direction
            arrows = {'-1': 'larrow', '1':'rarrow'}
            arrow = arrows[strand]
            # Choose hatch for partial genes
            # 00 - full gene
            # 10 - partial on left side
            # 01 - partial on right side
            hatches = {'00': None, '10': '\\', '01': '/', '11': 'x'} #'
            hatch = hatches[partial]
            lines = {
                '10': 'dashed', '01': 'dashed', '00': 'solid', '11': 'dashed'
                }
            line = lines[partial]
            # Plot the arrow
            x = contig_start + gene_start
            gene_width = gene_stop - gene_start - arrowhalf # arrow head length
            y = contig_y - arrowhalf # center arrow height
            patch = mpatches.FancyBboxPatch(
                    xy=(x,y), width=gene_width, height=arrowhieght,
                    boxstyle=arrow, facecolor=cluster_color, lw=lw,
                    edgecolor='#000000', hatch=hatch, linestyle=line,
                    transform=axs[contig_ax].transData, mutation_aspect=None,
                    mutation_scale=1, bbox_transmuter=None, zorder=3
                        )
            axs[contig_ax].add_patch(patch)
            # Plot cluster number above gene
            text_x = x + (gene_width / 2)
            text_y = contig_y + arrowhieght + arrownudge
            axs[contig_ax].text(
                            text_x, text_y,
                            cluster_number,
                            verticalalignment='center',
                            horizontalalignment='center',
                            transform=axs[contig_ax].transData,
                            fontsize=fs
                            )
    #####################################################################
    #####################################################################

    # Plot the legend in separate legend figure object
    ax_legend.legend(
            handles=legend_elements,
            ncol=1,
            loc='upper left',
            #bbox_to_anchor=(0.5, 0),
            #fancybox=True,
            frameon=False,
            fontsize=44
            )

    # Save the legend separately
    legendout = outfile.split('.')[0]
    fig_legend.set_tight_layout(True)
    fig_legend.savefig(f'{outfile}_legend.png')

    # adjust layout, save, and close
    fig.subplots_adjust(
            # the left side of the subplots of the figure
                    left  = 0.3,
            # the right side of the subplots of the figure
                    right = 0.9999,
            # the bottom of the subplots of the figure
                    bottom = 0.01,
            # the top of the subplots of the figure
                    top = 0.99,
            # the amount of width reserved for blank space between subplots
                    wspace = 0.2,
            # the amount of height reserved for white space between subplots
                    hspace = 0.0
                    )

    fig.savefig(outfile)
    plt.close() 

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-c', '--contig_fasta_file',
        help='Please specify the fasta file with contig sequence!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--prodigal_fasta_file',
        help='Please specify the fasta file with Prodigal gene predictions!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-r', '--cdhit_clstr_file',
        help='Please specify the CD-HIT clster file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='What would you like to call the output file?',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--annotation_file',
        help='(Optional) Annotation file for the legend.',
        #metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-y', '--contig_order',
        help='(Optional) Sort order of contigs on the y-axis by length. '
             'The default is to keep the same order as the contig fasta file',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '-x', '--contig_alignment',
        help='(Optional) Align contigs on the x-axis by cluster provided'
             ' (ex: Cluster_2). The default is to align the genes of '
             'the largest gene cluster. !! All contigs should have a '
             'gene in the selected cluster !!',
        #metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-fw', '--figure_width',
        help='(Optional) Use to change the width of the figure.',
        #metavar='',
        type=int,
        required=False
        )
    parser.add_argument(
        '-fh', '--figure_height',
        help='(Optional) Use to change the height of the figure.',
        #metavar='',
        type=int,
        required=False
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nRunning Script...\n\n')

    # Read in contig file
    contigs, cntg_order = parse_contig_file(args['contig_fasta_file'])

    # Read in prodigal gene predictions
    genes = parse_gene_file(args['prodigal_fasta_file'])

    # Read in cd-hit clusters
    clusters, clstr_order = parse_cluster_file(args['cdhit_clstr_file'])

    # if annotations file provided, need to parse it.
    if args['annotation_file']:
        annotations = parse_annotation_file(args['annotation_file'])
    else:
        annotations = None

    # Build the plot
    _ = contig_synteny_plot(
                        contigs,
                        cntg_order,
                        genes,
                        clusters,
                        clstr_order,
                        args['output_file'],
                        annotations,
                        args['contig_order'],
                        args['contig_alignment'],
                        args['figure_width'],
                        args['figure_height']
                        )

    print('Complete success space cowboy, cowgirl or cowfolk!\n')
    print('Script seems to have finished successfully!\n\n')

if __name__ == "__main__":
    main()
