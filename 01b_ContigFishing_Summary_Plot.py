#!/usr/bin/env python

'''Plots swarm plot from tabular blast file.

Assuming a select set of genes or contigs are used as blast query
sequences, and assembled metagenome or draft genome contigs are used as
the database or subject sequences, this script plots the percent IDs of
the query matches along the x-axis for each query along the y-axis to
show the distribution of similar sequences that are present in an
assembly.

It is assumed the blast table has been filtered for quality matches.

Optionally, an annotation metadata file can be given to match the query
sequence names to a short or long gene name and to specify a custom 
color. The metadata file format is a comma separated file with a header
and three columns. A single metadata value may be given by leaving the 
other column blank. Example:

query name, annotation, color
name_one, geneX, #006837
name_two, , #0868ac
name_three, geneY, 

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 3rd, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_meta(meta):

    # Read meta file and store key: value pairs of QueryID: [Annotation, color]
    meta_dict = {} # initialize dict
    # open file
    with open(meta, 'r') as f:
        # skip/set header
        header = f.readline()
        # iterate through lines
        for l in f:
            # split each line by tab into list X
            X = l.rstrip().split(', ')
            # select subjectID from the first column
            queryID = X[0]
            # select annotation
            if not X[1].isspace(): Annotation = X[1]
            else: Annotation = queryID
            # select color
            if len(X) == 3: color = X[2]
            else: color = ''
            # populate dictionary
            meta_dict[queryID] = [Annotation, color]

    return meta_dict


def parse_data(blst, meta):
    ''' reads blast table and matches query sequence names to the 
        anno file (annotation file). The anno file should be a comma
        separated file with two or three columns. The third column is for
        user defined colors but can be left blank for default colors.'''

    # initialize dict of Column Names: List of values 
    # this will be converted to a pandas dataframe to build the plot
    blast_dict = {'Gene': [], 'Percent ID': []}

    # if meta data is specified read in the info
    if meta:
        meta_dict = parse_meta(meta)
        blast_dict['Colors'] = []

    # open file
    with open(blst, 'r') as f:
        # iterate through lines
        for l in f:
            # split each line by tab into list X
            X = l.rstrip().split('\t')
            # select subjectID from 2nd column
            queryID = X[0]
            # select pID from 3rd column
            pID = float(X[2])
            # if meta data is given match queryID to meta data
            if meta:
                Annotation = meta_dict[queryID][0]
                color = meta_dict[queryID][1]
                # populate dictionary
                blast_dict['Gene'].append(Annotation)
                blast_dict['Percent ID'].append(pID)
                blast_dict['Colors'].append(color)
            # otherwise use what is given
            blast_dict['Gene'].append(queryID)
            blast_dict['Percent ID'].append(pID)

    # convert blast_dict to pandas dataframe
    df = pd.DataFrame(blast_dict)

    # return dataframe
    return df


def contig_fishing_plot(df, outfile):
    ''' Assumes a blast table has been loaded as a pandas dataframe.
        Builds a swarm plot of percent ID for each unique subject ID '''

    # Check for colors

    # initialize the plot
    fig, ax = plt.subplots(figsize=(10,5))

    # build the plot
    ax = sns.swarmplot(x='Percent ID', y='Gene', data=df, ax=ax)

    # set grid style
    ax.minorticks_on()
    ax.tick_params(axis='y', which='minor', left=False)
    ax.xaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=1
        )
    ax.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1.5
        )
    ax.set_axisbelow(True)

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(outfile)
    plt.close() 

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--tab_blast_file',
        help='Please specify the tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--metadata_file',
        help='(Optional) metadata file containing annotation and color.',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')

    blst = args['tab_blast_file']
    meta = args['metadata_file']
    outfile = blst.split('.')[0] + '.png'

    # Read input files and parse the data
    df = parse_data(blst, meta)
    # Build the plot
    contig_fishing_plot(df, outfile)

if __name__ == "__main__":
    main()
