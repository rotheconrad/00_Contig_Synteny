#!/usr/bin/env python

'''Plots swarm plot from tabular blast file.

Assuming a select set of genes or contigs are used as a blast database,
and assembled metagenome contigs are searched against it, this script
plots the percent IDs of the metagenome contig matches to each selected
gene or contig in the database.

It is assumed the blast table has filter for quality matches.
The swarm plot will show a dot for each match and it's corresponding 
precent ID of the sequence aligment.

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

def parse_data(blst, anno):
    ''' reads blast table and matches subject sequence names to the 
        anno file (annotation file). This file is currently named
        04_ClstrRepSeq_PanCat_ANATAD_Annotation and is the output of
        Step 02 of section 08x in the Sruber 00_Working_Notes.txt file.
        It is the output of the following script:
        08c_Match_ANA-TADs_toAnnotations.py '''

    # Read anno file and store key: value pairs of SubjectID: Annotation
    anno_dict = {} # initialize dict
    # open file
    with open(anno, 'r') as f:
        # skip/set header
        header = f.readline()
        # iterate through lines
        for l in f:
            # split each line by tab into list X
            X = l.rstrip().split('\t')
            # select subjectID from the first column
            subjectID = X[0]
            # select long gene name form annotation from 11th column
            Annotation = X[10]
            # populate dictionary
            anno_dict[subjectID] = Annotation

    # Read blst file, match subjectID to Annotation from anno_dict
    blast_dict = {'Annotation': [], 'Percent ID': []} # initialize dict
    # open file
    with open(blst, 'r') as f:
        # iterate through lines
        for l in f:
            # split each line by tab into list X
            X = l.rstrip().split('\t')
            # select subjectID from 2nd column
            subjectID = X[1]
            # select pID from 3rd column
            pID = float(X[2])
            # match subjectID to annotation
            Annotation = anno_dict[subjectID]
            # populate dictionary
            blast_dict['Annotation'].append(Annotation)
            blast_dict['Percent ID'].append(pID)

    # convert blast_dict to pandas dataframe
    df = pd.DataFrame(blast_dict)

    # return dataframe
    return df


def contig_fishing_plot(df, outfile):
    ''' Assumes a blast table has been loaded as a pandas dataframe.
        Builds a swarm plot of percent ID for each unique subject ID '''

    # initialize the plot
    fig, ax = plt.subplots(figsize=(10,5))

    # build the plot
    ax = sns.swarmplot(x='Percent ID', y='Annotation', data=df, ax=ax)

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
        '-b', '--blast_file',
        help='Please specify the tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--annotation_file',
        help='Please specify the annotation file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')

    blst = args['blast_file']
    anno = args['annotation_file']
    outfile = blst.split('.')[0] + '.png'

    # Read input files and parse the data
    df = parse_data(blst, anno)
    # Build the plot
    contig_fishing_plot(df, outfile)

if __name__ == "__main__":
    main()
