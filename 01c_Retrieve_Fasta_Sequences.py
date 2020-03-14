#!/usr/bin/env python

'''Retrieve Fasta Sequences for list from assembled contigs fasta file

This script takes a list of contig sequence names and returns the 
sequences in fasta format.

This tool takes the following input parameters:

    * contig_list_file : str - file_name.txt
        File containing list of contig names to retrieve the fasta
        sequence for. Should contain one sequence name per line.
    * fasta_assembly_file : str - file_name.fasta
        File containing the fasta sequence for the contig names in 
        cntg_list.
    * output_file_name : str - file_name.fasta
        The name for the output fasta formatted file.

This script returns the following files:

    * output_file_name.fasta
        Fasta file containing the fasta sequences matching the sequence
        names given in the contig_list_file.

This script requires the following packages:

    * argparse

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 13th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def read_fasta(fp):
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


def retrieve_fasta_sequence(cntg_list, assembly, outfile):
    """Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    cntg_list : str - file_name.txt
        File containing list of contig names to retrieve the fasta
        sequence for. Should contain one sequence name per line.
    assembly : str - file_name.fasta
        File containing the fasta sequence for the contig names in 
        cntg_list.

    Returns
    -------
    True - writes two files.
        Lines passing the filter to written to infile.fltrdBstHts.blst
        and a list of subject sequence (contig) names is written to
        infile.
    """

    cntgs = {}

    with open(cntg_list, 'r') as file:
        for line in file:
            cntg = line.rstrip()
            if line.startswith('>'): cntgs[cntg] = ''
            else: cntgs[f'>{cntg}'] = ''

    with open(assembly, 'r') as fa, open(outfile, 'w') as o:
        for name, seq in read_fasta(fa):
            if name in cntgs:
                o.write(f'{name}\n{seq}')


    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-c', '--contig_list_file',
        help=
            'Specify a text file containing a list of contig sequence names.'
            'There should be one contig name per line.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--fasta_assembly_file',
        help=
            'Specify the metagenome or draft genome assembly file '
            'in fasta format that contains the contig sequences to extract.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help=
            'What would you like to call the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    retrieve_fasta_sequence(
                    args['contig_list_file'],
                    args['fasta_assembly_file']
                    args['output_file_name']
                    )


if __name__ == "__main__":
    main()
