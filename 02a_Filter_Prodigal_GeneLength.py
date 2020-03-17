#!/usr/local/pacerepov1/python/2.7/bin/python

''' Filter fasta file by a minimum sequence length

This script takes an input file of nucleotide or amino acid sequence
and removes all sequences below a user specified chacter length or a 
default of 300 characters. Characters are base pairs if using nucleotide
sequence or amino acids if using protein sequence.

Outputs the filtered file as input_file.lenfltr

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
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def Fasta_rename_sequences(infile, len_filter):

    outfile = infile + '.lenfltr'

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        i = 0
        c = 0
        for name, seq in read_fasta(f):
            c += 1
            if len(seq) >= len_filter:
                o.write(f'{name}\n{seq}\n')
                i += 1

    print(
        f'\nRemoved {c-i} sequences shorter than {len_filter} characters '
        f'out of {c} total sequences.\n {i} Sequences passing filter written '
        f'to file: {outfile}. Happy sciencing super sciencer!\n\n'
        )


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Specify a the input fasta file name.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--minimum_sequence_length',
        help=
            '(Optional) Specify the minimum sequence length to return. '
            '(Default = 300 characters).',
        metavar='',
        type=int,
        required=False,
        default=300
        )
    args=vars(parser.parse_args())

    print('\n\nRunning Script ...')
    # Do what you came here to do
    Fasta_rename_sequences(args['input_file'], args['minimum_sequence_length'])



if __name__ == "__main__":
    main()

