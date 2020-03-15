#!/usr/bin/env python

'''Best Hit Filter for Tabular Blast Output.

This script filters tabular blast output for best hit based on bitscore,
as well as user defined percent match length, and percent identity.

This tool takes the following input parameters:

    * tabular blast input file that includes query and subject lengths
      in columns 13 and 14.
    * percent_match_length to filter for in base pairs (bps) (ex: 70.0).
    * percent_identity to filter for (ex: 70.0).
    * minimum_subject_length to filter for in bps (ex: 3000).

This script returns the following files:

    * input_file.fltrdBstHts.blst
        Contains entries from original file that passed the filters.
    * input_file.fltrdBstHtsLst.txt
        Contains list of subject sequence names (contigs) for entries
        that passed the filter.
    * input_file.metatemplate.txt
        Contains header and gene names for meta file to use for 01b

This script requires the following packages:

    * argparse
    * collections.defaultdict

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
from collections import defaultdict


def best_hits(query_subject, bitscore, d, line, dups):
    """ Filters the besthit based on bitscore """

    v = ''

    if query_subject in d:
        v = f', DUPLICATE {len(query_subject)}'
        dups += 1
        old_bitscore = float(d[query_subject].split('\t')[11])

        if bitscore > old_bitscore:
            v += '- KEPT - BETTER HIT THAN PREVIOUS'
            d[query_subject] = [line]

        else: v += '- REMOVED - LOWER HIT THAN PREVIOUS'

    else:
        d[query_subject] = [line]

    return d, dups, v


def filter_blast(infile, pml, pid, msl, verbose):
    """Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    infile : str - file_name.txt
        Tabular blast output file.
    pml : float - ex: 70.0
        Minimum Percent Match Length to keep
        (alignment length / query sequence length).
    pid : float - ex: 70.0
        Minimum Percent Sequence Identity to keep.
    msl : int - ex: 3000
        Minimum length of subject sequence (contig) to keep.

    Returns
    -------
    True - writes file.
        writes out tabular blast lines passing the filter to 
        a new file infile.fltrdBstHts.blst
    """

    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    low_pml = 0 # counter for number of matches failing filters
    low_pid = 0 # 
    low_msl = 0
    passes = 0 # counter for number of matches passing filters
    total = 0 # counter for total blast entries in file

    if verbose:
        print('query, subject, pid, aLen, qLen, pml, sLen, decision')

    with open(infile, 'r') as f:
        for l in f:
            total += 1
            X = l.rstrip().split('\t')
            query = X[0]
            subject = X[1] 
            bitscore = float(X[11]) # bitscore
            pIdent = float(X[2]) # percent sequence identity
            aLen = int(X[3]) # read alignment length
            qLen = int(X[12]) # full length of read
            sLen = int(X[13])
            pMatch = aLen / qLen # percent match length of read length

            verbose_line = (
                f'{query}, {subject}, {pIdent}, {aLen}, {qLen}, '
                f'{pMatch*100:.2f}, {sLen}'
                )

            if pMatch < (pml/100):
                low_pml += 1
                verbose_line += ', FAILED pml'
            elif pIdent < pid:
                low_pid += 1
                verbose_line += ', FAILED pid'
            elif sLen < msl:
                low_msl += 1
                verbose_line += ', FAILED msl'
            elif pMatch >= (pml/100) and pIdent >= pid and sLen >= msl:
                verbose_line += ', PASSED'
                qs = f'{query}_{subject}'
                d, dups, v = best_hits(qs, bitscore, d, l, dups)
                passes += 1
                verbose_line += v
            else:
                print(
                    '\n\nSomething is wrong with line:\n'
                    f'{l}\n'
                    )

            if verbose: print(verbose_line)

    ibase = infile.split('.')[0]
    query_sequence_set = defaultdict(list)
    outfile_blast = ibase + '.fltrdBstHts.blst'
    outfile_meta = ibase + '.metatemplate.txt'
    with open(outfile_blast, 'w') as ob:
        for k,v in d.items():
            ob.write(v[0])
            X = v[0].split('\t')
            subject_sequence_name = X[1]
            query_sequence_name = X[0]
            query_sequence_set[query_sequence_name].append(subject_sequence_name)

    with open(outfile_meta, 'w') as om:
        header = 'query name, annotation, color\n'
        om.write(header)
        for query, contigs in query_sequence_set.items():
            om.write(f'{query}, \n')
            outfile_list = ibase + f'_{query}.fltrdBstHtsLst.txt'
            with open(outfile_list, 'a') as ol:
                for contig in contigs:
                    ol.write(contig + '\n')

    print('\nTotal number of entries in blast file:', total)
    print('Number of matches lower than pml:', low_pml)
    print('Number of matches >= pml but lower than pid:', low_pid)
    print('Number of matches >= pml & pid but lower than msl:', low_msl)
    print('Number of matchs  passing the filters:', passes)
    print(
        'Number of passing matches with multiple hits to same contig '
        'that were removed:', dups
        )
    print(
        'Number of best hit entries passing the filter and written to '
        'new file:', len(d)
        )

    if len(d) == 0:
        print(
            '\n\nNo viable sequence alignments passing the filters.\n'
            'Try changing some parameters or go fish (for more contigs).\n'
            'You can use the -v option to see which filter removes the match.'
            )

    print('\n\n')

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the tabular blast file with added sequence length',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pml', '--percent_match_length',
        help=
            'Minimum Percent Match Length (in bps) to keep '
            '(alignment length / query sequence length). . . . '
            '(Default = 50)',
        metavar='',
        type=float,
        required=False,
        default=50
        )
    parser.add_argument(
        '-pid', '--percent_sequence_identity',
        help=
            'Minimum percent sequence identity to keep. . . . '
            '(Default = 70)',
        metavar='',
        type=float,
        required=False,
        default=70
        )
    parser.add_argument(
        '-msl', '--minimum_subject_length',
        help=
            'Minimum length (in bps) of subject (contig)                '
            'sequence to keep (Default = 1000).',
        metavar='',
        type=int,
        required=False,
        default=1000
        )
    parser.add_argument(
        '-v', '--verbose',
        help='Print decision for each match to screen.',
        required=False,
        action='store_true'
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')
    filter_blast(
        args['input_file'],
        args['percent_match_length'],
        args['percent_sequence_identity'],
        args['minimum_subject_length'],
        args['verbose']
        )


if __name__ == "__main__":
    main()
