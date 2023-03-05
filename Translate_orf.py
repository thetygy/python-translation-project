#! /usr/bin/env python3

import re
import sys
from find_orf import find_first_orf, parse_sequence_from_path
from translate import translate_sequence

def main():
    import argparse

    # Create a command-line parser object
    parser = argparse.ArgumentParser()

    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']

    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A start codon. This option can be used multiple times '
                    'if there are multiple start codons. '
                    'Default: {0}.'.format(" ".join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A stop codon. This option can be used multiple times '
                    'if there are multiple stop codons. '
                    'Default: {0}.'.format(" ".join(default_stop_codons))))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    # Check to see if start/stop codons were provided by the caller. If not,
    # use the defaults.
    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons

    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)
    sys.stdout.write('{}\n'.format(orf))



first_orf=find_first_orf(sequence)
translate_sequence(first_orf)


