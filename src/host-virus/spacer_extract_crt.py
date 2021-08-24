#!/usr/bin/env python

"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: August 17, 2021
Last updated: August 23, 2021

Purpose: Parse CRISPR Recognition Tool (CRT) output from IMG 5.0 and
output a FASTA file of spacer sequences.
"""
import argparse
import os
import sys
#------------------------------------------------------------------------------
def read_crt(crt_file):
    """ Read CRISPR Recognition Tool output and return a dictionary."""
    try:
        crt_handle = open(crt_file, 'r')
    except FileNotFoundError:
        print('cannot open file {}'.format(crt_file))

    crispr_dict = {}
    repeat_no = 0
    spacer_no = 0

    for feature in crt_handle:
        feature_list = feature.strip().split()
        contig_id = feature_list[0]
        crispr_no = feature_list[1]
        repeat_seq = feature_list[3]
        spacer_seq = feature_list[4]


        crispr_id = '{}|CRISPR-{}'.format(contig_id, crispr_no)

        if crispr_id in crispr_dict.keys():
            repeat_no += 1
            spacer_no += 1

        else:
            crispr_dict[crispr_id] = {}
            repeat_no = 0
            spacer_no = 0


        repeat_id = 'repeat-{}'.format(repeat_no)
        spacer_id = 'spacer-{}'.format(spacer_no)

        crispr_dict[crispr_id][repeat_id] = repeat_seq
        #Exclude entries where the spacer is missing
        if not spacer_seq == 'c':
            crispr_dict[crispr_id][spacer_id] = spacer_seq
        else:
            pass

    return crispr_dict

#------------------------------------------------------------------------------
def crispr_quality(crispr_dict, min_spc, min_rep, n_rep, min_rep_warn = 7, n_rep_warn = 2):
    """ Filter out poor-quality CRISPRs """

    rm_crispr = []
    spurious = 'Output may be sourced from spurious CRISPRs.'
    if min_rep < min_rep_warn:
        print('WARNING: a minimum repeat length >= {} is recommended. {}'.format(min_rep_warn, spurious))
    if n_rep < n_rep_warn:
        print('WARNING: CRISPRs should contain >= {} repeats. {}'.format(n_rep_warn, spurious))

    for crispr, feature in crispr_dict.items():
        if isinstance(feature, dict):
            repeats = [r for r in feature.keys() if r.startswith('repeat')]
            spacers = [s for s in feature.keys() if s.startswith('spacer')]

            if len(repeats) < n_rep:
                print('{} will be ignored; number of repeats < {}'.format(crispr, n_rep))
                rm_crispr.append(crispr)

            elif any([len(feature[spacer]) < min_spc for spacer in spacers]):
                print('{} will be ignored; contains spacer(s) of length < {} bp'.format(crispr, min_spc))
                rm_crispr.append(crispr)

            elif any([len(feature[repeat]) < min_rep for repeat in repeats]):
                print('{} will be ignored; contains repeat(s) of length < {} bp'.format(crispr, min_rep))
                rm_crispr.append(crispr)
            else:
                pass

    rm_crispr = set(rm_crispr)
    print('{} CRISPR did not pass quality check'.format(len(rm_crispr)))

    for rc in rm_crispr:
        crispr_dict.pop(rc, None)

    return crispr_dict

#------------------------------------------------------------------------------
def write_spacer_fasta(crispr_dict, spacer_fasta):
    """ Write an output FASTA file of CRISPR spacer sequences """
    try:
        fasta_handle = open(spacer_fasta, 'w')
    except IOError:
        print('could not open output FASTA file {}'.format(spacer_fasta))

    nspacers = 0

    for crispr, feature in crispr_dict.items():
        if isinstance(feature, dict):
            spacers = [s for s in feature.keys() if s.startswith('spacer')]
            for spacer in spacers:
                spacer_header = '>{}|{}'.format(crispr, spacer)
                spacer_entry = '{}\n{}\n'.format(spacer_header, feature[spacer])
                fasta_handle.write(spacer_entry)
            nspacers += len(spacers)
        else:
            print('ERROR - CRISPR dictionary does not contain nested feature dict')
            print('Cannot write to output FASTA, exiting...')
            sys.exit(1)
    print('{} spacers written to {}'.format(nspacers, spacer_fasta))

    return None

#------------------------------------------------------------------------------
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--crt_file', required=True, type=str, help='Input file in tsv format containing CRISPR Recognition Tool output')
    parser.add_argument('-o', '--output_fasta', required=True, type=str, help='Output FASTA file containing CRISPR spacer sequences')
    parser.add_argument('-s', '--min_spacer_length', required=False, type=int, default=23, help='minimum length for a CRISPR spacer. Default = 23')
    parser.add_argument('-r', '--min_repeat_length', required=False, type=int, default=11, help='minimum length for a CRISPR repeat. Default = 11')
    parser.add_argument('-m', '--min_repeats', required=False, type=int, default=3, help='minimum number of repeats required to retain CRISPR. Default = 3')
    parser.add_argument('-q', '--quality_off', required=False, action='store_true', help='use this option to skip quality control of CRISPRs')
    args = parser.parse_args()

    print('reading CRT input...')
    crispr_dict = read_crt(args.crt_file)
    ncrispr = len(crispr_dict.keys())
    print('{} CRISPR arrays identified'.format(ncrispr))
    if not args.quality_off:
        print('filtering poor-quality CRISPRs:')
        print('miminum spacer length: {}'.format(args.min_spacer_length))
        print('minimum repeat length: {}'.format(args.min_repeat_length))
        print('minimum number of repeats: {}'.format(args.min_repeats))
        crispr_dict_filter = crispr_quality(crispr_dict, min_spc = args.min_spacer_length,
        min_rep = args.min_repeat_length, n_rep = args.min_repeats)


        ncrispr_retained = len(crispr_dict_filter.keys())
        print('{} of {} CRISPRs retained'.format(ncrispr_retained, ncrispr))

        print('writing CRISPR spacers to output FASTA')
        write_spacer_fasta(crispr_dict_filter, args.output_fasta)


    else:
        print('--quality_off option selected: CRISPR quality control will not be performed')
        print('writing CRISPR spacers to output FASTA')
        write_spacer_fasta(crispr_dict, args.output_fasta)
    print('Finished!')
#------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
