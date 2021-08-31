#!/usr/bin/env python

"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: August 17, 2021
Last updated: August 30, 2021
License: GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

Purpose:
Parse CRISPR Recognition Tool (CRT) output from IMG 5.0
Perform optional CRISPR quality filtering
Output a FASTA file of spacer sequences.
"""
import argparse
import os
import sys
import logging
from datetime import datetime
#------------------------------------------------------------------------------
# def read_bin_map(bin_map):
#     """
#     Read in the bin - scaffold - sample map file.
#     """
#     bin_map_dict = dict()
#     try:
#         bin_map_handle = open(bin_map, 'r')
#         bin_map_handle.readline()
#         for i in bin_map_handle:
#             i = i.strip()
#             i_list = i.split('\t')
#             original_contig_name = i_list[0]
#             bin = i_list[1]
#             sampling_site = i_list[2]
#
#             bin_map_dict[original_contig_name] = {'bin':bin,
#                 'sampling_site':sampling_site}
#
#         return bin_map_dict
#
#         bin_map_handle.close()
#     except:
#         logging.exception('could not open {}'.format(bin_map))
#         return None
#------------------------------------------------------------------------------
# def read_contig_map(contig_map):
#     contig_map_dict = dict()
#     try:
#         contig_map_handle = open(contig_map, 'r')
#         for i in contig_map_handle:
#             i = i.strip()
#             i_list = i.split('\t')
#             img_contig_id = i_list[0]
#             gold_contig_id = i_list[1]
#             gold_oid = gold_contig_id.split('_')[0]
#             contig_map_dict[gold_contig_id] = {'img_contig_id':img_contig_id,
#             'gold_oid':gold_oid}
#         return contig_map_dict
#         contig_map_handle.close()
#     except:
#         logging.exception('could not read {}'.format(contig_map))
#         return None
# #------------------------------------------------------------------------------
# def read_gold_map(gold_map):
#     gold_map_dict = dict()
#     try:
#         gold_map_handle = open(gold_map, 'r')
#         for i in gold_map_handle:
#             i = i.strip()
#             i_list = i.split('\t')
#             gold_oid = i_list[0]
#             sample_id = i_list[1]
#             gold_dict[gold_oid] = sample_id
#
#         return gold_map_dict
#         gold_map_handle.close()
#     except:
#         logging.exception('could not read {}'.format(gold_map))
#         return None
#------------------------------------------------------------------------------
# def unified_map(bin_map, contig_map, gold_map):
#     bin_contig_dict = dict()
#     for key, value in contig_map:
#         if key in
#
#         sample = bin_contig_dict[gold_oid]


def map_gold(scaffold_map, sample_map, scaffold_samples=False):
    import pandas as pd
    """
    Create a unified mapping file for original IMG and Gold scaffolds.
    """
    try:
        print('reading IMG --> GOLD scaffold map to Pandas DF')
        gold_scaffold_map = pd.read_csv(scaffold_map, sep = '\t',
            names = ['Original_Contig_Name', 'IMG_Contig_Name'])
    except IOError as e:
        print('ERROR: could not open {}'.format(scaffold_map))

    try:
        print('reading GOLD --> Sample ID map to Pandas DF')
        gold_sample_map = pd.read_csv(sample_map, sep = '\t')
    except IOError as e:
        print('ERROR: could not open {}'.format(sample_map))

    #Extract the IMG GOLD analysis ID
    gold_scaffold_map['GOLD_OID'] = gold_scaffold_map['IMG_Contig_Name'].str.split('_').str[0].str.strip()

    gold_scaffold_sample_df = pd.merge(gold_scaffold_map, gold_sample_map, on = 'GOLD_OID', how = 'right')

    if scaffold_samples:
        gold_scaffold_sample_df['Original_Contig_Name'] = gold_scaffold_sample_df['SAMPLE_ID'].astype(str) + '_' + gold_scaffold_sample_df['Original_Contig_Name'].astype(str)
    else:
        pass

    return gold_scaffold_sample_df
#------------------------------------------------------------------------------
def map_bins(scaffold_id, bin_map_dict):
    """
    Map bins to scaffold IDs
    """
    if scaffold_id in bin_map_dict:
        bin = bin_map_dict[scaffold_id]['bin']

    else:
        bin = 'NoBin'
    return bin
#------------------------------------------------------------------------------
def read_crt(crt_file, bin_map_dict=None):
    """
    Read CRISPR Recognition Tool output and return a dictionary.
    """
    try:
        crt_handle = open(crt_file, 'r')
        logging.info('opening CRT input...')

        crispr_dict = {}
        repeat_no = 0
        spacer_no = 0

        if bin_map_dict:
            logging.info('Bin IDs will be mapped to contigs')
        else:
            pass

        for feature in crt_handle:
            feature_list = feature.strip().split()
            gold_contig_id = feature_list[0]
            crispr_no = feature_list[1]
            repeat_seq = feature_list[3]
            spacer_seq = feature_list[4]


            if bin_map_dict:
                bin = map_bins(gold_contig_id, bin_map_dict)
                if bin:
                    crispr_id = '{}|{}|CRISPR-{}'.format(bin, gold_contig_id, crispr_no)
            else:
                crispr_id = '{}|CRISPR-{}'.format(gold_contig_id, crispr_no)


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

        crt_handle.close()


        return crispr_dict

    except IOError as e:
        logging.exception('cannot open file {}'.format(crt_file))
        print('ERROR: cannot open file {}'.format(crt_file))
        return None
#------------------------------------------------------------------------------
def crispr_quality(crispr_dict, min_spc, min_rep, n_rep, min_rep_warn = 7, n_rep_warn = 2):
    """
    Filter out poor-quality CRISPRs and return a dictionary.
    """

    logging.info('filtering poor-quality CRISPRs:')
    logging.info('miminum spacer length: {}'.format(min_spc))
    logging.info('minimum repeat length: {}'.format(min_rep))
    logging.info('minimum number of repeats: {}'.format(n_rep))

    rm_crispr = []
    spurious = 'Output may be sourced from spurious CRISPRs.'
    if min_rep < min_rep_warn:
        #print('WARNING: a minimum repeat length >= {} is recommended. {}'.format(min_rep_warn, spurious))
        logging.warning('WARNING: a minimum repeat length >= {} is recommended. {}'.format(min_rep_warn, spurious))
    if n_rep < n_rep_warn:
        #print('WARNING: CRISPRs should contain >= {} repeats. {}'.format(n_rep_warn, spurious))
        logging.warning('WARNING: CRISPRs should contain >= {} repeats. {}'.format(n_rep_warn, spurious))

    for crispr, feature in crispr_dict.items():
        if isinstance(feature, dict):
            repeats = [r for r in feature.keys() if r.startswith('repeat')]
            spacers = [s for s in feature.keys() if s.startswith('spacer')]

            if len(repeats) < n_rep:
                logging.info('{} will be ignored; number of repeats < {}'.format(crispr, n_rep))
                rm_crispr.append(crispr)

            elif any([len(feature[spacer]) < min_spc for spacer in spacers]):
                logging.info('{} will be ignored; contains spacer(s) of length < {} bp'.format(crispr, min_spc))
                rm_crispr.append(crispr)

            elif any([len(feature[repeat]) < min_rep for repeat in repeats]):
                logging.info('{} will be ignored; contains repeat(s) of length < {} bp'.format(crispr, min_rep))
                rm_crispr.append(crispr)
            else:
                pass

    rm_crispr = set(rm_crispr)
    logging.info('{} CRISPR did not pass quality check'.format(len(rm_crispr)))

    for rc in rm_crispr:
        crispr_dict.pop(rc, None)

    return crispr_dict

#------------------------------------------------------------------------------
def write_spacer_fasta(crispr_dict, spacer_fasta):
    """
    Write an output FASTA file of CRISPR spacer sequences.
    """
    try:
        fasta_handle = open(spacer_fasta, 'w')
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
                logging.error('ERROR: CRISPR dictionary does not contain nested feature dict')
                logging.warning('Cannot write to output FASTA, exiting...')
                sys.exit(1)

        spacer_msg = '{} spacers written to {}'.format(nspacers, spacer_fasta)

        print(spacer_msg)
        logging.info(spacer_msg)

        fasta_handle.close()

    except IOError:
        logging.error('ERROR: could not open output FASTA file {}'.format(spacer_fasta))

    return None
#------------------------------------------------------------------------------
def main():

    logfile_default = 'spacer_extract_crt_{}.log'.format(str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S')))

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--crt_file', required=True, type=str, help='Input file in tsv format containing CRISPR Recognition Tool output')
    parser.add_argument('-o', '--output_fasta', required=True, type=str, help='Output FASTA file containing CRISPR spacer sequences')
    parser.add_argument('-s', '--min_spacer_length', required=False, type=int, default=23, help='minimum length for a CRISPR spacer. Default = 23')
    parser.add_argument('-r', '--min_repeat_length', required=False, type=int, default=11, help='minimum length for a CRISPR repeat. Default = 11')
    parser.add_argument('-m', '--min_repeats', required=False, type=int, default=3, help='minimum number of repeats required to retain CRISPR. Default = 3')
    parser.add_argument('-q', '--quality_off', required=False, action='store_true', help='use this option to skip quality control of CRISPRs')
    parser.add_argument('-b', '--bin_map', required=False, action='store', type=str, help='bin-contig-sample map file')
    parser.add_argument('-c', '--contig_map', required=False, action='store', type=str, help='IMG contig')
    parser.add_argument('-l', '--logfile', required=False, action='store', default=logfile_default, help='path to logfile')
    args = parser.parse_args()

    logging_format = '%(name)s :: %(levelname)s :: %(message)s :: %(asctime)s'

    print('Logfile written to: {}'.format(args.logfile))

    logging.basicConfig(filename = args.logfile,
        level = logging.DEBUG,
        format = logging_format)

    if args.bin_map:
        logging.info('reading in scaffold - bin - sample map')
        bin_map_dict = read_bin_map(args.bin_map)
        crispr_dict = read_crt(args.crt_file, bin_map_dict = args.bin_map)
    else:
        crispr_dict = read_crt(args.crt_file)

    fasta_write_exception = 'could not write output spacer FASTA {}'.format(args.output_fasta)

    if isinstance(crispr_dict, dict):
        ncrispr = len(crispr_dict.keys())
        if not args.quality_off:
            try:
                crispr_dict_filter = crispr_quality(crispr_dict, min_spc = args.min_spacer_length,
                    min_rep = args.min_repeat_length, n_rep = args.min_repeats)

                ncrispr_retained = len(crispr_dict_filter.keys())
                logging.info('{} of {} CRISPRs retained'.format(ncrispr_retained, ncrispr))

            except:
                logging.exception('could not filter CRISPR dictionary')

            #Write to output FASTA
            try:
                logging.info('writing CRISPR spacers to output FASTA {}'.format(args.output_fasta))
                write_spacer_fasta(crispr_dict_filter, args.output_fasta)
            except:
                logging.exception(fasta_write_exception)

        else:
            try:
                write_spacer_fasta(crispr_dict, args.output_fasta)
            except:
                logging.exception(fasta_write_exception)
        print('Finished!')
        logging.info('Finished')
    else:
        logging.error('invalid CRISPR dictionary!')
        sys.exit(1)
#------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
