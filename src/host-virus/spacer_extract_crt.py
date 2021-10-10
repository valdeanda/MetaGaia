#!/usr/bin/env python

"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: August 17, 2021
Last updated: September 10, 2021
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
import pandas as pd
#------------------------------------------------------------------------------
# def map_gold(scaffold_map, sample_map, bin_map, scaffold_samples=False):
#     """
#     Create a unified mapping dataframe from of IMG and GOLD scaffolds to bins.
#     """
#     """
#     Create a unified mapping file for original IMG and Gold scaffolds.
#     """
#     try:
#         logging.info('reading IMG --> GOLD scaffold map to Pandas DF')
#         gold_scaffold_df = pd.read_csv(scaffold_map, sep = '\t',
#             names = ['Original_Contig_Name', 'IMG_Contig_Name'])
#
#     except IOError as e:
#         logging.exception('could not open {}'.format(scaffold_map))
#     #Extract the IMG GOLD analysis ID
#     try:
#         gold_scaffold_df['GOLD_OID'] = gold_scaffold_df['IMG_Contig_Name'].str.split('_').str[0].str.strip()
#         print(gold_scaffold_df.head())
#     except:
#         logging.exception('could not parse GOLD_OID from IMG_Contig_Name')
#
#     try:
#         logging.info('reading GOLD --> Sample ID map to Pandas DF')
#         gold_sample_df = pd.read_csv(sample_map, sep = '\t')
#         print(gold_sample_df.head())
#     except IOError as e:
#         logging.exception('could not open {}'.format(sample_map))
#
#     try:
#         logging.info('reading renamed contig --> Bin --> Sample ID map to Pandas DF')
#         bin_map_df = pd.read_csv(bin_map, sep = '\t')
#     except IOError as e:
#         logging.exception('could not open {}'.format(bin_map))
#
#     #Extract the IMG GOLD analysis ID
#     #gold_scaffold_df['GOLD_OID'] = gold_scaffold_df['IMG_Contig_Name'].str.split('_').str[0].str.strip()
#
#     gold_scaffold_sample_df = pd.merge(gold_scaffold_df, gold_sample_df, on = 'GOLD_OID', how = 'left')
#
#     #gold_sample_bin_df = pd.merge(gold_sample_df, bin_map_df, on = 'SAMPLE_ID')
#     print('gold_scaffold_sample_df')
#     print(gold_scaffold_sample_df.head())
#
#     if scaffold_samples:
#         gold_scaffold_sample_df['Original_Contig_Name'] = gold_scaffold_sample_df['Sampling_Site'].astype(str) + '_' + gold_scaffold_sample_df['Original_Contig_Name'].astype(str)
#     else:
#         pass
#     print('gold_scaffold_sample_df_renamed')
#     print(gold_scaffold_sample_df.head())
#
#
#     bin_gold_sample_df = pd.merge(gold_scaffold_sample_df, bin_map_df,
#         on = 'Original_Contig_Name', how = 'left')
#     print('bin_gold_sample_df')
#     print(bin_gold_sample_df.head())
#
#     return bin_gold_sample_df
#------------------------------------------------------------------------------
def read_crt(crt_file):
    """
    Read CRISPR Recognition Tool output and return a dictionary.
    """
    try:
        crt_handle = open(crt_file, 'r')
        logging.info('opening CRT input...')

        crispr_dict = dict()
        repeat_no = 0
        spacer_no = 0

        for feature in crt_handle:
            feature_list = feature.strip().split()
            gold_contig_id = feature_list[0]
            crispr_no = feature_list[1]
            repeat_seq = feature_list[3]
            spacer_seq = feature_list[4]

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

            crispr_dict[crispr_id]['gold_contig_id'] = gold_contig_id
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
    Filter out poor-quality CRISPRs, return a dictionary of retained CRISPR info.
    """

    logging.info('filtering poor-quality CRISPRs:')
    logging.info('miminum spacer length: {}'.format(min_spc))
    logging.info('minimum repeat length: {}'.format(min_rep))
    logging.info('minimum number of repeats: {}'.format(n_rep))

    rm_crispr = []
    spurious = 'Output may be sourced from spurious CRISPRs.'
    if min_rep < min_rep_warn:
        logging.warning('WARNING: a minimum repeat length >= {} is recommended. {}'.format(min_rep_warn, spurious))
    if n_rep < n_rep_warn:
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
    if rm_crispr:
        rm_crispr = set(rm_crispr)
        logging.info('{} CRISPR did not pass quality check'.format(len(rm_crispr)))

        for rc in rm_crispr:
            crispr_dict.pop(rc, None)
    else:
        logging.info('no CRISPRs will be removed due to quality check')

    return crispr_dict
#------------------------------------------------------------------------------
def img_bin_map(img_map, crispr_contigs):
    """
    Return a dictionary of contig IDs mapped to a Bin/MAG
    containing a CRISPR array.
    """
    try:
        bm_df = pd.read_csv(img_map, sep = '\t', compression = 'infer')
    except:
        logging.error('could not read IMG bin map input file {}'.format(img_map))

    bm_df['bin_contig'] = bm_df['Bin'] + '|' + bm_df['Original_Contig_Name']
    bm_df = bm_df[bm_df['IMG_Contig_Name'].isin(crispr_contigs)]

    bm_dict = bm_df[['bin_contig', 'IMG_Contig_Name']].set_index('IMG_Contig_Name').T.to_dict('list')

    return bm_dict
#------------------------------------------------------------------------------
def write_spacer_fasta(crispr_dict, spacer_fasta, id_change=None):
    """
    Write an output FASTA file of CRISPR spacer sequences.
    """
    if id_change:
        changed = []
    else:
        pass

    try:
        fasta_handle = open(spacer_fasta, 'w')
        nspacers = 0

        for crispr, feature in crispr_dict.items():
            if id_change:
                gold_oid = crispr.split('|')[0]
                if gold_oid in id_change.keys():
                    crispr = '{}|{}'.format(id_change[gold_oid][0], crispr)
                    changed.append(crispr)
                else:
                    pass
            else:
                pass

            if isinstance(feature, dict):
                spacers = [s for s in feature.keys() if s.startswith('spacer')]
                for spacer in spacers:
                    spacer_header = '>{}|{}'.format(crispr, spacer)
                    spacer_entry = '{}\n{}\n'.format(spacer_header, feature[spacer])
                    fasta_handle.write(spacer_entry)
                nspacers += len(spacers)
            else:
                logging.error('ERROR: CRISPR dictionary does not contain nested feature dict')
                logging.error('Cannot write to output FASTA, exiting...')
                sys.exit(1)

        spacer_msg = '{} spacers written to: {}'.format(nspacers, spacer_fasta)

        print(spacer_msg)
        logging.info(spacer_msg)
        
        if changed:
            logging.info('{} CRISPR IDs changed'.format(len(changed)))
        else:
            pass



        fasta_handle.close()


    except IOError:
        logging.error('ERROR: could not open output FASTA file {}'.format(spacer_fasta))

    return None
#------------------------------------------------------------------------------
def main():

    logfile_default = 'spacer_extract_crt_{}.log'.format(str(datetime.now().strftime('%d-%m-%Y_%H-%M-%S')))

    #epilog = 'spacer_extract_crt.py -i IMG CRT CRISPR file'
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--crt_file', required=True, type=str, help='Input file in tsv format containing CRISPR Recognition Tool output. Required.')
    parser.add_argument('-o', '--output_fasta', required=True, type=str, help='Path to output FASTA file containing CRISPR spacer sequences. Required.')
    parser.add_argument('-s', '--min_spacer_length', required=False, type=int, default=23, help='minimum length for a CRISPR spacer. Default = 23')
    parser.add_argument('-r', '--min_repeat_length', required=False, type=int, default=11, help='minimum length for a CRISPR repeat. Default = 11')
    parser.add_argument('-m', '--min_repeats', required=False, type=int, default=3, help='minimum number of repeats required to retain CRISPR. Default = 3')
    parser.add_argument('-q', '--quality_off', required=False, action='store_true', help='use this option to skip quality control of CRISPRs')
    parser.add_argument('-c', '--img_bin_map', required=False, action='store', type=str, help='bin map output from img_bin_map.py')
    parser.add_argument('-l', '--logfile', required=False, action='store', default=logfile_default, help='path to logfile')
    args = parser.parse_args()

    logging_format = '%(name)s :: %(levelname)s :: %(message)s :: %(asctime)s'

    logging.basicConfig(filename = args.logfile,
        level = logging.DEBUG,
        format = logging_format)

    print('Logfile written to: {}'.format(args.logfile))

    logging.info('SCRIPT: {}'.format(os.path.basename(__file__)))

    crispr_dict = dict()


    crispr_dict = read_crt(args.crt_file)

    fasta_write_exception = 'could not write output spacer FASTA {}'.format(args.output_fasta)


    if crispr_dict:
        bcd = dict()
        if args.img_bin_map:
            try:
                gold_contigs_crispr = list(set([crispr_dict[key]['gold_contig_id'] for key in crispr_dict.keys()]))
                bcd = img_bin_map(args.img_bin_map, gold_contigs_crispr)

            except:
                logging.warning('could not get crispr bin ids')
            finally:
                logging.info('img bin-crispr try/except block is finished')

        else:
            pass


        ncrispr = len(crispr_dict.keys())
        if not args.quality_off:
            try:
                logging.info('CRISPR quality control')
                crispr_dict_filter = crispr_quality(crispr_dict, min_spc = args.min_spacer_length,
                    min_rep = args.min_repeat_length, n_rep = args.min_repeats)

                ncrispr_retained = len(crispr_dict_filter.keys())
                logging.info('{} of {} CRISPRs retained'.format(ncrispr_retained, ncrispr))

            except:
                logging.exception('could not filter CRISPR dictionary')

            #Write to output FASTA
            try:
                logging.info('writing CRISPR spacers to output FASTA: {}'.format(args.output_fasta))
                write_spacer_fasta(crispr_dict = crispr_dict_filter,
                spacer_fasta = args.output_fasta,
                id_change = bcd)
            except:
                logging.exception(fasta_write_exception)

        else:
            logging.info('skipping CRISPR quality control')
            try:
                write_spacer_fasta(crispr_dict = crispr_dict,
                spacer_fasta = args.output_fasta,
                id_change = bcd)
            except:
                logging.exception(fasta_write_exception)

        print('Finished!')
        logging.info('Finished')
    else:
        logging.error('CRISPR dictionary is empty!')
        sys.exit(1)
#------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
