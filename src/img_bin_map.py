#!/usr/bin/env python

"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: September 10, 2021
Last updated: September 12, 2021
License: GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

Purpose:
Create a file mapping IMG, GOLD, original Contig IDs to Bin(s) for downstream use.
"""

import argparse
import logging
import os
from datetime import datetime
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------
def map_gold(scaffold_map, sample_map, bin_map,
    scaffold_samples=False, only_bin=False):
    """
    Create a unified mapping dataframe.
    """

    #Read in the component dataframes
    try:
        logging.info('reading IMG --> GOLD scaffold map to Pandas DF')
        gold_scaffold_df = pd.read_csv(scaffold_map, sep = '\t',
            names = ['Original_Contig_Name', 'IMG_Contig_Name'],
            compression = 'infer')

    except IOError as e:
        logging.exception('could not open {}'.format(scaffold_map))
    #Extract the IMG GOLD analysis ID
    try:
        gold_scaffold_df['GOLD_OID'] = gold_scaffold_df['IMG_Contig_Name'].str.split('_').str[0].str.strip()
    except:
        logging.exception('could not parse GOLD_OID from IMG_Contig_Name')

    try:
        logging.info('reading GOLD --> Sample ID map to Pandas DF')
        gold_sample_df = pd.read_csv(sample_map, sep = '\t', compression = 'infer')
    except IOError as e:
        logging.exception('could not open {}'.format(sample_map))

    try:
        logging.info('reading renamed contig --> Bin --> Sample ID map to Pandas DF')
        bin_map_df = pd.read_csv(bin_map, sep = '\t', compression = 'infer')
    except IOError as e:
        logging.exception('could not open {}'.format(bin_map))

    ###MERGE START
    gold_scaffold_sample_df = pd.merge(gold_scaffold_df, gold_sample_df,
        on = 'GOLD_OID', how = 'left')

    if scaffold_samples:
        #Add Sampling_Site IDs to Original_Contig_Name - used to map to Bin
        og_name = gold_scaffold_sample_df['Sampling_Site'].astype(str) + '_' + gold_scaffold_sample_df['Original_Contig_Name'].astype(str)
        gold_scaffold_sample_df['Original_Contig_Name'] = og_name
    else:
        pass

    bin_gold_sample_df = pd.merge(gold_scaffold_sample_df,
        bin_map_df[['Original_Contig_Name', 'Bin']],
        on = 'Original_Contig_Name', how = 'left')

    if only_bin:
        bin_gold_sample_df = bin_gold_sample_df.dropna(subset=['Bin']).sort_values(by='Bin')
    else:
        pass

    ###MERGE END

    return bin_gold_sample_df
#------------------------------------------------------------------------------
def main():

    script_basename = os.path.basename(__file__)

    logfile_default = '{}_{}.log'.format(os.path.splitext(script_basename)[0],
        str(datetime.now().strftime('%d-%m-%Y_%H-%M-%S')))

    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--output', required=True, type=str, help='Path to output map file. Required.')
    parser.add_argument('-b', '--bin_map', required=True, action='store', type=str, help='bin-contig-sample map file')
    parser.add_argument('-c', '--contig_map', required=True, action='store', type=str, help='IMG contig to GOLD map')
    parser.add_argument('-g', '--gold_sample_map', required=True, help='Mapping file of GOLD IDs to Sample IDs')
    parser.add_argument('-l', '--logfile', required=False, action='store', default=logfile_default, help='path to logfile')

    args = parser.parse_args()

    logging_format = '%(name)s :: %(levelname)s :: %(message)s :: %(asctime)s'

    logging.basicConfig(filename = args.logfile,
        level = logging.DEBUG,
        format = logging_format)

    print('Logfile written to: {}'.format(args.logfile))

    logging.info('SCRIPT: {}'.format(script_basename))

    assert os.path.isfile(args.bin_map)
    assert os.path.isfile(args.contig_map)
    assert os.path.isfile(args.gold_sample_map)

    try:
        bin_map_df = map_gold(scaffold_map = args.contig_map,
            sample_map = args.gold_sample_map,
            bin_map = args.bin_map,
            scaffold_samples=True,
            only_bin=True)
        #Write the mapping dataframe to TSV file
        try:
            bin_map_df.to_csv(args.output, index=False, sep='\t', encoding='utf-8')
            write_success = 'Writing to output tsv: {}'.format(args.output)
            print(write_success)
            logging.info(write_success)
        except IOError as e:
            write_error = 'Unable to write output tsv file {}'.format(args.output)
            print('ERROR: {}'.format(write_error))
            logging.error(write_error)
        finally:
            logging.info('dataframe write try/except block is finished')
    except:
        logging.error('Could not create bin map')
    finally:
        logging.info('read/write try/except block is finished')


#------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
