"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: August 23, 2021
Last updated: September 12, 2021

Purpose: Parse IMG 5.0 GFF3 output for MetaGaia metabolic analysis.
"""

import argparse
import logging
import pandas as pd
import os
import sys
from datetime import datetime
#------------------------------------------------------------------------------
def gff3_to_df(gff_file, attribute_select, format='img5'):
    """
    Parse a GFF3 file and return a Pandas dataframe.
    """

    valid_formats = ['standard', 'img5']
    format_error_msg = 'ERROR - gff3_to_df() requires a valid format, choose from: {}'.format(' '.join(valid_formats))
    gff_data = []

    try:
        gff_handle = open(gff_file, 'r', encoding = 'utf-8')
    except IOError as e:
        print('Could not open {}'.format(gff_file))

    if not format in valid_formats:
        print(format_error_msg)
        exit(1)
    else:
        print('Parsing {}, expecting {} GFF3 format...'.format(gff_file, format))

    for record in gff_handle:
        record_dict = {}
        if not record.startswith('#'):
            record = record.rstrip()
            record_list = record.split('\t')
            if len(record_list) < 9:
                print('Only {} column(s) in {}  - not consistent with GFF3 standard!'.format(len(record_list), gff_file))
                return None
            else:
                seqid = record_list[0]
                source = record_list[1]
                ftype = record_list[2]
                start = record_list[3]
                end = record_list[4]
                score = record_list[5]
                strand = record_list[6]
                phase = record_list[7]
                attributes = record_list[8]
                attributes_list = attributes.split(';')

                attributes_dict = {}
                if attributes_list:
                    for attribute in attributes_list:
                        alist = attribute.split('=')
                        att_key = alist[0]
                        att_value = alist[1]
                        attributes_dict[att_key] = att_value
                else:
                    pass

                if format == 'standard':
                    gff_dict = {'Contig_Name':seqid, 'source':source, 'ftype':ftype,
                    'start':start, 'end':end, 'score':score, 'strand':strand,
                    'phase':phase, 'attributes':attributes_dict}
                elif format == 'img5':
                    record_dict['IMG_Contig_Name'] = seqid
                    record_dict['Gene_Type'] = ftype

                else:
                    print(format_error_msg)

                for key in attribute_select:
                    if key in attributes_dict.keys():
                    #if key in attribute_select:
                        record_dict[key] = attributes_dict[key]

                    else:
                        record_dict[key] = 'NaN'

        gff_data.append(record_dict)

    gff_handle.close()

    print('Converting GFF3 dict to pandas df...')
    gff_df = pd.DataFrame.from_dict(gff_data)

    return gff_df
#------------------------------------------------------------------------------
# def map_gold(scaffold_map, sample_map, scaffold_samples=False):
#     """
#     Create a unified mapping file for original IMG and Gold scaffolds.
#     """
#     try:
#         print('reading IMG --> GOLD scaffold map to Pandas DF')
#         gold_scaffold_map = pd.read_csv(scaffold_map, sep = '\t',
#             names = ['Original_Contig_Name', 'IMG_Contig_Name'])
#     except IOError as e:
#         print('ERROR: could not open {}'.format(scaffold_map))
#
#     try:
#         print('reading GOLD --> Sample ID map to Pandas DF')
#         gold_sample_map = pd.read_csv(sample_map, sep = '\t')
#     except IOError as e:
#         print('ERROR: could not open {}'.format(sample_map))
#
#     #Extract the IMG GOLD analysis ID
#     gold_scaffold_map['GOLD_OID'] = gold_scaffold_map['IMG_Contig_Name'].str.split('_').str[0].str.strip()
#
#     gold_scaffold_sample_df = pd.merge(gold_scaffold_map, gold_sample_map, on = 'GOLD_OID', how = 'right')
#
#     if scaffold_samples:
#         gold_scaffold_sample_df['Original_Contig_Name'] = gold_scaffold_sample_df['SAMPLE_ID'].astype(str) + '_' + gold_scaffold_sample_df['Original_Contig_Name'].astype(str)
#     else:
#         pass
#
#     return gold_scaffold_sample_df
#------------------------------------------------------------------------------
def gff_pd_multiread(path_file, attr_sel, colsep = '\t', format = 'img5'):
    """
    Read in a file with paths of multiple GFF files,
    read them in as pandas dataframes, and concatenate them.
    """
    pd_df_list = []

    try:
        path_handle = open(path_file, 'r')
    except IOError as e:
        print('could not open path file {}')

    for path in path_handle:
        path = path.strip()
        if os.path.isfile(os.path.abspath(path)):
            try:
                df = gff3_to_df(gff_file = path,
                    attribute_select = attr_sel,
                    format = format)

                pd_df_list.append(df)

            except IOError as e:
                print('WARNING - could not convert GFF file {} to pandas dataframe'.format(path))

    path_handle.close()

    gff_df_concat = pd.concat(pd_df_list)

    return gff_df_concat
#------------------------------------------------------------------------------
def main():

    script_basename = os.path.basename(__file__)

    logfile_default = '{}_{}.log'.format(os.path.splitext(script_basename)[0],
        str(datetime.now().strftime('%d-%m-%Y_%H-%M-%S')))

    field_default = 'ID,ko,pfam,ec_number,cog,locus_tag'

    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--input_gff', required=False, help='Path to input GFF3 file from IMG 5.0.')
    parser.add_argument('-p','--path_file', required=False, help='File contaning paths to input GFF3 file(s) from IMG 5.0 for processing multiple files.')
    parser.add_argument('-o', '--output', required=True, help='Output tsv file.')
    parser.add_argument('-f', '--fields', required=False, type=str, default=field_default, help='GFF attributes to parse. Input as a comma-separated string. Default: {}'.format(field_default))
    parser.add_argument('-m', '--img_map', required=False, help='Contig - IMG contig - GOLD OID - Sample - Bin map generated by img_bin_map.py.')
    parser.add_argument('-l', '--logfile', required=False, action='store', default=logfile_default, help='path to logfile')

    args = parser.parse_args()

    #Set up logger
    logging_format = '%(name)s :: %(levelname)s :: %(message)s :: %(asctime)s'

    logging.basicConfig(filename = args.logfile,
        level = logging.DEBUG,
        format = logging_format)

    print('Logfile written to: {}'.format(args.logfile))

    logging.info('SCRIPT: {}'.format(os.path.basename(__file__)))

    key_fields = list()
    #If --fields is selected but no options are specified, use the defaults
    if args.fields:
        key_fields = args.fields.split(',')
        if all([not k for k in key_fields]):
            key_fields = field_default
            print('Option --fields selected but empty. Using default values.')
        else:
            pass
    else:
        key_fields = field_default.split(',')

    if args.input_gff and not args.path_file:
        try:
            img_df = gff3_to_df(gff_file = args.input_gff,
            format = 'img5',
            attribute_select = key_fields)
        except IOError as e:
            print('could not parse GFF3 file {} - ensure this is not a path file'.format(args.input_gff))

    elif args.path_file and os.path.isfile(os.path.abspath(args.path_file)) and not args.input_gff:
        multi_gff_msg = 'Reading in multiple GFF3 files based on paths in {}'.format(args.path_file)
        print(multi_gff_msg)
        logging.info(multi_gff_msg)
        img_df = gff_pd_multiread(path_file = args.path_file, attr_sel = key_fields)

    elif args.path_file and args.input_gff:
        input_error_msg = 'Specify either a single GFF3 file, or a file containing multiple GFF3 paths for --input_gff'
        print('ERROR: {}'.format(input_error_msg))
        logging.error(input_error_msg)
        sys.exit(1)
    else:
        pass

    if args.img_map:
            try:
                print('reading GOLD --> Sample ID map to Pandas DF')
                img_map = pd.read_csv(args.img_map, sep = '\t')
            except IOError as e:
                img_map_error_msg = 'Could not open {}'.format(args.img_map)
                print('ERROR: {}'.format(img_map_error_msg))
                logging.error(img_map_error_msg)
            try:
                scaffold_anno_df = pd.merge(img_df, img_map,
                    on = 'IMG_Contig_Name', how = 'left')
            except:
                logging.error('could not merge img/gff3 dataframe with img_map')

            print(scaffold_anno_df.head())
    else:
        logging.info('No --img_map, ignoring merge...')

    #Column names for IMG Annotation output file for use with metabolic_profile.py
    accession_colname_dict = {'cog':'COG_ID', 'ko':'KO_Term',
    'pfam':'PFAM_ID', 'ec_number':'EC_Number', 'ID':'IMG_Gene_ID',
    'locus_tag':'Locus_Tag'}

    scaffold_anno_df.rename(columns=accession_colname_dict, inplace=True)

    #Write to the output TSV file
    try:
        scaffold_anno_df.to_csv(args.output, index=True, sep='\t', encoding='utf-8')
        print('Writing to output tsv {}'.format(args.output))
    except IOError as e:
        print('WARNING - unable to write output tsv file {}'.format(args.output))

#------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
