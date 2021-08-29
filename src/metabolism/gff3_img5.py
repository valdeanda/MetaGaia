"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: August 23, 2021
Last updated: August 25, 2021

Purpose: Parse IMG 5.0 GFF3 output for MetaGaia metabolic analysis.
"""

import argparse
import logging
import re
import pandas as pd
#------------------------------------------------------------------------------
def gff3_to_df(gff_file, format='img5', key_select = ['ID', 'ko', 'pfam', 'cog', 'ec_number']):
    """
    Parse a GFF3 file and return a Pandas dataframe.
    """

    valid_formats = ['standard', 'img5']
    gff_data = []
    try:
        gff_handle = open(gff_file, 'r', encoding = 'utf-8')
    except IOError as e:
        print('Could not open {}'.format(gff_file))

    for record in gff_handle:
        record_dict = {}
        if not record.startswith('#'):
            record = record.rstrip()
            recordList = record.split()
            seqid = recordList[0]
            source = recordList[1]
            ftype = recordList[2]
            start = recordList[3]
            end = recordList[4]
            score = recordList[5]
            strand = recordList[6]
            phase = recordList[7]
            attributes = recordList[8]
            attributes_list = attributes.split(';')

            attributes_dict = {}
            for attribute in attribute_list:
                alist = attribute.split('=')
                att_key = alist[0]
                att_value = alist[1]
                attributes_dict[att_key] = att_value

            if not format in valid_formats:
                print('ERROR - requires a valid format, choose from: {}'.format(' '.join(valid_formats)))
                exit(1)

            elif format == 'standard':
                gff_dict = {'Contig_Name':seqid, 'source':source, 'ftype':ftype,
                'start':start, 'end':end, 'score':score, 'strand':strand,
                'phase':phase, 'attributes':attributes_dict}
            elif format == 'img5':
                record_dict['IMG_Contig_Name'] = seqid
                #GOLD annotation keys used in metabolic analyses
                for key in key_select:
                    record_dict[key] = "NaN"

                for key in attributes_dict.keys():
                    if key in key_select:
                        record_dict[key] = attributes_dict[key]
                    else:
                        pass

            else:
                pass

            gff_data.append(record_dict)

        else:
            pass

    gff_handle.close()

    gff_df = pd.DataFrame.from_dict(gff_data)

    return gff_df
#------------------------------------------------------------------------------
def map_gold(scaffold_map, sample_map, scaffold_samples=False):
    """
    Create a unified mapping file for original IMG and Gold scaffolds.
    """
    try:
        print('reading IMG - GOLD scaffold map to Pandas DF')
        gold_scaffold_map = pd.read_csv(scaffold_map, sep = '\t',
            names = ['Original_Contig_Name', 'IMG_Contig_Name'])
    except IOError as e:
        print('ERROR: could not open {}'.format(scaffold_map))

    try:
        print('reading GOLD - Sample ID map to Pandas DF')
        gold_sample_map = pd.read_csv(sample_map, sep = '\t')
    except IOError as e:
        print('ERROR: could not open {}'.format(sample_map))

    gold_scaffold_map['GOLD_OID'] = gold_scaffold_map['IMG_Contig_Name'].str.split('_').str[0].str.strip()

    gold_scaffold_sample_df = pd.merge(scaffold_map, sample_map, on = 'GOLD_OID', how = 'right')

    if scaffold_samples:
        gold_scaffold_sample_df['Original_Contig_Name'] = gold_scaffold_sample_df['SAMPLE_ID'].astype(str) + '_' + gold_scaffold_sample_df['Original_Contig_Name'].astype(str)
    else:
        pass

    return gold_scaffold_sample_df
#------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_gff', required=True, help="Input GFF3 file from IMG 5.0")
    parser.add_argument('-o', '--output', required=True, help="Output tsv file.")
    parser.add_argument('-f', '--fields', required=False, type='str', help="GFF attributes to parse. Input as a comma-separated string: KEGG,EC_NUMBER,PFAM")
    parser.add_argument('-m', '--gold_scaffold_map', required=True, help="Scaffold mapping file from IMG GOLD, containing original IMG scaffold ID and GOLD scaffold ID")
    parser.add_argument('-g', '--gold_sample_map', required=True, help="Mapping file of GOLD IDs to Sample IDs")
    parser.add_argument('-s', '--samples_to_id', required=False, action='store_true', help='adds sample IDs to the beginning of the IMG original scaffold IDs')
    args = parser.parse_args()

if args.fields:
    key_fields = args.fields.upper().split(',')
    if key_fields:
        print('extracting selected fields from GFF: {}'.format(args.fields))
        img_df = gff3_to_df(gff_file = args.input_gff, format = 'img5', key_select = args.fields)
    else:
        print('ERROR: could not parse desired fields from --fields argument, ensure input is comma delimited, e.g. pfam,ko,ec_number')
        exit(1)
else:
    img_df = gff3_to_df(gff_file = args.input_gff, format = 'img5')

if args.samples_to_id:
    gold_scaffold_sample_df = map_gold(scaffold_map = args.gold_scaffold_map,
    sample_map = arg.gold_sample_map, scaffold_samples=True)

else:
    gold_scaffold_sample_df = map_gold(scaffold_map = args.gold_scaffold_map,
    sample_map = arg.gold_sample_map)

scaffold_anno_df = pd.merge(img_df, gold_scaffold_sample_df,
    on='IMG_Contig_Name', how='left')

print(scaffold_anno_df.head())

#scaffold_anno_df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()
