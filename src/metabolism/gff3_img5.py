"""
Author: Ian Rambo - ian.rambo@utexas.edu
Created: August 23, 2021
Last updated: August 25, 2021

Purpose: Parse IMG 5.0 GFF3 output for metabolism and virus-host analyses.
"""

import argparse
import logging
import re
#------------------------------------------------------------------------------
def gff3_to_dict(gff_file, format='standard'):
    """
    Parse a GFF3 file and return a dictionary.
    """

    gff_dict = {}
    with open(gff_file, 'r', encoding = 'utf-8') as gff:
        for record in gff:
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

                if format == 'standard':
                    gff_dict[seqid] = {'source':source, 'ftype':ftype,
                    'start':start, 'end':end, 'score':score, 'strand':strand,
                    'phase':phase, 'attributes':attributes_dict}
                elif format == 'img':
                    

                else:
                    pass

    return gff_dict
#------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, help="Input GFF3 file from IMG 5.0")
    parser.add_argument('-o', '--output', required=True, help="Output tsv file.")
    args = parser.parse_args()

if __name__ == "__main__":
    main()
