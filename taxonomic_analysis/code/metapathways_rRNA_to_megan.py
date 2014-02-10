#!/usr/bin/python
"""
metapathways_rRNA_to_megan.py

Created by Niels Hanson on 2013-08-16.
Copyright (c) 2013 Steven J. Hallam Laboratory. All rights reserved.
"""
from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright 2013"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Niels W Hanson"
__status__ = "Release"

try:
     import os
     import re
     import sys
     import argparse
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)


# Example: python metapathways_rRNA_to_megan.py -i output/HOT_Sanger/*/results/rRNA/*greengenes*txt -o .
what_i_do = """Parses processed rRNA result files <sample>.<database>.rRNA.stats.txt, e.g. my_sample.greengenes.rRNA.stats.txt, 
files from the <sample>/results/rRNA/ directory and formats them as csv files for import into MEGAN.
"""
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_files', type=str, nargs='+',
                required=True, help='glob of <sample>.<database>.rRNA.stats.txt files from the MetaPathways output folder <sample>/results/rRNA/ to be parsed', default=None)                
parser.add_argument('-o', dest='output_dir', type=str, nargs='?',
                required=False, help='directory where <sample>.<database>.megan.csv.txt files will be put', default=os.getcwd())
parser.add_argument('--dsv', dest='dsv', action='store_true',
                required=False, help='flag to output a .dsv instead of a .csv file', default=False)


def main(argv):
    args = vars(parser.parse_args())
    
    # setup input and output file_names
    input_files = args['input_files']
    output_dir = os.path.abspath(args['output_dir'])
    
    for f in input_files:
        file_handle = open(f, "r")
        lines = file_handle.readlines()
        file_handle.close()
        
        end = ".csv.txt"
        if args['dsv']:
            end = ".dsv"
        # clean up filename
        sample_db = re.sub(".rRNA.stats.txt", "", os.path.basename(f), re.I)
        output_file = [output_dir, os.sep, sample_db, ".megan", end]
        output_handle = open("".join(output_file), "w")

        for l in lines:
           fields = l.split("\t")
           if len(fields) == 7 and len(fields[0].strip()) > 2:
               # probably a result line
               read = fields[0].strip()
               last_score = fields[5].strip()
               taxa = fields[6].strip("\n").strip()
               taxa = re.sub(r"[a-z]__", "", taxa)
               taxa = re.sub(r";+", ";", taxa)
               out_line = read + ", " + taxa + ", " + last_score + "\n"
               output_handle.write(out_line)
        
        output_handle.close()
    
    exit()

# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
