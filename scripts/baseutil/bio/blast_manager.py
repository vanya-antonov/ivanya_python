#!/usr/bin/env python3

# $Id: blast_manager.py 2908 2018-08-10 13:29:30Z antonov $

###
# Ivan Antonov (vanya.antonov@gmail.com)
#

import argparse
import fileinput
import os
import re
import subprocess
import sys
from pprint import pprint, pformat

import logging    # https://www.youtube.com/watch?v=-RcDmGNSuvU
logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

###
# SUBROUTINES
def main(args):
    if args.todo == 'tblastn2bed':
        tblastn2bed(args.arg1, args)

def tblastn2bed(input_fn, args):
    fields = ['query_id', 'hit_id', 'identity', 'ali_length', 'mismatches', 'gap_opens',
              'q_start', 'q_end', 'h_start', 'h_end', 'evalue', 'score']
    for line in fileinput.input(files=input_fn):
        vals = line.split()
        if len(fields) != len(vals):
            Exception('Input file has wrong format!')
        dd = dict(zip(fields, vals))
        
        dd['strand'] = '+'
        if int(dd['h_end']) < int(dd['h_start']):
            dd['h_start'], dd['h_end'], dd['strand'] = dd['h_end'], dd['h_start'], '-'
        
        bed_keys = ['hit_id', 'h_start', 'h_end', 'query_id', 'evalue', 'strand']
        print("\t".join(str(dd[k]) for k in bed_keys))

def parse_args():
    # Parse command line arguments: https://stackoverflow.com/a/30493366/310453
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
requirements:
    bedtools (v2.27.1)

description:
    <TODO>
        tblastn2bed  <BLAST_OUT>  -- ''')
    all_todo = ['tblastn2bed']
    parser.add_argument('todo', metavar='TODO', choices=all_todo, help=', '.join(all_todo))
    parser.add_argument('arg1', nargs='?', help='optional argument (see specific TODO for details)')
    parser.add_argument('arg2', nargs='?', help='optional argument (see specific TODO for details)')
    #parser.add_argument('--file', metavar='FN', help='additional file (see specific TODO for details)')
    #parser.add_argument('--bed', metavar='FN.bed', help='additional bed-file (see specific TODO for details)')
    
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())

