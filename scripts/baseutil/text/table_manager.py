#!/usr/bin/env python3

# $Id: table_manager.py 2913 2018-08-16 15:38:12Z antonov $

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
    if args.todo == 'cut':
        better_cut(args.input, args.cols)
    elif args.todo == 'select_inline':
        select_inline(args.input, args)

def better_cut(input_fn, cols):
    if len(cols) == 0:
        Exception('--cols option is required!')
    
    col_i = []
    for num in cols:
        col_i += [num-1] if num > 0 else [num]
    
    for line in fileinput.input(files=input_fn):
        fields = line.split()
        subset = [fields[i] for i in col_i]
        print("\t".join(subset))

def select_inline(input_fn, args, inline_delim=','):
    if len(args.cols) != 1:
        Exception('--cols option is required!')
    tar_col_i = args.cols[0] - 1
    
    max_col_num = max(args.cols + args.other_cols)
    for line in fileinput.input(files=input_fn):
        fields = line.split()
        if len(fields) < max_col_num:
            Exception('Input file has wrong format (%d < %d)!' %(len(fields), max_col_num))
        
        tar_vals = [float(x) for x in fields[tar_col_i].split(inline_delim)]
        if args.ops == 'min':
            tar_i = tar_vals.index(min(tar_vals))
        else:
            Exception('Unknown operation = "%s"!' % args.ops)
        fields[tar_col_i] = str(tar_vals[tar_i])
        
        for col_i in [num - 1 for num in args.other_cols]:
            vals = fields[col_i].split(inline_delim)
            if len(vals) != len(tar_vals):
                Exception('Column %d has wrong number of values: "%s"!' % (col_i+1, fields[col_i]))
            fields[col_i] = vals[tar_i]
        
        print("\t".join(fields))

def parse_args():
    # Parse command line arguments: https://stackoverflow.com/a/30493366/310453
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
description:
    <TODO>
        cut  --  an improved version of the UNIX cut. -c option is required!
        select_inline  -- similat to the 'groupby' or 'merge' bedtools TODOs:
103858730	7213	8736	YP_353360.1,ABL66375.1,ABU57630.1,BAA10114.1	3.72e-89,6.45e-101,5.01e-117,0.0	+,+,+,+
-----------------------------------------------------------
|   select_inline  --op min  --cols 5  --other_cols 4 6   |
-----------------------------------------------------------
103858730	7213	8736	BAA10114.1	0.0	+
        ''')
    all_todo = ['cut', 'select_inline']
    parser.add_argument('todo', metavar='TODO', choices=all_todo, help=', '.join(all_todo))
    parser.add_argument('input', nargs='?', default='-', help='table file name or "-" to use stdin')
    #parser.add_argument('arg2', nargs='?', help='optional argument (see specific TODO for details)')
    parser.add_argument('-c', '--cols', metavar='LIST', nargs='+', type=int, default=[],
                        help='Specify the target columns (1-based)')
    parser.add_argument('-o', '--ops', metavar='STR',
                        help="Specify the operation that should be applied to target columns.\nValid operations: min")
    parser.add_argument('--other_cols', metavar='LIST', nargs='+', type=int, default=[],
                        help='A list of other column indexes (1-based)')
    
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())

