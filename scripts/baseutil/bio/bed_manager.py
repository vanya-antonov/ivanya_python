#!/usr/bin/env python3

# $Id: bed_manager.py 2945 2018-10-24 10:49:05Z antonov $

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
    if args.todo == 'bkg4bed':
        bkg4bed(args.arg1, args.arg2, args.bed)

def bkg4bed(bed_fn, chrom_sizes_fn, exclude_bed=None):
    len_dict = _count_lens(bed_fn=bed_fn)
    while len_dict:
        logging.info('Generating %d regions for %d different lengths...' %
                     (sum(len_dict.values()), len(len_dict.keys())))
        bkg_bed_txt = _sample_regions_for_len_dict(len_dict, chrom_sizes_fn)
        rejected_txt = ''
        if exclude_bed is not None:
            bkg_bed_txt, rejected_txt = _process_exclude_bed(bkg_bed_txt, exclude_bed)
        print(bkg_bed_txt, end='')
        len_dict = _count_lens(bed_txt=rejected_txt)

def _process_exclude_bed(input_bed_txt, exclude_bed):
    good_txt = subprocess.run(
        'bedtools intersect  -v  -a -  -b ' + exclude_bed,
        shell=True, stdout=subprocess.PIPE, input=input_bed_txt, universal_newlines=True
    ).stdout
    
    bad_txt = subprocess.run(
        'bedtools intersect  -wa  -a -  -b ' + exclude_bed,
        shell=True, stdout=subprocess.PIPE, input=input_bed_txt, universal_newlines=True
    ).stdout
    
    return good_txt, bad_txt

def _sample_regions_for_len_dict(len_dict, chrom_sizes_fn):
    if chrom_sizes_fn is None:
        raise Exception('chrom_sizes_fn if required!')
    
    bed_list = []
    for (reg_len, reg_num) in len_dict.items():
        bed_list.append(subprocess.getoutput(
            'bedtools random -l ' + str(reg_len) +
            ' -n ' + str(reg_num) + ' -g ' + chrom_sizes_fn))
    
    return "\n".join(bed_list)

def _count_lens(bed_fn=None, bed_txt=None):
    if bed_fn is not None:
        iterator = fileinput.input(files=bed_fn)
    elif bed_txt is not None:
        iterator = bed_txt.splitlines()
    else:
        Exception('Either bed_fn or bed_txt is required!')
    #for line in fileinput.input(files=bed_fn):
    len_dict = {}
    for line in iterator:
        if re.compile(r'^\s*$').match(line):
            continue
        vals = line.split()
        num = int(vals[2]) - int(vals[1])
        if num not in len_dict:
            len_dict[num] = 1
        else:
            len_dict[num] += 1
    
    return len_dict

def parse_args():
    # Parse command line arguments: https://stackoverflow.com/a/30493366/310453
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
requirements:
    bedtools (v2.27.1)

description:
    <TODO>
     - bkg4bed  <FN.bed>  <FN.chrom.sizes>  --  randomly select the same number of regions with the same length
                as regions in the input .bed file. <FN.chrom.sizes> can be obtain by:
                fetchChromSizes hg38 > hg38.chrom.sizes
                --bed  <STOP.bed>  --  the generated BKG-regions and STOP-regions must not overlap''')
    all_todo = ['bkg4bed']
    parser.add_argument('todo', metavar='TODO', choices=all_todo, help=', '.join(all_todo))
    parser.add_argument('arg1', nargs='?', help='optional argument (see specific TODO for details)')
    parser.add_argument('arg2', nargs='?', help='optional argument (see specific TODO for details)')
    #parser.add_argument('--file', metavar='FN', help='additional file (see specific TODO for details)')
    parser.add_argument('--bed', metavar='FN.bed', help='additional bed-file (see specific TODO for details)')
    
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())

