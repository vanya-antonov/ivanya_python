#!/usr/bin/env python3

# $Id: gbk_manager.py 2905 2018-08-07 15:42:08Z antonov $

###
# Ivan Antonov (vanya.antonov@gmail.com)
#

import sys, os, re, argparse
from pprint import pprint
from warnings import warn

from Bio import SeqIO

###
# SUBROUTINES
def run(args):
	if(args.todo == 'summary') : print_summary(args)


#---
def print_summary(args):
	head = ['file','id','id_ver','len','type','BioProject','BioSample','Assembly','Assembly_ver','all_ids','org','kingdom','phylum','genera','taxonomy']
	if not args.no_header: print("\t".join(head))
	for seq_record in SeqIO.parse(args.gbk_fn, "genbank"):
		dbx_ll   = [str.split(':') for str in seq_record.dbxrefs]   # list of lists
		dbx_dict = {l[0] : l[1] for l in dbx_ll}
		asmbl    = dbx_dict.get('Assembly',   '-')
		taxa_l   = seq_record.annotations.get('taxonomy',[])
		info = {
			'file'         : os.path.basename(args.gbk_fn),
			'id'           : re.compile(r'\.\d+$').sub('', seq_record.id),   # NC_000916.1  =>  NC_000916
			'id_ver'       : seq_record.id,
			'org'          : seq_record.annotations['organism'],
			'kingdom'      : taxa_l[0]  if len(taxa_l) > 0 else '-',
			'phylum'       : taxa_l[1]  if len(taxa_l) > 2 else '-',
			'genera'       : taxa_l[-1] if len(taxa_l) > 2 else '-',
			'len'          : len(seq_record),
			'type'         : seq_record.annotations['molecule_type'],
			'BioProject'   : dbx_dict.get('BioProject', '-'),
			'BioSample'    : dbx_dict.get('BioSample',  '-'),
			'Assembly'     : re.compile(r'\.\d+$').sub('', asmbl),           # GCF_000007565.2  =>  GCF_000007565
			'Assembly_ver' : asmbl,
			'all_ids'      : ';'.join(seq_record.annotations['accessions']),
			'taxonomy'     : ';'.join(taxa_l),
		}
		print("\t".join([str(info.get(k,'-')) for k in head]))


###
# Parse command line arguments: https://stackoverflow.com/a/30493366/310453
parser = argparse.ArgumentParser(description='Parsing GenBank files with Biopython')
parser.add_argument('todo',              metavar='TODO',       help="supported values: 'summary'", choices=['summary'])
parser.add_argument('gbk_fn',            metavar='SEQ.gbk',    help='input file in GenBank format')
parser.add_argument('-s', '--silent',    action ='store_true', help='do not print progress messages')
parser.add_argument('-H', '--no_header', action ='store_true', help='do not print table header')
args = parser.parse_args()

args.gbk_fn = os.path.abspath(args.gbk_fn)

###
run( args )
###

