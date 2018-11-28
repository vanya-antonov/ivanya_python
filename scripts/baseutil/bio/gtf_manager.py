#!/usr/bin/env python3

# $Id: gtf_manager.py 2905 2018-08-07 15:42:08Z antonov $

###
# Ivan Antonov (vanya.antonov@gmail.com)
#

import sys, os, re
from pprint import pprint
from subprocess import Popen, PIPE

###
# CONSTANTS

###
# SUBROUTINES
def run(opts):
	if opts['todo'] == 'make_long_chimerome':
		make_long_chimerome(opts['arg1'], opts)
	else:
		usage("Unknown TODO = '%s'" % opts['todo'])


def make_long_chimerome(gtf_fn, opts={}):
	keys = ['chr', 'c2', 'c3', 'start', 'end', 'c6',  'strand',  'c8',  'info']
	(cur_exons, processed_genes) = ([], {})
	with open(gtf_fn) as f:
		for line in f:
			if line == '\n':
				continue
			line.rstrip()   # remove \n
			
			vals = line.split('\t')
			if len(vals) != len(keys):
				sys.stderr.write("Wrong line: '%s'", line)
				continue
			
			# chr1  FANTOM6 exon    91421   91629   .   -   .   gene_id "ENSG00000225880"; transcript_id "FTMT20100027365.C1";
			exon = dict(zip(keys, vals))
			gene_mo = re.compile(r'gene_id\s*"(.+?)"').search(exon['info'])
			if gene_mo == None:
				sys.stderr.write("Can't find gene_id in string '%s'", exon['info'])
				continue
			exon['gene_id'] = gene_mo.group(1)
			
			if (len(cur_exons) == 0) or (cur_exons[0]['gene_id'] == exon['gene_id']):
				cur_exons.append( exon )
			else:
				# Exons of the next gene have begun
				sys.stderr.write("\rGenerating chimera for '%s' with %d exons...      " % (cur_exons[0]['gene_id'], len(cur_exons)) )
				_print_long_chimeric_trx(cur_exons)
				processed_genes[ cur_exons[0]['gene_id'] ] = 1
				
				assert exon['gene_id'] not in processed_genes.keys(), "The input file is not sorted: gene '%s' is among already processed genes!" % exon['gene_id']
				cur_exons = [ exon ]
		# END: for line in f
		_print_long_chimeric_trx(cur_exons)
	# END: with open(gtf_fn) as f
	sys.exit()

def _print_long_chimeric_trx(all_exons):
	if( len(all_exons) == 0 ):
		return
	
	bed_txt = ''
	for exon in all_exons:
		# Print something like 'chr1	1018110	1018979	CATG00000000002	.	+'
		vals = [exon[k] for k in ['chr',  'start',  'end',  'gene_id',  'c8',  'strand']]
		bed_txt += "\t".join(vals) + "\n"
	
	# https://stackoverflow.com/a/8475367/310453
	proc = Popen(['bedtools', 'merge', '-s'], stdin=PIPE, stdout=PIPE)
	out_bytes = proc.communicate( bytes(bed_txt, 'utf-8') )[0]
	
	# bed2gff
	gene_id  = all_exons[0]['gene_id']
	info_str = 'gene_id "%s"; transcript_id "%s"' % (gene_id, gene_id)
	for line in out_bytes.decode("utf-8").split("\n"):
		vals = line.split()
		if len(vals) == 0:
			continue
		elif len(vals) != 4:
			sys.stderr.write("Wrong line in the bedtools output: '%s'" % line)
			continue
		# line = 'chr1	4873173	4873320	+'
		(chrom, left, right, strand) = vals
		print("\t".join([chrom, 'chimerome', 'exon', left, right, '.', strand, '.', info_str]))

def usage(msg = ''):
	script = os.path.basename(sys.argv[0])
	sys.stderr.write('''%(msg)s
DESCRIPTION:
    <TODO>
    * make_long_chimerome   <ANNOTATION.gtf>   >   <CHIMEROME.gtf>
        - Requirements: bedtools (v2.26.0) and gffread (0.9.9)
        - The input file must be sorted by chrom, the gene_id, then start:
          sort -k1,1 -k10,10 -k4,4n F6_CAT.transcript.gtf > F6_CAT.transcript.SORTED.gtf

USAGE:
    %(script)s   [OPTIONS]   <TODO>   <ARG1>   <ARG2> ...
    
OPTIONS:
    --silent\n''' % locals() )

###
# Parse command line arguments
if len(sys.argv) < 2:
	usage()
	sys.exit()

###
#my $START_TIME = time;
run({
    'todo' : sys.argv[1],
	'arg1' : sys.argv[2] if len(sys.argv) >=2 else '',
});
#warn "\nElapsed time: ".(time-$START_TIME)." sec\n" if !$SILENT;
###

