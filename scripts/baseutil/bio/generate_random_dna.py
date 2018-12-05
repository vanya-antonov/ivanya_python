#!/usr/bin/env python3

# Copyright 2018 by Ivan Antonov. All rights reserved.

"""Generate random DNA sequences."""

import argparse
import logging
import random
import sys

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC


def main(args):
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    else:
        logging.getLogger().setLevel(logging.INFO)

    rand_seqs = make_random_dna_records(
        args.dna_len, args.gc, args.num, prefix=args.prefix,
        max_diff=args.max_diff)
    SeqIO.write(rand_seqs, sys.stdout, args.format)

def make_random_dna_records(dna_len, dna_gc, num_records, prefix='', max_diff=1):
    """Returs a list of SeqRecord objects containing random DNA sequences."""
    all_records = []
    for i in range(num_records):
        record = SeqRecord(
            seq=make_random_dna_seq(dna_len, dna_gc, max_diff),
            id=prefix + str(i+1),
            name="",
            description="")
        all_records.append(record)
    return all_records

def make_random_dna_seq(dna_len, dna_gc, max_diff):
    """Generates random DNA sequence with a specific length and GC content.

    Arguments:
     - dna_len - integer
     - dna_gc - floating number between 0 and 100
     - max_diff - maximum allowed difference between the desired GC content
       and the actual GC content of the generated random sequence (in %)
    """
    if dna_gc < 1:
        raise Exception("Use percents instead of fractions for GC content "
                       "(the provided dna_gc = '%.2f')" % dna_gc)

    num_tries = 1
    while True:
        rand_dna = make_random_dna(dna_len, dna_gc)
        rand_gc = GC(rand_dna)
        if abs(rand_gc - dna_gc) < max_diff:
            # The difference between the GC-content of the generated seq and
            # the desired value is acceptable
            return Seq(rand_dna, generic_dna)
        elif num_tries > 1000:
            raise Exception("Can't generate random DNA with length = '%s' and "
                           "GC content = '%.1f'" % (dna_len, dna_gc))
        logging.info("[%s] Random DNA is discarded because its GC=%.2f%%" %
                     (num_tries, rand_gc))
        num_tries += 1

def make_random_dna(dna_len, dna_gc):
    # https://stackoverflow.com/a/21205929/310453
    rand_nt = []
    for _ in range(dna_len):
        if random.uniform(0, 100) < dna_gc:
            letters = 'CG'
        else:
            letters = 'AT'
        rand_nt.append(random.choice(letters))
    return ''.join(rand_nt)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate random DNA sequences with fixed length "
        "and GC content.")
    parser.add_argument('dna_len', metavar='DNA_LEN', type=int,
                        help='DNA sequence length')
    parser.add_argument('gc', metavar='DNA_GC', type=float,
                        help='DNA sequence GC content (between 0 and 100)')
    parser.add_argument('num', metavar='NUM_SEQS', type=int,
                        help="number of random sequences to generate")
    parser.add_argument('--prefix', metavar='str', default='',
                        help="prefix for the sequence names (e.g. 'rand_')")
    parser.add_argument('--max_diff', metavar='GC', type=float, default=5.0,
                        help="the maximum acceptable difference between the "
                       "DNA_GC and the actual GC content of random sequence; "
                       "default is 5.0%")
    parser.add_argument('--format', metavar='str', default='fasta-2line',
                        help="output format, default is 'fasta-2line'. "
                        "See the list of other formats at "
                        "https://biopython.org/wiki/SeqIO")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="do not write info messages to stderr")
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())

