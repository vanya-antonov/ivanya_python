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


def main(args):
    check_input(args)
    rand_seqs = make_random_dna_records(
        args.dna_len, args.gc, args.num, prefix=args.prefix)
    SeqIO.write(rand_seqs, sys.stdout, args.format)

def check_input(args):
    "Sets the logging level and validates the user provided parameters."
    if args.quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)

    if args.dna_len < 1:
        raise ValueError("Wrong DNA length = %s" % args.dna_len)

    if args.gc < 0 or args.gc > 100:
        raise ValueError("Wrong DNA GC content = %.1f%%" % args.gc)

    if 0 < args.gc < 1:
        logging.warning("Didn't you mean GC content = %.0f%% (not %.2f%%)?" %
                        (100*args.gc, args.gc))

def make_random_dna_records(dna_len, dna_gc, num_records, prefix=''):
    "Returs a generator of SeqRecord objects containing random DNA sequences."

    # Estimate the number of AT and GC letters in the seq
    num_gc = int(round(dna_len*dna_gc/100, 0))
    if 100*num_gc/dna_len != dna_gc:
        logging.warning("The GC content of the generated sequences is %.2f%% "
                       "(instead of %.2f%%)" % (100*num_gc/dna_len, dna_gc))

    for i in range(num_records):
        yield SeqRecord(
            seq=make_random_dna_seq(dna_len, num_gc),
            id=prefix + str(i+1),
            name="",
            description="")

def make_random_dna_seq(dna_len, num_gc):
    """Returns a Bio.Seq object containing random DNA sequence with a
    specific length and number of G/C letters.
    """
    nt_list = [random.choice('GC') for _ in range(num_gc)]

    num_at = dna_len - num_gc
    nt_list += [random.choice('AT') for _ in range(num_at)]

    # Shuffle list in place: https://stackoverflow.com/a/2668325/310453
    random.shuffle(nt_list)

    return Seq(''.join(nt_list), generic_dna)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate random DNA sequences with fixed length "
        "and GC content.")
    parser.add_argument('dna_len', metavar='DNA_LEN', type=int,
                        help='DNA sequence length')
    parser.add_argument('-g', '--gc', metavar='GC', type=float, default=50.0,
                        help='DNA sequence GC content (default 50)')
    parser.add_argument('-n', '--num', metavar='NUM', type=int, default=1,
                        help="number of random sequences to generate (default 1)")
    parser.add_argument('-p', '--prefix', metavar='STR', default='',
                        help="prefix for the sequence names (e.g. 'rand_')")
    parser.add_argument('-f', '--format', metavar='STR', default='fasta-2line',
                        help="output format, default is 'fasta-2line'. "
                        "See the list of other formats at "
                        "https://biopython.org/wiki/SeqIO")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="do not write warning messages to stderr")
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())

