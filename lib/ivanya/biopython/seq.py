# Copyright 2018 by Ivan Antonov. All rights reserved.


import random

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def _make_random_dna(dna_len, dna_gc):
    # https://stackoverflow.com/a/21205929/310453
    rand_nt = []
    for _ in range(dna_len):
        if random.uniform(0, 100) < dna_gc:
            letters = 'CG'
        else:
            letters = 'AT'
        rand_nt.append(random.choice(letters))
    return ''.join(rand_nt)

def make_random_dna_seq(dna_len, dna_gc, max_diff=1):
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

    num_tries = 0
    while True:
        rand_dna = _make_random_dna(dna_len, dna_gc)
        if abs(GC(rand_dna) - dna_gc) < max_diff:
            # The difference between the GC-content of the generated seq and
            # the desired value is acceptable
            return Seq(rand_dna, generic_dna)
        elif num_tries > 1000:
            raise Exception("Can't generate random DNA with length = '%s' and "
                           "GC content = '%.1f'" % (dna_len, dna_gc))
        num_tries += 1

