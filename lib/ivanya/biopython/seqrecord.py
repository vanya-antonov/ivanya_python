# Copyright 2018 by Ivan Antonov. All rights reserved.


from Bio.SeqRecord import SeqRecord

from ivanya.biopython.seq import make_random_dna_seq


def make_random_dna_records(dna_len, num_records, gc=50.0, prefix='rand_'):
    """Returs a list of SeqRecord objects containing random DNA sequences."""
    all_records = []
    for i in range(num_records):
        record = SeqRecord(
            seq=make_random_dna_seq(dna_len, gc),
            id=prefix + str(i+1),
            name="",
            description="")
        all_records.append(record)
    return all_records

