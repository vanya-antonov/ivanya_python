
import types
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from generate_random_dna import make_random_dna_seq, make_random_dna_records


class Test(unittest.TestCase):
    def test_make_random_dna_records(self):
        generator = make_random_dna_records(
            dna_len=100, dna_gc=30, num_records=10, prefix='test_')
        all_recs = list(generator)

        # Check the types
        self.assertEqual(type(generator), types.GeneratorType)
        self.assertEqual(type(all_recs), list)
        self.assertEqual(type(all_recs[0]), SeqRecord)
        with self.assertRaises(TypeError):
            # TypeError: 'generator' object is not subscriptable
            x = generator[0]

        # Check length and GC content
        self.assertEqual(len(all_recs), 10)
        self.assertTrue(all([len(r.seq) == 100 for r in all_recs]))
        self.assertTrue(all([GC(r.seq) == 30 for r in all_recs]))

        # Random sequences should be different
        self.assertNotEqual(all_recs[0].seq, all_recs[1].seq)

        # Check record attributes
        self.assertEqual(all_recs[0].id, 'test_1')
        self.assertEqual(all_recs[9].id, 'test_10')
        self.assertEqual(all_recs[0].name, '')
        self.assertEqual(all_recs[0].description, '')
        
    def test_make_random_dna_records_impossible_gc(self):
        """Tests the function behaviour when the sequence with the exact GC
        content can't be generated.
        """
        generator51 = make_random_dna_records(
            dna_len=10, dna_gc=51, num_records=10)
        generator59 = make_random_dna_records(
            dna_len=10, dna_gc=59, num_records=10)

        # The first generator call should produce a warning
        with self.assertLogs():
            all_recs_51 = list(generator51)
        with self.assertLogs():
            all_recs_59 = list(generator59)

        # Seqs should have GC=50% instead of 51% and 60% instead of 59%
        self.assertFalse(all([GC(r.seq) == 51 for r in all_recs_51]))
        self.assertTrue(all([GC(r.seq) == 50 for r in all_recs_51]))
        self.assertFalse(all([GC(r.seq) == 59 for r in all_recs_59]))
        self.assertTrue(all([GC(r.seq) == 60 for r in all_recs_59]))

    def test_make_random_dna_seq(self):
        seq = make_random_dna_seq(dna_len=10, num_gc=5)

        # Check the type, length and GC content
        self.assertEqual(type(seq), Seq)
        self.assertEqual(len(seq), 10)
        self.assertEqual(GC(seq), 50.0)

        # num_gc must be integer
        with self.assertRaises(TypeError):
            make_random_dna_seq(dna_len=10, num_gc=5.1)

