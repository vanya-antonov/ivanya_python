
# $Id$

import os
import re
import subprocess
import sys
import logging    # https://www.youtube.com/watch?v=-RcDmGNSuvU
from pprint import pprint

from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Blast import NCBIXML
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def is_seq_dual_coding(seq, frame2, gencode):
    """The function checks if given seq is dual coding, i.e. does not have
    inner stop codons in the main as well as in the alternative (frame2)
    reading frame. Stop codons at sequence ends are allowed.
    """
    if frame2 == +1:
        seq2 = seq[2:]   # skip the first 2 nt
    elif frame2 == -1:
        seq2 = seq[1:]   # skip the first nt
    else:
        raise Exception("Unknown frame = '%d'" % frame2)
    
    # Make the whole number of codons to avoid Biopythong warnings
    seq1 = seq[0:(3*int(len(seq)/3))]
    seq2 = seq2[0:(3*int(len(seq2)/3))]
    
    prot1 = seq1.translate(gencode).strip('*')
    prot2 = seq2.translate(gencode).strip('*')
    if '*' in prot1 + prot2:
        return False
    else:
        return True

def get_dual_coding_stop_stop_location(
    fs_coord, fs_type, strand, chr_seq, gencode):
    """Extend the starting coordinate (fs_coord) downstream until the stop-codon 
    in the main frame, shift the frame (fs_type) at this stop-codon and then
    go upstream in the alternative until the second stop codon.
    
    For GeneTack prediction this algorithm should result in a sequence that does NOT
    have stop-codons in both the main and the alternative (fs_type) frames.
    
    Returns a Bio.SeqFeature.FeatureLocation object.
    
    >>> from Bio.Alphabet import IUPAC
    >>> from Bio.Seq import Seq
    >>> from mylib.baseutilbio import get_dual_coding_stop_stop_location
    >>>
    >>> seq = Seq("ATGATAATTGAAAAATGA", IUPAC.unambiguous_dna)
    >>> loc_minus = get_dual_coding_stop_stop_location(9, -1, 1, seq, 11)
    >>> loc_minus.extract(seq)   # returns "TGAAAAATGA"
    >>> loc_plus = get_dual_coding_stop_stop_location(9, +1, 1, seq, 11)
    >>> loc_plus.extract(seq)   # returns "TAATTGAAAAATGA"
    >>>
    >>> seq_rc = seq.reverse_complement()
    >>> loc_rc = get_dual_coding_stop_stop_location(9, -1, -1, seq_rc, 11)
    >>> loc_rc.extract(seq_rc)   # returns "TGAAAAATGA"
    """
    if strand == -1:
        # Yes, this is not the most efficient function:)
        chr_seq = chr_seq.reverse_complement()
        fs_coord = len(chr_seq) - fs_coord
    
    ss_right = _walk_to_the_right_till_stop(fs_coord, chr_seq, gencode)
    if fs_type == -1:
        shifted_right = ss_right - 1
    elif fs_type == +1:
        shifted_right = ss_right - 2
    else:
        raise Exception("Unknown fs_type = '%d'" % fs_type)
    
    ss_left = _walk_to_the_left_till_stop(shifted_right, chr_seq, gencode)
    if not ss_left <= fs_coord <= ss_right:
        logging.error("FS_coord '%d' is outside the stop-stop seq '%d-%d'" %
                     (fs_coord, ss_left, ss_right))
        return None
    
    ss_seq = chr_seq[ss_left:ss_right]
    if not is_seq_dual_coding(ss_seq, fs_type, gencode):
        logging.error("The sequence is NOT dual coding:\n%s" % ss_seq)
        return None
    
    if strand == -1:
        new_left = len(chr_seq) - ss_right
        new_right = len(chr_seq) - ss_left
        ss_left, ss_right = new_left, new_right
    return FeatureLocation(ss_left, ss_right, strand)

def _walk_to_the_right_till_stop(start_coord, chr_seq, gencode):
    """Positive strand is assumed."""
    ss_right = start_coord
    while ss_right <= len(chr_seq)-3:
        cur_codon = chr_seq[ ss_right : (ss_right+3) ]
        codon_type = get_codon_type(cur_codon, gencode)
        if codon_type == 'coding':
            ss_right += 3  # expand the stop-stop seq
        elif codon_type == 'stop':
            ss_right += 3  # stop-stop seq includes stops as well
            break
        else:
            break
    return ss_right
    
def _walk_to_the_left_till_stop(start_coord, chr_seq, gencode):
    """Positive strand is assumed."""
    ss_left = start_coord
    while ss_left >= 3:
        cur_codon = chr_seq[ (ss_left-3) : ss_left ]
        codon_type = get_codon_type(cur_codon, gencode)
        #print("\t".join([str(cur_codon.translate()), str(cur_codon), codon_type]))
        if codon_type == 'coding':
            ss_left -= 3  # expand the stop-stop seq
        elif codon_type == 'stop':
            ss_left -= 3  # stop-stop seq includes stops as well
            break
        else:
            break
    return ss_left

def get_codon_type(codon, gencode, strand=1):
    acgt_only_re = re.compile('^[ACGT]+$', re.IGNORECASE)
    if not acgt_only_re.match(str(codon)):
        return 'non_ACGT'
    
    if strand == -1:
        codon = codon.reverse_complement()
    
    if codon in CodonTable.unambiguous_dna_by_id[gencode].stop_codons:
        return 'stop'
    else:
        return 'coding'

def blastn_start_end_strand(start, end):
    """Make sure the start < end
    """
    if start < end:
        return start, end, 1
    else:
        return end, start, -1

def make_translation_SeqRecord_from_CDS_feature(f, id_qualifier='locus_tag'):
    if 'translation' not in f.qualifiers:
        logging.debug("Feature " + repr(f) + "doesn't have translation!")
        return None
    prot_seq = Seq(f.qualifiers['translation'][0], generic_protein)
    
    if id_qualifier not in f.qualifiers:
        logging.debug("Feature %s doesn't have id_qualifier '%s'!" %
                      (repr(f), id_qualifier))
        return None
    
    return SeqRecord(prot_seq, id=f.qualifiers[id_qualifier][0])

def get_overlapping_feats_from_record(
    record, start, end, strand, all_types=['CDS'], min_overlap=1, max_feats=None):
    """    _overlap_len attribute will be added to returned features
    """
    target_loc = FeatureLocation(start, end, strand)
    overlapping_feats = []
    for f in record.features:
        if f.type not in all_types:
            continue
        f._overlap_len = get_FeatureLocation_overlap_len(f.location, target_loc)
        if f._overlap_len >= min_overlap:
            overlapping_feats.append(f)
    
    if max_feats is not None and len(overlapping_feats) > max_feats:
        # Sort by overlap_len: https://docs.python.org/3.6/howto/sorting.html
        overlapping_feats = sorted(overlapping_feats, reverse=True,
                                   key=lambda f: f._overlap_len)
        # Take the longest feats only
        overlapping_feats = overlapping_feats[0:max_feats]
    
    return overlapping_feats

def get_FeatureLocation_overlap_len(f1, f2):
    # https://github.com/biopython/biopython/issues/896
    if not isinstance(f1, (FeatureLocation, CompoundLocation)) or \
       not isinstance(f2, (FeatureLocation, CompoundLocation)):
        raise ValueError("Wrong feature types: %s / %s" %
                         (type(f1), type(f2)))
    
    if f1.strand is not None and \
       f2.strand is not None and \
       f1.strand != f2.strand:
        return 0
    
    return len(set(f1).intersection(set(f2)))

def fasta2dict(fn, alphabet=None):
    """TODO: Use SeqIO instead!!
    https://biopython.org/wiki/SeqIO
    from Bio import SeqIO
    record_dict = SeqIO.to_dict(SeqIO.parse("example.fasta", "fasta"))
    print(record_dict["gi:12345678"])  # use any record ID.
    """
    seq_dict = {}
    for seq_record in SeqIO.parse(fn, "fasta", alphabet):
        if seq_record.id in seq_dict:
            logging.warning("Sequence name '%s' is duplicated in file '%s'" %
                            (seq_record.id, fn))
            continue
        seq_dict[seq_record.id] = seq_record
    return seq_dict

def read_blast_xml(fn):
    f = open( fn )
    all_res = list(NCBIXML.parse(f))
    if(len(all_res) != 1):
        raise Exception('File %s must contain results for a single query!!' % fn)
    all_hits = {}
    for alignment in all_res[0].alignments:
        hits = []
        for hsp in alignment.hsps:
            hit_l, hit_r = hsp.sbjct_start-1, hsp.sbjct_end
            hits.append({
                'h_name'  : alignment.hit_def,
                'h_left'  : hit_l,
                'h_right' : hit_r,
                'h_ali'   : hsp.sbjct,
                'h_frame' : hsp.frame[1],
                'strand'  : 1 if hsp.frame[1] >= 0 else -1,
                'evalue'  : hsp.expect,
                'ali_len' : hsp.align_length,
            })
        if(len(hits) > 0):
            all_hits[alignment.hit_def] = hits
    f.close()
    return all_hits

