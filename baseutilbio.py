
# $Id$

from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation

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

def get_overlapping_feats_from_record(record, start, end, strand,
                                      all_types=['CDS'], min_overlap=1,
                                      max_feats=None):
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
                'strand'  : '+' if hsp.frame[1] >= 0 else '-',
                'evalue'  : hsp.expect,
                'ali_len' : hsp.align_length,
            })
        if(len(hits) > 0):
            all_hits[alignment.hit_def] = hits
    f.close()
    return all_hits

