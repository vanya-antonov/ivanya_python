
# $Id$

from Bio.Blast import NCBIXML
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def get_FeatureLocation_overlap_len(f1, f2):
    # https://github.com/biopython/biopython/issues/896
    if not isinstance(f1, (FeatureLocation, CompoundLocation)) or \
       not isinstance(f2, (FeatureLocation, CompoundLocation)):
        raise ValueError("The function just works on FeatureLocation and CompoundLocation objects")
    
    if f1.strand is not None and \
       f2.strand is not None and \
       f1.strand != f2.strand:
        return 0
    
    return len(set(f1).intersection(set(f2)))

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

