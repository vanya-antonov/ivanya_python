
# $Id$

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

