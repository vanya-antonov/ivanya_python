#!/usr/bin/env python3

from mylib.genetackdb2 import Org, FSGene

class ChelOrg(Org):
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id)
        self.load_prm(prm_str=['chel_genotype_LMS', 'chel_genotype_NHDI',
                              'chel_genotype_NTS_HDI'])

    def get_true_M_chelatase_fsgene_ids(self):
        return [d['id'] for d in self.gtdb.exec_sql_ar(
            '''select distinct fs.id from seqs s, fsgenes fs
            where s.org_id=%s and fs.seq_id=s.id and fs.fs_type <> 0
            ''', self.id)]

class ChelFSGene(FSGene):
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id)
        self.load_prm(prm_str=['chel_gene', 'chel_subunit'])
    
