#!/usr/bin/env python3

from mylib.genetackdb2 import Org

class ChelOrg(Org):
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id)
        self.load_prm(prm_str=['chel_genotype_LMS', 'chel_genotype_NHDI',
                              'chel_genotype_NTS_HDI'])

