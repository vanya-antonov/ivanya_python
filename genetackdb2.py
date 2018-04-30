#!/usr/bin/env python3

import sys, os, re, shutil
from pprint import pprint
from warnings import warn

from Bio.SeqUtils import GC
from Bio.SeqFeature import FeatureLocation, CompoundLocation

from mylib.local import PATH2

class GeneTackDB:
    def __init__(self, budb):
        self.budb = budb
        self.gtdb_dir = PATH2['gtdb2']
    
    
    #---
    def create_new_seq_from_seq_record(self, record, org_id):   # record is Biopython's SeqRecord object: http://biopython.org/wiki/SeqRecord
        name     = record.name
        descr    = record.description
        ext_id   = record.id
        molecule = record.annotations['molecule_type']
        length   = len(record.seq)
        num_n    = len(re.compile('[^ACGTacgt]').findall(str(record.seq)))
        gc       = GC( re.compile('[^ACGTacgt]').sub('', str(record.seq)) )   # GC content without ambiguous letters
        
        fna_path = None
        warn('SAVE SEQ.fna!!!!!!!!!!!!!!!!!1')
        
        seq_id = self.budb.get_random_db_id('seqs', 'seq_id')
        self.budb.exec_sql_in('''INSERT INTO seqs (
            seq_id, org_id, name, descr, molecule, ext_id, len,    num_n, gc, fna_path)
        ''',seq_id, org_id, name, descr, molecule, ext_id, length, num_n, gc, fna_path
        )
        
        for f in record.features:
            self.create_new_sfeat(seq_id, f)
        sys.exit()
        
        return seq_id
    
    #---
    def create_new_sfeat(self, seq_id, f):
        feature   = f.type
        start     = int(f.location.start)
        end       = int(f.location.end)
        strand    = '+' if f.location.strand == 1 else '-'
        gene      = f.qualifiers.get('gene',      [None])[0]
        locus_tag = f.qualifiers.get('locus_tag', [None])[0]
        descr     = f.qualifiers.get('product',   [None])[0]
        slippage  = 1 if 'ribosomal_slippage' in f.qualifiers else 0
        pseudo    = 1 if 'pseudo'             in f.qualifiers else 0
        
        sfeat_id = self.budb.get_random_db_id('sfeats', 'sfeat_id')
        self.budb.exec_sql_in('''INSERT INTO sfeats (
            sfeat_id, seq_id, feature, start, end, strand, gene, locus_tag, descr, slippage, pseudo)
        ''',sfeat_id, seq_id, feature, start, end, strand, gene, locus_tag, descr, slippage, pseudo
        )
        
        # SFEAT_PARAMS
        prm_keys = ['note', 'transl_table', 'translation', 'protein_id',
            'EC_number', 'gene_synonym', 'function', 'regulatory_class', 'ncRNA_class'
        ]
        for k, vals in f.qualifiers.items():
            if k in prm_keys:
                for v in vals:
                    self.add_param_to('sfeat', sfeat_id, k, v)
        
        dbx_ll   = [s.split(':') for s in f.qualifiers.get('db_xref', [])]   # list of lists
        dbx_dict = {l[0] : l[1] for l in dbx_ll}
        for (name,value) in dbx_dict.items():
            self.add_param_to('sfeat', sfeat_id, name, value)
        
        if slippage == 1:
            coord = _get_slippage_coord(f)
            if coord is not None:
                self.add_param_to('sfeat', sfeat_id, 'ribosomal_slippage', str(f.location), coord)
    
    #---
    def add_param_to(self, unit, id, name, value, num=None):
        sql = 'INSERT INTO {}_params ({}_id, name, value, num)'.format(unit,unit)
        self.budb.exec_sql_in(sql, id, name, value, num)
    
    #---
    # The function returns new file path RELATIVE to GTDB dir
    def save_file(self, from_fn, to_dir, to_name):
        to_path = os.path.join(self.gtdb_dir, to_dir)
        if not os.path.exists(to_path):
            os.makedirs(to_path)
        
        to_fn = os.path.join(to_path, to_name)
        if os.path.exists(to_fn):
            raise Exception("File '{}' already exists!".format(to_fn))
        
        return os.path.join(to_dir, to_name)
    
    #---
    # Argument is the path RELATIVE to the GTDB dir
    def delete_file(self, fn):
        fn = os.path.join(self.gtdb_dir, fn)
        if os.path.exists(fn):
            os.remove(fn)
        else:
            warn("File '%s' doesn't exist!" % fn)
    
    #---
    # Argument is the path RELATIVE to the GTDB dir
    def delete_dir(self, path):
        path = os.path.join(self.gtdb_dir, path)
        if os.path.exists(path):
            shutil.rmtree(path)
        else:
            warn("Folder '%s' doesn't exist!" % path)
    


class Org:
    
    #---
    def __init__(self, gtdb, info):
        self.gtdb = gtdb
        
        # use dictionary keys as object arguments: https://stackoverflow.com/a/2466207/310453
        for k, v in info.items():
            setattr(self, k, v)
    
    #---
    # I want to be able to return None: https://stackoverflow.com/a/25200825/310453
    @classmethod
    def by_db_entry_id(cls, gtdb, org_id):
        res = gtdb.budb.exec_sql_ar('''
            SELECT org_id, name, genus, phylum, kingdom, dir_path
            FROM orgs WHERE org_id=%s
        ''', org_id)
        if len(res) == 0:
            return None
        else:
            return cls(gtdb, res[0])
    
    #---
    def delete_from_db(self):
        if self.dir_path is not None:
            self.gtdb.delete_dir(self.dir_path)
        self.gtdb.budb.exec_sql_nr('DELETE FROM org_params WHERE org_id=%s', self.org_id)
        self.gtdb.budb.exec_sql_nr('DELETE FROM orgs       WHERE org_id=%s', self.org_id)
    
    #---
    # https://stackoverflow.com/a/12179752/310453
    @staticmethod
    def get_db_entry_id_for_seq_record(gtdb, record):   # record is Biopython's SeqRecord object: http://biopython.org/wiki/SeqRecord
        dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        asm_id   = dbx_dict.get('Assembly')
        if asm_id is None:
            raise Exception("Can't determine species: sequence '%s' doesn't have Assembly ID!" % record.id)
        org_l = gtdb.budb.exec_sql_ar('select org_id from org_params where name=%s and value=%s', 'Assembly', asm_id)
        
        if len(org_l) == 0:
            return None
        if len(org_l) > 1:
            warn("Several different orgs were found for seq '%s' -- the first one is used" % record.id)
        return org_l[0]['org_id']
    
    #---
    @staticmethod
    def create_new_db_entry_from_seq_record(gtdb, record):   # record is Biopython's SeqRecord object: http://biopython.org/wiki/SeqRecord
        name = record.annotations['organism']
        
        genus, phylum, kingdom = _get_genus_phylum_kingdom(record)
        
        dir_path = Org._create_org_dir(gtdb, name)
        
        org_id = gtdb.budb.get_random_db_id('orgs', 'org_id')
        gtdb.budb.exec_sql_in('''INSERT INTO orgs (
            org_id, name, genus, phylum, kingdom, dir_path)
        ''',org_id, name, genus, phylum, kingdom, dir_path
        )
        
        # ORG_PARAMS
        dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        for (name,value) in dbx_dict.items():
            gtdb.add_param_to('org', org_id, name, value)
        
        num = 0
        for value in record.annotations.get('taxonomy',[]):
            gtdb.add_param_to('org', org_id, 'taxonomy', value, num)
            num += 1
        return org_id
    
    #---
    # The function returns new path RELATIVE to the GTDB root dir
    @staticmethod
    def _create_org_dir(gtdb, org_name):
        # 'Natranaerobius thermophilus JW/NM-WN-LF'  =>  'Natranaerobius_thermophilus_JW_NM_WN_LF'
        folder = re.compile('[^\w\d]+').sub('_', org_name).strip('_')
        dir_path = os.path.join(gtdb.gtdb_dir, 'orgs', folder)
        if os.path.exists(dir_path):
            raise Exception("Folder '%s' already exists!" % dir_path)
        os.makedirs(dir_path)
        return os.path.join('orgs', folder)

#---
def db_xref_list_to_dict(dbxrefs):
    dbx_ll = [str.split(':') for str in dbxrefs]   # list of lists
    return {l[0] : l[1] for l in dbx_ll}

#---
def _get_slippage_coord(f):
    if not isinstance(f.location, CompoundLocation):
        warn("Can't determine slippage coord: '{}' is not a CompoundLocation!".format(str(f.location)))
        return None
    if len(f.location.parts) != 2:
        warn("Can't determine slippage coord: CompoundLocation '{}' has '{}' parts!".format(str(f.location), len(f.location.parts)))
        return None
    
    if f.location.strand == 1:
        # join(3990432..3990762,3990762..3991240)
        return int(f.location.parts[0].end)
    else:
        # complement(join(3413273..3413751,3413751..3414081))
        return int(f.location.parts[1].start)

#---
def _get_genus_phylum_kingdom(seq_record):
    taxa_l = seq_record.annotations.get('taxonomy')
    if taxa_l is None or len(taxa_l) < 3:
        return None, None, None
    
    kingdom = taxa_l[0]
    if kingdom not in ['Bacteria', 'Archaea', 'Eukaryota']:
        warn("Unknown kingdom '%s'" % kingdom)
        return None, None, None
    
    species = seq_record.annotations['organism']
    genus   = taxa_l[-1]
    if genus not in species:
        warn("Can't find genus '%s' in species name '%s'" % (genus, species))
        genus = None
    
    phylum = taxa_l[1]
    return genus, phylum, kingdom


