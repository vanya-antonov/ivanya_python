#!/usr/bin/env python3

import logging
import sys
import os
import re
import shutil
from pprint import pprint
#from warnings import warn

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from mylib.baseutildb import BaseUtilDB
from mylib.local import PATH2


class GeneTackDB(BaseUtilDB):
    
    #---
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.gtdb_dir = PATH2['gtdb2']
    
    #---
    def add_param_to(self, unit, db_id, name, value=None, num=None):
        sql = 'INSERT INTO {}_params ({}_id, name, value, num)'.format(unit,unit)
        self.exec_sql_in(sql, db_id, name, value, num)
    
    #---
    # Argument is the path RELATIVE to the GTDB dir
    def delete_file(self, fn):
        fn = os.path.join(self.gtdb_dir, fn)
        if os.path.exists(fn):
            os.remove(fn)
        else:
            logging.warning("File '%s' doesn't exist!" % fn)
    
    #---
    # Argument is the path RELATIVE to the GTDB dir
    def delete_dir(self, path):
        path = os.path.join(self.gtdb_dir, path)
        if os.path.exists(path):
            shutil.rmtree(path)
        else:
            logging.warning("Folder '%s' doesn't exist!" % path)


class TableObject:
    
    def __init__(self, gtdb, db_id, main_sql, unit=None, add_prm=False,
                 prm_value=[], prm_int=[], prm_float=[], prm_list=[]):
        self.gtdb = gtdb
        self.prm = {}
        
        res = gtdb.exec_sql_ar(main_sql, db_id)
        if len(res) == 0:
            raise Exception("Unknown db_id = '{}'".format(db_id))
        
        # use dictionary keys as object arguments: https://stackoverflow.com/a/2466207/310453
        for k, v in res[0].items():
            setattr(self, k, v)
        
        # Load params if requested
        if add_prm:
            prm_sql = 'select name, value, num from {0}_params where {0}_id=%s'.format(unit)
            prm_res = gtdb.exec_sql_ar(prm_sql, db_id)
            for d in prm_res:
                if d['name'] in prm_list:
                    if d['name'] not in self.prm:
                        self.prm[d['name']] = []
                    self.prm[d['name']].append(d)
                elif d['name'] in self.prm:
                    raise Exception("Param name '{}' is duplicated, but not in prm_list".format(d['name']))
                
                if d['name'] in prm_value:
                    self.prm[d['name']] = d['value']
                elif d['name'] in prm_int:
                    self.prm[d['name']] = int(d['num'])
                elif d['name'] in prm_float:
                    self.prm[d['name']] = float(d['num'])
                else:
                    self.prm[d['name']] = d


class Org(TableObject):
    
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id, """
            SELECT id, c_date, name, genus, phylum, kingdom, dir_path
            FROM orgs WHERE id=%s
        """)
    
    def delete_from_db(self):
        for seq in self.gtdb.exec_sql_ar('select id from seqs where org_id=%s', self.id):
            Seq(self.gtdb, seq['id']).delete_from_db()
        self.gtdb.exec_sql_nr('DELETE FROM org_params WHERE org_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM orgs       WHERE     id=%s', self.id)
        if self.dir_path is not None:
            self.gtdb.delete_dir(self.dir_path)
    
    #---
    # https://stackoverflow.com/a/12179752/310453
    @staticmethod
    def get_db_id_by_SeqRecord(gtdb, record):
        dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        asm_id   = dbx_dict.get('Assembly')
        if asm_id is None:
            raise Exception("Can't determine species: sequence '%s' doesn't have Assembly ID!" % record.id)
        org_l = gtdb.exec_sql_ar('SELECT org_id FROM org_params WHERE name=%s AND value=%s', 'Assembly', asm_id)
        
        if len(org_l) == 0:
            return None
        return org_l[0]['org_id']
    
    #---
    @staticmethod
    def create_new_in_db_from_SeqRecord(gtdb, record, user_id):
        name = record.annotations['organism']
        
        genus, phylum, kingdom = _get_genus_phylum_kingdom(record)
        
        dir_path = Org._create_org_dir(gtdb, name)
        
        db_id = gtdb.get_random_db_id('orgs', 'id')
        gtdb.exec_sql_in("""INSERT INTO orgs (
                id, user_id, name, genus, phylum, kingdom, dir_path)
        """, db_id, user_id, name, genus, phylum, kingdom, dir_path
        )
        
        # ORG_PARAMS
        dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        for (name,value) in dbx_dict.items():
            gtdb.add_param_to('org', db_id, name, value)
        
        num = 0
        for value in record.annotations.get('taxonomy',[]):
            gtdb.add_param_to('org', db_id, 'taxonomy', value, num)
            num += 1
        return db_id
    
    #---
    # The function returns new path RELATIVE to the GTDB root dir
    @staticmethod
    def _create_org_dir(gtdb, org_name):
        # 'Natranaerobius thermophilus JW/NM-WN-LF'  =>  'Natranaerobius_thermophilus_JW_NM_WN_LF'
        folder = re.compile('[^\w\d]+').sub('_', org_name).strip('_')
        subdir = os.path.join(gtdb.gtdb_dir, 'orgs')
        full_path = make_new_fullpath_for_basename(subdir, folder)
        os.makedirs(full_path)
        return os.path.relpath(full_path, gtdb.gtdb_dir)


class Seq(TableObject):
    
    def __init__(self, gtdb, db_id):
        main_sql = """
            SELECT id, user_id, c_date, org_id, name, descr, type, ext_id, len
            FROM seqs WHERE id=%s
        """
        super().__init__(
            gtdb, db_id, main_sql, unit='seq', add_prm=True,
            prm_value = ['fna_path', 'gbk_path'],
            prm_float = ['gc'],
            prm_int = ['num_n'],
        )
    
    #---
    def delete_from_db(self):
        for d in self.gtdb.exec_sql_ar('select id from fsgenes where seq_id=%s', self.id):
            FSGene(self.gtdb, d['id']).delete_from_db()
        for d in self.gtdb.exec_sql_ar('select id from sfeats where seq_id=%s', self.id):
            SFeat(self.gtdb, d['id']).delete_from_db()
        self.gtdb.exec_sql_nr('DELETE FROM seq_params WHERE seq_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM seqs       WHERE     id=%s', self.id)
        logging.warning('Delete files is NOT implemented!!')
#        if self.fna_path is not None:
#            self.gtdb.delete_file(self.fna_path)
    
    #---
    # https://stackoverflow.com/a/12179752/310453
    @staticmethod
    def get_db_id_by_ext_id(gtdb, ext_id):
        res_l = gtdb.exec_sql_ar('SELECT id FROM seqs WHERE ext_id=%s', ext_id)
        
        if len(res_l) == 0:
            return None
        elif len(res_l) > 1:
            raise Exception("EXT_ID '{}' matches several records in DB!".format(record.id))
        else:
            return res_l[0]['id']
    
    #---
    @staticmethod
    def create_new_in_db_from_SeqRecord(gtdb, record, org, user_id):
        name   = record.name
        descr  = record.description
        ext_id = record.id
        m_type = record.annotations['molecule_type']
        length = len(record.seq)
        
        db_id = gtdb.get_random_db_id('seqs', 'id')
        
        gtdb.exec_sql_in("""INSERT INTO seqs (
               id, user_id, org_id, name, descr,   type, ext_id, len)
        """,db_id, user_id, org.id, name, descr, m_type, ext_id, length
        )
        
        # Seq params
        fna_path = Seq._save_fasta(gtdb, record, org, db_id)
        gtdb.add_param_to('seq', db_id, 'fna_path', value=fna_path)
        
        gbk_path = Seq._save_gbk(gtdb, record, org)
        gtdb.add_param_to('seq', db_id, 'gbk_path', value=gbk_path)
        
        num_n = len(re.compile('[^ACGTacgt]').findall(str(record.seq)))
        gtdb.add_param_to('seq', db_id, 'num_n', num=num_n)
        
        # GC content without ambiguous letters
        gc = GC(re.compile('[^ACGTacgt]').sub('', str(record.seq)))
        gc_str = '%.1f%%' % gc   # 72.1186611098208   => ' 72.1%'
        gtdb.add_param_to('seq', db_id, 'gc', value=gc_str, num=gc)
        
        return db_id
    
    @staticmethod
    def _save_gbk(gtdb, record, org):
        """The function returns new path RELATIVE to the GTDB root dir"""
        subdir = os.path.join(gtdb.gtdb_dir, org.dir_path, 'seq_gbk')
        full_path = make_new_fullpath_for_basename(subdir, record.id)
        SeqIO.write(record, full_path, "genbank")
        return os.path.relpath(full_path, gtdb.gtdb_dir)
    
    @staticmethod
    def _save_fasta(gtdb, record, org, seq_id):
        """The function returns new path RELATIVE to the GTDB root dir"""
        fna_dir = os.path.join(gtdb.gtdb_dir, org.dir_path, 'seq_fna')
        full_path = make_new_fullpath_for_basename(fna_dir, record.id)
        
        # Create tmp SeqRecord to have just '>seq_in' inside the fasta file
        tmp = SeqRecord(record.seq, id=str(seq_id), name='', description='')
        SeqIO.write(tmp, full_path, "fasta")
        return os.path.relpath(full_path, gtdb.gtdb_dir)


class SFeat(TableObject):
    
    #---
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id, """
            SELECT id, seq_id, type, start, end, strand, name, descr, ext_id
            FROM sfeats WHERE id=%s
        """)
    
    #---
    def delete_from_db(self):
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_sfeats WHERE sfeat_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM sfeat_params  WHERE sfeat_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM sfeats        WHERE       id=%s', self.id)
    
    #---
    @staticmethod
    def create_new_in_db_from_SeqFeature(gtdb, f, seq_id):
        start   = int(f.location.start)
        end     = int(f.location.end)
        strand  = f.location.strand
        name    = f.qualifiers.get('gene',      [None])[0]
        descr   = f.qualifiers.get('product',   [None])[0]
        ext_id  = f.qualifiers.get('locus_tag', [None])[0]
        
        db_id = gtdb.get_random_db_id('sfeats', 'id')
        gtdb.exec_sql_in("""INSERT INTO sfeats (
               id, seq_id,   type, start, end, strand, name, descr, ext_id)
        """,db_id, seq_id, f.type, start, end, strand, name, descr, ext_id
        )
        
        # SFEAT_PARAMS
        prm_keys = [
            'translation', 'protein_id', 'ribosomal_slippage', 'pseudo', 'experiment',
            'EC_number', 'gene_synonym', 'function', 'regulatory_class', 'ncRNA_class',
        ]
        for k, vals in f.qualifiers.items():
            if k in prm_keys:
                for v in vals:
                    gtdb.add_param_to('sfeat', db_id, k, v)
        
        dbx_dict = db_xref_list_to_dict(f.qualifiers.get('db_xref', []))
        for (name,value) in dbx_dict.items():
            gtdb.add_param_to('sfeat', db_id, name, value)
        
        return db_id


class FSGene(TableObject):
    
    #---
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id, """
            SELECT id, user_id, seq_id, fs_coord, type, start, end, strand, source
            FROM fsgenes WHERE id=%s
        """)
    
    #---
    def delete_from_db(self):
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_sfeats WHERE fsgene_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_params WHERE fsgene_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM fsgenes       WHERE        id=%s', self.id)
    
    #---
    @staticmethod
    def create_new_in_db_from_SeqFeature(gtdb, f, seq_id, user_id, source=None):
        error = None
        if not isinstance(f.location, CompoundLocation):
            error = "Feature location is not a CompoundLocation object: '{}'".format(str(f.location))
        if len(f.location.parts) != 2:
            error = "CompoundLocation '{}' has '{}' parts!".format(str(f.location), len(f.location.parts))
        if f.location.operator != 'join':
            error = "Wrong feature location operator = '{}'".format(f.location.operator)
        if f.location.strand is None:
            error = "The parts of CompoundLocation have different strands = '{}'".format(f.location)
        if error is not None:
            logging.warning(error)
            return None
        
        if f.location.parts[0].start < f.location.parts[1].start:
            left_part, right_part = f.location.parts[0], f.location.parts[1]
        else:
            left_part, right_part = f.location.parts[1], f.location.parts[0]
        
        start1, end1 = int(left_part.start),  int(left_part.end)
        start2, end2 = int(right_part.start), int(right_part.end)
        fs_coord = end1 if f.strand == 1 else start2
        fs_type = start2 - end1
        if fs_type != 0 and fs_type % 3 == 0:
            logging.warning("fs_type = '{}' for feature {} is divisible by 3!".format(fs_type, f))
        
        # make name like: 'NC_010002.1:1922457:-1'
        seq = Seq(gtdb, seq_id)
        name = '{}:{}:{:+d}'.format(seq.ext_id, fs_coord, fs_type)
        
        db_id = gtdb.get_random_db_id('fsgenes', 'id')
        gtdb.exec_sql_in("""INSERT INTO fsgenes (
               id, user_id, seq_id, name, fs_coord,    type, start,  end,    strand, source)
        """,db_id, user_id, seq_id, name, fs_coord, fs_type, start1, end2, f.strand, source
        )
        
        return db_id
    
    #---
    @staticmethod
    def set_fsgene_sfeats(gtdb, fsgene_id, sfeat_list):
        
        # Remove old sfeat links
        gtdb.exec_sql_nr('delete from fsgene_sfeats where fsgene_id=%s', fsgene_id)
        
        sfeat_descr_l = []
        for sfeat_id in sfeat_list:
            gtdb.exec_sql_in('INSERT INTO fsgene_sfeats (fsgene_id, sfeat_id)', fsgene_id, sfeat_id)
            sfeat = SFeat(gtdb, sfeat_id)
            if sfeat.descr is not None:
                sfeat_descr_l.append(sfeat.descr)
        
        # Update fsgene descr
        fsgene_descr = None
        if len(sfeat_descr_l) > 0:
            fsgene_descr = '|'.join(sfeat_descr_l)
        gtdb.exec_sql_nr('update fsgenes set descr=%s where id=%s', fsgene_descr, fsgene_id)



def make_new_fullpath_for_basename(subdir, basename):
    """Returns a full path to create a new file or folder"""
    # create subdir if needed
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    
    # Make sure the path doesn't exist yet
    num, suffix = 1, ''
    full_path = os.path.join(subdir, basename + suffix)
    while os.path.exists(full_path):
        logging.warning("Path '{}' already exists!".format(full_path))
        suffix = '.' + str(num)
        full_path = os.path.join(subdir, basename + suffix)
        num += 1
    return full_path

#---
def db_xref_list_to_dict(dbxrefs):
    dbx_ll = [str.split(':') for str in dbxrefs]   # list of lists
    return {l[0] : l[1] for l in dbx_ll}

#---
def _get_genus_phylum_kingdom(seq_record):
    taxa_l = seq_record.annotations.get('taxonomy')
    if taxa_l is None or len(taxa_l) < 3:
        return None, None, None
    
    kingdom = taxa_l[0]
    if kingdom not in ['Bacteria', 'Archaea', 'Eukaryota']:
        logging.warning("Unknown kingdom '%s'" % kingdom)
        return None, None, None
    
    species = seq_record.annotations['organism']
    if taxa_l[-1] in species:
        genus = taxa_l[-1]
    elif taxa_l[-2] in species:
        genus = taxa_l[-2]
    else:
        logging.warning("Can't determine genus for species name '{}' and taxonomy {}".format(species, taxa_l))
        genus = None
    
    phylum = taxa_l[1]
    return genus, phylum, kingdom


