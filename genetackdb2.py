#!/usr/bin/env python3

import logging
import sys
import os
import re
import shutil
from pprint import pprint
from collections import Counter

import Bio.Alphabet
import Bio.Seq
from Bio import SeqIO
#from Bio.Seq import Seq as BioSeq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from mylib.baseutildb import BaseUtilDB
from mylib.local import PATH2


class GeneTackDB(BaseUtilDB):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.gtdb_dir = PATH2['gtdb2']
    
    def add_param_to(self, unit, db_id, name, value=None, num=None):
        sql = 'INSERT INTO {0}_params ({0}_id, name, value, num)'.format(unit)
        self.exec_sql_in(sql, db_id, name, value, num)
    
    def delete_param(self, unit, db_id, name):
        sql = 'delete from {0}_params where {0}_id=%s and name=%s'.format(unit)
        self.exec_sql_nr(sql, db_id, name)
    
    # Argument is the path RELATIVE to the GTDB dir
    def delete_file(self, fn):
        fn = os.path.join(self.gtdb_dir, fn)
        if os.path.exists(fn):
            os.remove(fn)
        else:
            logging.warning("File '%s' doesn't exist!" % fn)
    
    # Argument is the path RELATIVE to the GTDB dir
    def delete_dir(self, path):
        path = os.path.join(self.gtdb_dir, path)
        if os.path.exists(path):
            shutil.rmtree(path)
        else:
            logging.warning("Folder '%s' doesn't exist!" % path)


class TableObject:
    
    def __init__(self, gtdb, db_id, main_sql, add_prm=False,
                 prm_str=[], prm_int=[], prm_float=[], prm_list=[]):
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
            prm_sql = 'select name, value, num from {0}_params where {0}_id=%s'.format(self._unit)
            prm_res = gtdb.exec_sql_ar(prm_sql, db_id)
            for d in prm_res:
                if d['name'] in prm_list:
                    if d['name'] not in self.prm:
                        self.prm[d['name']] = []
                    self.prm[d['name']].append(d)
                elif d['name'] in self.prm:
                    raise Exception("Param name '{}' is duplicated, but not in prm_list".format(d['name']))
                
                if d['name'] in prm_str:
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


class Seq(TableObject, Bio.Seq.Seq):
    
    def __init__(self, gtdb, db_id, read_seq=False, use_seq=None):
        """Create Seq corresponding to the SEQ GeneTack DB table

        Arguments:
         - gtdb - [reqiured] GeneTackDB object
         - db_id - [required] valid ID from the SEQ table
         - read_seq - [optional] read the sequence from a DB file
         - use_seq - [optional] take sequence from this Bio.Seq.Seq object
        """
        self._unit = 'seq'
        self._has_seq = False
        main_sql = """
            SELECT id, user_id, c_date, org_id, name, descr, type, ext_id, len
            FROM seqs WHERE id=%s
        """
        TableObject.__init__(
            self, gtdb, db_id, main_sql,
            add_prm = True,
            prm_str = ['fna_path', 'gbk_path'],
            prm_int = ['num_n', 'transl_table'],
            prm_float = ['gc'],
        )
        
        # Initialize the Bio.Seq if requested
        if use_seq is not None:
            Bio.Seq.Seq.__init__(self, str(use_seq), use_seq.alphabet)
        elif read_seq:
            fasta_fn = os.path.join(gtdb.gtdb_dir, self.prm['fna_path'])
            record = SeqIO.read(fasta_fn, "fasta")
            alphabet = Bio.Alphabet.SingleLetterAlphabet
            if self.type == 'DNA' and self.prm['num_n'] == 0:
                alphabet = Bio.Alphabet.IUPAC.unambiguous_dna
            elif self.type == 'DNA' and self.prm['num_n'] > 0:
                alphabet = Bio.Alphabet.IUPAC.ambiguous_dna
            Bio.Seq.Seq.__init__(self, str(record.seq), alphabet)
        
        # self._data is created by the Bio.Seq.Seq constructor
        if hasattr(self, '_data'):
            # make sure that len in DB and the actual sequence length match
            if self.len != len(self):
                raise Exception("Sequence length ('%d') != value in DB ('%d')", (self.len, len(self)))
            else:
                self._has_seq = True
    
    
    def delete_from_db(self):
        for d in self.gtdb.exec_sql_ar('select id from fsgenes where seq_id=%s', self.id):
            FSGene(self.gtdb, d['id']).delete_from_db()
        for d in self.gtdb.exec_sql_ar('select id from sfeats where seq_id=%s', self.id):
            SFeat(self.gtdb, d['id']).delete_from_db()
        self.gtdb.exec_sql_nr('DELETE FROM seq_params WHERE seq_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM seqs       WHERE     id=%s', self.id)
        
        all_files = [self.prm.get('fna_path'), self.prm.get('gbk_path')]
        for fn in filter(lambda fn: fn is not None, all_files):
            self.gtdb.delete_file(fn)
    
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
        
        # SEQ_PARAMS
        fna_path = Seq._save_fasta(gtdb, record, org, db_id)
        gtdb.add_param_to('seq', db_id, 'fna_path', value=fna_path)
        
        gbk_path = Seq._save_gbk(gtdb, record, org)
        gtdb.add_param_to('seq', db_id, 'gbk_path', value=gbk_path)
        
        # Get translation table
        transl_table = _get_genetic_code_from_SeqRecord(record)
        if transl_table is not None:
            gtdb.add_param_to('seq', db_id, 'transl_table', num=transl_table)
        
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
    
    def __init__(self, gtdb, db_id):
        super().__init__(gtdb, db_id, """
            SELECT id, seq_id, type, start, end, strand, name, descr, ext_id
            FROM sfeats WHERE id=%s
        """)
    
    def delete_from_db(self):
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_sfeats WHERE sfeat_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM sfeat_params  WHERE sfeat_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM sfeats        WHERE       id=%s', self.id)
    
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
    def __init__(self, gtdb, db_id):
        self._unit = 'fsgene'
        main_sql = """
            SELECT id, user_id, seq_id, fs_coord, type, start, end, strand, source
            FROM fsgenes WHERE id=%s
        """
        TableObject.__init__(
            self, gtdb, db_id, main_sql,
            add_prm = True,
            prm_str = ['prot_seq_n', 'prot_seq_c', 'nt_seq_n', 'nt_seq_c'],
            prm_int = [],
            prm_float = [],
        )
    
    def delete_from_db(self):
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_sfeats WHERE fsgene_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_params WHERE fsgene_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM fsgenes       WHERE        id=%s', self.id)
    
    def set_sfeats(self, sfeat_list):
        # Remove old sfeat links
        self.gtdb.exec_sql_nr('delete from fsgene_sfeats where fsgene_id=%s', self.id)
        
        sfeat_descr_l = []
        for sfeat_id in sfeat_list:
            self.gtdb.exec_sql_in('INSERT INTO fsgene_sfeats (fsgene_id, sfeat_id)', self.id, sfeat_id)
            sfeat = SFeat(self.gtdb, sfeat_id)
            if sfeat.descr is not None:
                sfeat_descr_l.append(sfeat.descr)
        
        # Update fsgene descr
        fsgene_descr = None
        if len(sfeat_descr_l) > 0:
            fsgene_descr = ' | '.join(sfeat_descr_l)
        self.gtdb.exec_sql_nr('update fsgenes set descr=%s where id=%s', fsgene_descr, self.id)
    
    def make_all_params(self, seq=None):
        if seq is None:
            seq = Seq(self.gtdb, self.seq_id, read_seq=True)
        self.make_prm_seqs(seq)
        
    def make_prm_seqs(self, seq):
        seq_dict = FSGene._get_prot_seq_parts(
            self.start, self.end, self.strand, self.fs_coord, self.type, seq
        )
        
        prm_names = ['prot_seq_n', 'prot_seq_c', 'nt_seq_n', 'nt_seq_c']
        for name in prm_names:
            self.gtdb.delete_param(self._unit, self.id, name)
            self.gtdb.add_param_to(self._unit, self.id, name, value=seq_dict[name])
    
    def _get_prot_seq_parts(start, end, strand, fs_coord, fs_type, seq):
        up_len_nt = fs_coord - start
        down_len_nt = end - fs_coord
        if strand == -1:
            up_len_nt, down_len_nt = down_len_nt, up_len_nt
        
        # Backward frameshiting increases the length of the C-terminal part of the fs-protein
        down_len_nt -= fs_type
        
        if up_len_nt % 3 != 0 or down_len_nt % 3 != 0:
            fmt = """
            The length (%d) of the upstream / downstream fsgene part (%d / %d) is not divisible by 3:
            start=%d, fs_coord=%d, end=%d, strand=%d
            """
            raise Exception(fmt % (up_len_nt, down_len_nt, start, fs_coord, end, strand))
        
        chunk_nt = seq[start:end]
        if strand == -1:
            chunk_nt = chunk_nt.reverse_complement()
        
        up_chunk_nt = chunk_nt[:up_len_nt]
        down_chunk_nt = chunk_nt[-down_len_nt:]
        prot_seq_n = up_chunk_nt.translate(table = seq.prm['transl_table'])
        prot_seq_c = down_chunk_nt.translate(table = seq.prm['transl_table'])
        prot_seq_c = prot_seq_c.rstrip('*')  # remove possible stop codon at the end
        
        if '*' in prot_seq_n or '*' in prot_seq_c:
            raise Exception("FS-prot seq contains in-frame stop codon:\n%s\n\n%s", (prot_seq_n, prot_seq_c))
        
        return {
            'nt_seq_n': up_chunk_nt.upper(), 'nt_seq_c': down_chunk_nt.upper(),
            'prot_seq_n': prot_seq_n.upper(), 'prot_seq_c': prot_seq_c.upper(),
        }
    
    @staticmethod
    def create_new_in_db_from_GT_FS(gtdb, fs_id, source='genetack'):
        fs, seq = FSGene._get_info_about_GT_FS(gtdb, fs_id)
        
        # If up_len_nt is not divisible by 3 => the FS was predicted in a middle of codon.
        # Move the fs-coord upstream (not downstrem) to avoid stop codon in cases like 'aaa_tgA_CCC'.
        adjust_len = fs['up_len_nt'] % 3
        if fs['strand'] == 1:
            fs['fs_coord'] -= adjust_len
        else:
            fs['fs_coord'] += adjust_len
        
        fsgene_id = FSGene.create_new_in_db(
            gtdb, seq, fs['user_id'], fs['start'], fs['end'], fs['strand'], fs['fs_coord'], fs['type'],
            source=source, cof_id=fs['cof_id'], c_date=fs['c_date'], db_id=fs_id
        )
        
        FSGene(gtdb, fsgene_id).make_all_params(seq)
        
        return fsgene_id
    
    @staticmethod
    def create_new_in_db(gtdb, seq, user_id, start, end, strand, fs_coord, fs_type,
                         source=None, cof_id=None, c_date=None, db_id=None):
        # make name like: 'NC_010002.1:1922457:-1'
        name = '%s:%d:%+d' % (seq.ext_id, fs_coord, fs_type)
        
        if db_id is None:
            db_id = gtdb.get_random_db_id('fsgenes', 'id')
        
        gtdb.exec_sql_in("""INSERT INTO fsgenes (
               id, c_date, user_id, seq_id, name, fs_coord,    type, start, end, strand, source, cof_id)
        """,db_id, c_date, user_id, seq.id, name, fs_coord, fs_type, start, end, strand, source, cof_id
        )
        
        return db_id
    
    @staticmethod
    def _get_info_about_GT_FS(gtdb, fs_id):
        res = gtdb.exec_sql_ar("""
            select f.fs_id, j.user_id, j.c_date, s.id AS seq_id, s.ext_id AS seq_ext_id,
                   f.fs_coord, f.type, f.init_gene_seq,
                   IF(f.strand = "+", 1, -1) AS strand,
                   (select cof_id from cof_gtfs cg where cg.fs_id = f.fs_id limit 1) AS cof_id
            from gt_fs f, jobs j, seqs s
            where j.job_id=f.job_id
            and s.name = j.name
            and f.fs_id = %s
        """, fs_id
        )
        if len(res) == 0:
            logging.warning("Can't get info about FS_ID = '%s'" % fs_id)
            return None
        elif not res[0]['seq_ext_id'].endswith('.1'):
            logging.warning("Sequence '%s' has wrong version!" % res[0]['seq_ext_id'])
            return None
        else:
            fs = res[0]
            fs['type'] = int(fs['type'])
        
        if fs['strand'] == -1:
            fs['fs_coord'] -= 1   # switch to zero-based coordinate system?
        
        up_match = re.compile('^[a-z]+').search(fs['init_gene_seq'])
        down_match = re.compile('[A-Z]+$').search(fs['init_gene_seq'])
        if up_match is None or down_match is None:
            logging.warning("Wrong fsgene sequence '%s'" % fs['init_gene_seq'])
            return None
        fs['up_len_nt'] = up_match.end() - up_match.start()
        fs['down_len_nt'] = down_match.end() - down_match.start()
        
        if fs['strand'] == 1:
            fs['start'] = fs['fs_coord'] - fs['up_len_nt']
            fs['end'] = fs['fs_coord'] + fs['down_len_nt']
        else:
            fs['start'] = fs['fs_coord'] - fs['down_len_nt']
            fs['end'] = fs['fs_coord'] + fs['up_len_nt']
        
        seq = Seq(gtdb, fs['seq_id'], read_seq=True)
        fsgene_seq = seq[fs['start']:fs['end']]
        if fs['strand'] == -1:
            fsgene_seq = fsgene_seq.reverse_complement()
        
        if fsgene_seq.upper() != fs['init_gene_seq'].upper():
            logging.warning("GT_FS and fsgene sequence do not match:\n%s\n%s" % (fs['init_gene_seq'], fsgene_seq))
            return None
        
        return fs, seq
    
    @staticmethod
    def create_new_in_db_from_SeqFeature(gtdb, f, seq, user_id, source=None):
        error = FSGene._validate_SeqFeature(f)
        if error is not None:
            logging.warning(error)
            return None
        
        # Determine fs_coord and fs_type
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
        
        # Validate fs_type and fs_coord, if possible
        if 'translation' in f.qualifiers:
            seq_dict = FSGene._get_prot_seq_parts(start1, end2, f.strand, fs_coord, fs_type, seq)
            my_prot = seq_dict['prot_seq_n'].lower() + seq_dict['prot_seq_c'].upper()
            true_prot = f.qualifiers['translation'][0]
            if my_prot.upper() != true_prot.upper():
                fmt = "\nThe true and reconstructed fs-prot seqs do not match:\n%s\n%s\n"
                raise Exception(fmt % (true_prot, my_prot))
        
        return db_id
    
    @staticmethod
    def _validate_SeqFeature(f):
        if not isinstance(f.location, CompoundLocation):
            return "Feature location is not a CompoundLocation object: '{}'".format(str(f.location))
        if len(f.location.parts) != 2:
            return "CompoundLocation '{}' has '{}' parts!".format(str(f.location), len(f.location.parts))
        if f.location.operator != 'join':
            return "Wrong feature location operator = '{}'".format(f.location.operator)
        if f.location.strand is None:
            return "The parts of CompoundLocation have different strands = '{}'".format(f.location)
        return None


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

def _get_genetic_code_from_SeqRecord(record):
    all_transl_f = list(filter(lambda f: 'transl_table' in f.qualifiers, record.features))
    if len(all_transl_f) == 0:
        logging.warning("Sequence '%s' doesn't have transl_table annotations", record.name)
        return None
    
    all_transl_tables = [f.qualifiers['transl_table'][0] for f in all_transl_f]
    transl_counts = Counter(all_transl_tables).most_common()  # https://stackoverflow.com/a/10797913/310453
    if len(transl_counts) > 1:
        logging.warning("Different genetic codes are annotated for seq: '%s'" % transl_counts)
    
    return transl_counts[0][0]

def db_xref_list_to_dict(dbxrefs):
    dbx_ll = [str.split(':') for str in dbxrefs]   # list of lists
    return {l[0] : l[1] for l in dbx_ll}

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


