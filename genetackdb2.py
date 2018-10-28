#!/usr/bin/env python3

import logging
import sys
import os
import re
import shutil
import datetime
from pprint import pprint
from collections import Counter

import Bio.Alphabet
import Bio.Seq
from Bio import SeqIO
#from Bio.Seq import Seq as BioSeq
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from mylib.baseutilbio import get_overlapping_feats_from_record
from mylib.baseutildb import BaseUtilDB
from mylib.local import PATH2


class GeneTackDB(BaseUtilDB):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.gtdb_dir = PATH2['gtdb2']
    
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
    def __init__(self, gtdb, db_id, main_sql, unit, add_prm=False, **kwargs):
        self.gtdb = gtdb
        self._unit = unit
        
        res = gtdb.exec_sql_ar(main_sql, db_id)
        if len(res) == 0:
            raise Exception("Unknown db_id = '{}'".format(db_id))
        
        # use dictionary keys as object arguments:
        # https://stackoverflow.com/a/2466207/310453
        for k, v in res[0].items():
            setattr(self, k, v)
        
        # Load params if requested
        self.prm = {}
        if add_prm:
            self.load_prm(**kwargs)
    
    def load_prm(self, prm_str=[], prm_int=[], prm_float=[], prm_list=[]):
        prm_sql = 'select name, value, num from {0}_params where {0}_id=%s'.format(
            self._unit)
        prm_res = self.gtdb.exec_sql_ar(prm_sql, self.id)
        for d in prm_res:
            if d['name'] in prm_list:
                if d['name'] not in self.prm:
                    self.prm[d['name']] = []
                self.prm[d['name']].append(d)
            elif d['name'] in prm_str:
                self.prm[d['name']] = d['value']
            elif d['name'] in prm_int:
                self.prm[d['name']] = int(d['num'])
            elif d['name'] in prm_float:
                self.prm[d['name']] = float(d['num'])
    
    def set_param(self, name, value=None, num=None):
        self.delete_param(name)
        self.add_param(name, value, num)
    
    def add_param(self, name, value=None, num=None):
        sql = 'INSERT INTO {0}_params ({0}_id, name, value, num)'.format(self._unit)
        self.gtdb.exec_sql_in(sql, self.id, name, value, num)
    
    def delete_param(self, name):
        sql = 'delete from {0}_params where {0}_id=%s and name=%s'.format(self._unit)
        self.gtdb.exec_sql_nr(sql, self.id, name)
    
    def get_param_with_full_path(self, name):
        return os.path.join(self.gtdb.gtdb_dir, self.prm[name])
    
    def get_myself_and_all_parents(self):
        '''Returns a dictionary with objects
        '''
        all_objs = {self._unit: self}
        if 'fsgene' in all_objs.keys():
            all_objs['seq'] = Seq(self.gtdb, all_objs['fsgene'].seq_id)
        elif 'sfeat' in all_objs.keys():
            all_objs['seq'] = Seq(self.gtdb, all_objs['sfeat'].seq_id)
        
        if 'seq' in all_objs.keys():
            all_objs['org'] = Org(self.gtdb, all_objs['seq'].org_id)
        
        return all_objs
    
    @staticmethod
    def object_for_unit(gtdb, unit, db_id):
        unit = unit.lower()
        if unit == 'org':
            return Org(gtdb, db_id)
        elif unit == 'seq':
            return Seq(gtdb, db_id)
        elif unit == 'sfeat':
            return SFeat(gtdb, db_id)
        elif unit == 'fsgene':
            return FSGene(gtdb, db_id)
        else:
            raise Exception('Unknown unit: %s'  % unit)

class Org(TableObject):
    def __init__(self, gtdb, db_id):
        TableObject.__init__(
            self, gtdb, db_id,
            main_sql = '''SELECT id, c_date, name, genus, phylum, kingdom, dir_path
                FROM orgs WHERE id=%s''',
            unit = 'org',
            add_prm = True,
            prm_str = ['BioSample', 'BioProject', 'Assembly', 'source_fn', 'short_name',
                      'blastdb_genome'],
            prm_list = ['taxonomy'])
        
        # Modify prm['taxonomy']: list of dicts  =>  list of strings
        if 'taxonomy' in self.prm:
            self.prm['taxonomy'] = [d['value'] for d in sorted(
                self.prm['taxonomy'], key=lambda d: d['num'])]
    
    def make_all_params(self):
        short_name = self.get_short_name()
        if short_name is None:
            logging.warning("Can't generate short name from '%s'" % self.name)
        else:
            self.set_param('short_name', value = short_name)
    
    def get_short_name(self):
        '''e.g. 'Mycobacterium tuberculosis H37Rv'  => 'M.tuberculosis'
        '''
        if re.compile(r'^(\w)\w*\s+(\w+).*$').match(self.name):
            return re.compile(r'^(\w)\w*\s+(\w+).*$').sub(r'\1.\2', self.name)
        else:
            return None
    
    def get_all_seq_ids(self):
        return [d['id'] for d in self.gtdb.exec_sql_ar(
            'select id from seqs where org_id=%s', self.id)]
    
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
        gtdb.exec_sql_in(
            "INSERT INTO orgs (id, user_id, name, genus, phylum, kingdom, dir_path)",
            db_id, user_id, name, genus, phylum, kingdom, dir_path)
        org = Org(gtdb, db_id)
        
        # ORG_PARAMS
        dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        for (name,value) in dbx_dict.items():
            org.add_param(name, value)
        
        num = 0
        for value in record.annotations.get('taxonomy',[]):
            org.add_param('taxonomy', value, num)
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


class BaseSeq(TableObject):
    def __init__(self, gtdb, db_id):
        """Basic Seq corresponding to the SEQ GeneTack DB table
        
        Arguments:
         - gtdb - [reqiured] GeneTackDB object
         - db_id - [required] valid ID from the SEQ table
        """
        TableObject.__init__(
            self, gtdb, db_id,
            main_sql = '''SELECT id, user_id, c_date, org_id, name, descr, type, ext_id, len
                FROM seqs WHERE id=%s''',
            unit = 'seq',
            add_prm = True,
            prm_str = ['fna_path', 'gbk_path'],
            prm_int = ['num_n', 'transl_table'],
            prm_float = ['gc'],
        )
    
    def get_all_fsgene_ids(self):
        return [d['id'] for d in self.gtdb.exec_sql_ar(
            'select id from fsgenes where seq_id=%s', self.id)]
    
    def read_genbank_file(self):
        gb_fn = os.path.join(self.gtdb.gtdb_dir, self.prm['gbk_path'])
        return SeqIO.read(gb_fn, "genbank")
    
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
        gtdb.exec_sql_in(
            "INSERT INTO seqs (id, user_id, org_id, name, descr, type, ext_id, len)",
            db_id, user_id, org.id, name, descr, m_type, ext_id, length)
        seq = BaseSeq(gtdb, db_id)
        
        # SEQ_PARAMS
        fna_path = Seq._save_fasta(gtdb, record, org, db_id)
        seq.add_param('fna_path', value=fna_path)
        
        gbk_path = Seq._save_gbk(gtdb, record, org)
        seq.add_param('gbk_path', value=gbk_path)
        
        # Get translation table
        transl_table = _get_genetic_code_from_SeqRecord(record)
        if transl_table is not None:
            seq.add_param('transl_table', num=transl_table)
        
        num_n = len(re.compile('[^ACGTacgt]').findall(str(record.seq)))
        seq.add_param('num_n', num=num_n)
        
        # GC content without ambiguous letters
        gc = GC(re.compile('[^ACGTacgt]').sub('', str(record.seq)))
        gc_str = '%.1f%%' % gc   # 72.1186611098208   => ' 72.1%'
        seq.add_param('gc', value=gc_str, num=gc)
        
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


class Seq(BaseSeq, Bio.Seq.UnknownSeq):
    def __init__(self, gtdb, db_id):
        BaseSeq.__init__(self, gtdb, db_id)
        Bio.Seq.UnknownSeq.__init__(self, self.len)
        self._has_seq = False


class SeqWithSeq(BaseSeq, Bio.Seq.Seq):
    def __init__(self, gtdb, db_id, use_seq=None):
        """Create BaseSeq object and read the actual sequence from file
        
        Arguments:
         - use_seq - [optional] take sequence from this Bio.Seq.Seq object
        """
        BaseSeq.__init__(self, gtdb, db_id)
        
        if use_seq is None:
            use_seq = self._read_seq_from_file()
        
        # make sure that len in DB and the actual sequence length match
        if self.len != len(use_seq):
            raise Exception("Sequence length ('%d') != value in DB ('%d')",
                            (self.len, len(use_seq)))
        
        Bio.Seq.Seq.__init__(self, str(use_seq), use_seq.alphabet)
        self._has_seq = True
    
    def _read_seq_from_file(self):
        fasta_fn = os.path.join(self.gtdb.gtdb_dir, self.prm['fna_path'])
        record = SeqIO.read(fasta_fn, "fasta")
        alphabet = Bio.Alphabet.SingleLetterAlphabet
        if self.type == 'DNA' and self.prm['num_n'] == 0:
            record.seq.alphabet = Bio.Alphabet.IUPAC.unambiguous_dna
        elif self.type == 'DNA' and self.prm['num_n'] > 0:
            record.seq.alphabet = Bio.Alphabet.IUPAC.ambiguous_dna
        return record.seq


class SFeat(TableObject):
    def __init__(self, gtdb, db_id):
        super().__init__(
            gtdb, db_id,
            main_sql = '''SELECT id, seq_id, type, start, end, strand, name, descr, ext_id
                FROM sfeats WHERE id=%s''',
            unit = 'sfeat',
            add_prm = True,
            prm_str = ['protein_id', 'ribosomal_slippage', 'pseudo', 'experiment',
                       'translation', 'nt_seq'],
            prm_list = ['EC_number'])
    
    def delete_from_db(self):
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_sfeats WHERE sfeat_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM sfeat_params  WHERE sfeat_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM sfeats        WHERE       id=%s', self.id)
    
    @staticmethod
    def create_new_in_db(gtdb, user_id, seq_id, start, end, strand, f_type,
                         name=None, descr=None, ext_id=None):
        db_id = gtdb.get_random_db_id('sfeats', 'id')
        gtdb.exec_sql_in("""INSERT INTO sfeats (
               id, user_id, seq_id,   type, start, end, strand, name, descr, ext_id)
        """,db_id, user_id, seq_id, f_type, start, end, strand, name, descr, ext_id
        )
        
        logging.info("New SFEAT has beed created with id '%d'" % db_id)
        return SFeat(gtdb, db_id)
    
    @staticmethod
    def get_db_id_by_SeqFeature(gtdb, f):
        """Checks if this feature has already been loaded to GTDB.
        """
        start = int(f.location.start)
        end = int(f.location.end)
        strand = f.location.strand
        res_l = gtdb.exec_sql_ar(
            'SELECT id FROM sfeats WHERE start=%s and end=%s and strand=%s',
            start, end, strand)
        if len(res_l) > 0:
            return res_l[0]['id']
        
        if 'locus_tag' in f.qualifiers:
            res_l = gtdb.exec_sql_ar('SELECT id FROM sfeats WHERE ext_id=%s',
                                     f.qualifiers['locus_tag'])
            if len(res_l) > 0:
                return res_l[0]['id']
        return None
    
    @staticmethod
    def create_new_in_db_from_SeqFeature(gtdb, user_id, seq_id, f):
        sfeat_id = SFeat.get_db_id_by_SeqFeature(gtdb, f)
        if sfeat_id is not None:
            logging.info('SeqFeature %s already present in DB (%d)' %
                         (f.location, sfeat_id))
            return SFeat(gtdb, sfeat_id)
        
        start = int(f.location.start)
        end = int(f.location.end)
        strand = f.location.strand
        return SFeat.create_new_in_db(
            gtdb, user_id, seq_id, start, end, strand, f.type,
            name = f.qualifiers.get('gene', [None])[0],
            descr = f.qualifiers.get('product', [None])[0],
            ext_id = f.qualifiers.get('locus_tag', [None])[0])
    
    @staticmethod
    def create_new_in_db_from_CDS_SeqFeature(gtdb, user_id, seq_id, f):
        sfeat_id = SFeat.get_db_id_by_SeqFeature(gtdb, f)
        if sfeat_id is not None:
            logging.info('SeqFeature %s already present in DB (%d)' %
                         (f.location, sfeat_id))
            return SFeat(gtdb, sfeat_id)
        
        sfeat = SFeat.create_new_in_db_from_SeqFeature(gtdb, user_id, seq_id, f)
        
        # SFEAT_PARAMS
        prm_keys = [
            'protein_id', 'ribosomal_slippage', 'pseudo', 'experiment',
            'EC_number', 'gene_synonym', 'function', 'regulatory_class', 'ncRNA_class']
        for k, vals in f.qualifiers.items():
            if k in prm_keys:
                for v in vals:
                    sfeat.add_param(k, v)
        
        for prot_seq in f.qualifiers.get('translation', []):
            sfeat.add_param('translation', value=prot_seq, num=len(prot_seq))
        
        dbx_dict = db_xref_list_to_dict(f.qualifiers.get('db_xref', []))
        for (name, value) in dbx_dict.items():
            sfeat.add_param(name, value)
        
        return sfeat
    
    @staticmethod
    def create_new_in_db_from_gbk_for_region(gtdb, user_id, seq_id, start, end, strand,
                                             record=None, min_overlap=1, max_sfeats=-1):
        '''Load features from GenBank file that overlap with the given target region.
        Arguments:
            - record - Biopython GenBank record object
            - max_sfeats - maximum number of sfeats with longest overlaps to create
        '''
        if record is None:
            record = Seq(gtdb, seq_id).read_genbank_file()
        
        overlapping_feats = get_overlapping_feats_from_record(
            record, start, end, strand, all_types=['CDS'],
            min_overlap=min_overlap, max_feats=max_sfeats)
        
        all_sfeats = []
        for f in overlapping_feats:
            all_sfeats.append(
                SFeat.create_new_in_db_from_CDS_SeqFeature(gtdb, user_id, seq_id, f))
        
        return all_sfeats

class FSGene(TableObject):
    def __init__(self, gtdb, db_id):
        self._prm_seq_names = ['prot_seq', 'prot_seq_n', 'prot_seq_c',
                               'nt_seq_corr', 'nt_seq_n', 'nt_seq_c']
        super().__init__(gtdb, db_id,
            main_sql = '''SELECT id, user_id, c_date, seq_id, cof_id, name, descr,
                fs_coord, fs_type, start, end, strand, source
                FROM fsgenes WHERE id=%s''',
            unit = 'fsgene',
            add_prm = True,
            prm_str = self._prm_seq_names,
            prm_int = [],
            prm_float = [])

    def delete_from_db(self):
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_sfeats WHERE fsgene_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM fsgene_params WHERE fsgene_id=%s', self.id)
        self.gtdb.exec_sql_nr('DELETE FROM fsgenes       WHERE        id=%s', self.id)
    
    def set_sfeats(self, all_sfeats):
        '''make DESCR like "magnesium chelatase | VWA domain-containing protein"
        '''
        # Sort in the order of translation
        all_sfeats = sorted(all_sfeats, key=lambda s: s.start)
        if self.strand < 0:
            all_sfeats.reverse()
        
        # Remove old sfeat links
        self.gtdb.exec_sql_nr('delete from fsgene_sfeats where fsgene_id=%s', self.id)
        
        sfeat_descr_l = []
        for sfeat in all_sfeats:
            self.gtdb.exec_sql_in('INSERT INTO fsgene_sfeats (fsgene_id, sfeat_id)', self.id, sfeat.id)
            if sfeat.descr is not None:
                sfeat_descr_l.append(sfeat.descr)
        
        # Update fsgene descr
        fsgene_descr = None
        if len(sfeat_descr_l) > 0:
            fsgene_descr = ' | '.join(sfeat_descr_l)
        self.gtdb.exec_sql_nr('update fsgenes set descr=%s where id=%s', fsgene_descr, self.id)
    
    def make_all_params(self, seq=None):
        if seq is None:
            seq = SeqWithSeq(self.gtdb, self.seq_id)
        
        # make NAME like: 'NC_010002.1:1922457:-1'
        if self.fs_type == 0:
            # to avoid 'NC_005085.1:1684916:+0'
            name = '%s:%d:%d'  % (seq.ext_id, self.fs_coord, self.fs_type)
        else:
            name = '%s:%d:%+d' % (seq.ext_id, self.fs_coord, self.fs_type)
        self.gtdb.exec_sql_nr('update fsgenes set name=%s where id=%s', name, self.id)
        
        try:
            self.make_prm_seqs(seq)
        except Exception as e:
            logging.error("Couldn't generate prm_seqs for fsgene '%d':\n%s" % (self.id, e))
        
    def make_prm_seqs(self, seq):
        seq_dict = FSGene._get_prot_seq_parts(self.start, self.end, self.strand,
                                              self.fs_coord, self.fs_type, seq)
        for name in self._prm_seq_names:
            self.set_param(name, value=seq_dict[name], num=len(seq_dict[name]))
    
    def _get_prot_seq_parts(start, end, strand, fs_coord, fs_type, seq):
        up_len_nt = fs_coord - start
        down_len_nt = end - fs_coord
        if strand == -1:
            up_len_nt, down_len_nt = down_len_nt, up_len_nt
        
        # Backward frameshiting increases the length of the C-terminal part of the fs-protein
        down_len_nt -= fs_type
        
        if up_len_nt % 3 != 0 or down_len_nt % 3 != 0:
            fmt = """
            The length of the upstream or downstream fsgene part (%d / %d) is not divisible by 3:
            start=%d, end=%d, fs_coord=%d, fs_type=%+d, strand=%d
            """
            raise Exception(fmt % (up_len_nt, down_len_nt, start, end, fs_coord, fs_type, strand))
        
        chunk_nt = seq[start:end]
        if strand == -1:
            chunk_nt = chunk_nt.reverse_complement()
        
        up_chunk_nt = chunk_nt[:up_len_nt]
        down_chunk_nt = chunk_nt[-down_len_nt:]
        prot_seq_n = up_chunk_nt.translate(table = seq.prm['transl_table'])
        prot_seq_c = down_chunk_nt.translate(table = seq.prm['transl_table'])
        prot_seq_c = prot_seq_c.rstrip('*')  # remove possible stop codon at the end
        
        if '*' in prot_seq_n or '*' in prot_seq_c:
            raise Exception("FS-prot seq contains in-frame stop codon:\n%s\n\n%s" % (prot_seq_n, prot_seq_c))
        
        return {
            'nt_seq_n': up_chunk_nt.upper(), 'nt_seq_c': down_chunk_nt.upper(),
            'nt_seq_corr': up_chunk_nt.lower() + down_chunk_nt.upper(),
            'prot_seq_n': prot_seq_n.upper(), 'prot_seq_c': prot_seq_c.upper(),
            'prot_seq': prot_seq_n.lower() + prot_seq_c.upper(),
        }
    
    @staticmethod
    def create_new_in_db(gtdb, seq_id, user_id, start, end, strand, fs_coord, fs_type,
                         source=None, cof_id=None, c_date=None, db_id=None):
        # check if fsgene already exists in this location
        existing_fsgenes = FSGene.get_ids_for_seq_coords(
            gtdb, seq_id, start, end)
        if len(existing_fsgenes) > 0:
            logging.warning("FSGene '%d' already exists for coords %d-%d" %
                            (existing_fsgenes[0], start, end))
            return None
        
        if db_id is None:
            db_id = gtdb.get_random_db_id('fsgenes', 'id')
        if c_date is None:
            c_date = datetime.datetime.now()
        
        gtdb.exec_sql_in("""
            INSERT INTO fsgenes (
            id, c_date, user_id, seq_id, fs_coord, fs_type,
            start, end, strand, source, cof_id)
            """, db_id, c_date, user_id, seq_id, fs_coord, fs_type,
            start, end, strand, source, cof_id)
        
        logging.info("New FSGENE has beed created with id '%d'" % db_id)
        return FSGene(gtdb, db_id)
    
    @staticmethod
    def create_new_in_db_from_Seq_Feature(gtdb, f, seq, user_id, source=None):
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
    
    @staticmethod
    def get_db_id_by_sfeat_id(gtdb, sfeat_id):
        res_l = gtdb.exec_sql_ar(
            'SELECT fsgene_id FROM fsgene_sfeats WHERE sfeat_id=%s', sfeat_id)
        
        if len(res_l) == 0:
            return None
        elif len(res_l) > 1:
            raise Exception("SFEAT_ID '%d' matches several FSGenes!" % sfeat_id)
        else:
            return res_l[0]['fsgene_id']
    
    @staticmethod
    def get_ids_for_seq_coords(gtdb, seq_id, left, right):
        return [d['id'] for d in gtdb.exec_sql_ar(
            """select distinct id from fsgenes
            where seq_id={S} and (
                ({L} <= start AND start <= {R}) OR
                ({L} <= end AND end <= {R})
            )""".format(S=seq_id, L=left, R=right))]

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
        logging.warning("Can't determine genus for species name '%s' and taxonomy %s" %
                        (species, taxa_l))
        genus = None
    
    if genus is not None and ' ' in genus:
        # Genus must be a single word!
        logging.warning("Wrong genus '%s' for species name '%s'" % (genus, species))
        genus = None
    
    phylum = taxa_l[1]
    return genus, phylum, kingdom


