# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Dec 21, 2012

@author: Allis Tauri <allista@gmail.com>
'''
import os
import shutil
import sqlite3
import numpy as np
import multiprocessing
from time import time
from datetime import timedelta
from array import array
from StringIO import StringIO
from scipy.fftpack import fft, ifft
from SecStructures import Duplex

class SeqDB(object):
    '''Create and manage a database of sequences with fast approximate match 
    searching.'''
    
    _sequences_schema = \
    '''CREATE TABLE "sequences" (
           "id"         INTEGER PRIMARY KEY AUTOINCREMENT,
           "name"       TEXT NOT NULL,
           "sequence"   TEXT NOT NULL,
           "fwd_AT_map" BLOB NOT NULL,
           "rev_AT_map" BLOB NOT NULL,
           "fwd_GC_map" BLOB NOT NULL,
           "rev_GC_map" BLOB NOT NULL)'''
    
#    _strands_schema = \
#    '''CREATE TABLE "strands" (
#           "id"        INTEGER PRIMARY KEY AUTOINCREMENT,
#           "sequence"  BLOB NOT NULL,
#           "strand"    INTEGER NOT NULL,
#           "AT_map_id" INTEGER NOT NULL,
#           "GC_map_id" INTEGER NOT NULL)'''
#    
#    _maps_schema = \
#    '''CREATE TABLE "maps" (
#           "id" INTEGER PRIMARY KEY AUTOINCREMENT,
#           "letters"   TEXT NOT NULL,
#           "map"       BLOB NOT NULL)'''
    
    
    _w_3_0 = 1
    _w_3_1 = (-1/2.0+np.sqrt(3)/2.0j)
    _w_3_2 = (-1/2.0-np.sqrt(3)/2.0j)
    
    _unambiguous = array('b','ATGC')
    _ambiguous   = _unambiguous+array('b', 'RYSWKMBDHVN') 
    
    _T_AT_mapping = dict(zip(_unambiguous, (_w_3_1,_w_3_2,0,0))) 
#                    {'A':_w_3_1,
#                     'T':_w_3_2,
#                     'G':0,'C':0}
    
    _T_GC_mapping = dict(zip(_unambiguous, (0,0,_w_3_1,_w_3_2)))
#                    {'A':0,'T':0,
#                     'G':_w_3_1,
#                     'C':_w_3_2}
    
#    _P_AT_mapping = dict(zip(_ambiguous, 
#                            (_w_3_2,
#                            _w_3_1,
#                            _w_3_0,
#                            _w_3_0,
#                            _w_3_2,
#                            _w_3_1,
#                            _w_3_0,
#                            _w_3_2+_w_3_1,
#                            _w_3_1,
#                            _w_3_2,
#                            _w_3_1,
#                            _w_3_2+_w_3_1,
#                            _w_3_2+_w_3_1,
#                            _w_3_2,
#                            _w_3_2+_w_3_1)))
    _P_AT_mapping = {'A':_w_3_2,
                     'T':_w_3_1,
                     'G':_w_3_0,
                     'C':_w_3_0,
                     
                     'R':_w_3_2,
                     'Y':_w_3_1,
                     'S':_w_3_0,
                     'W':_w_3_2+_w_3_1,
                     'K':_w_3_1,
                     'M':_w_3_2,
                     'B':_w_3_1,
                     'D':_w_3_2+_w_3_1,
                     'H':_w_3_2+_w_3_1,
                     'V':_w_3_2,
                     'N':_w_3_2+_w_3_1}
    
#    _P_GC_mapping = dict(zip(_ambiguous, 
#                            (_w_3_0,
#                             _w_3_0,
#                             _w_3_2,
#                             _w_3_1,
#                             _w_3_2,
#                             _w_3_1,
#                             _w_3_2+_w_3_1,
#                             _w_3_0,
#                             _w_3_2,
#                             _w_3_1,
#                             _w_3_2+_w_3_1,
#                             _w_3_2,
#                             _w_3_1,
#                             _w_3_2+_w_3_1,
#                             _w_3_2+_w_3_1)))
    _P_GC_mapping = {'A':_w_3_0,
                     'T':_w_3_0,
                     'G':_w_3_2,
                     'C':_w_3_1,
                     
                     'R':_w_3_2,
                     'Y':_w_3_1,
                     'S':_w_3_2+_w_3_1,
                     'W':_w_3_0,
                     'K':_w_3_2,
                     'M':_w_3_1,
                     'B':_w_3_2+_w_3_1,
                     'D':_w_3_2,
                     'H':_w_3_1,
                     'V':_w_3_2+_w_3_1,
                     'N':_w_3_2+_w_3_1}
    
    _chunk_size = 100000

    def __init__(self):
        self._db = None
        self._cursor = None
    #end def
    
    
    @classmethod
    def _map_letter(cls, letter, _map):
        try: return _map[letter]
        except KeyError: return 0
    #end def
    
    
    @classmethod
    def _map_template(cls, template, map_len):
        fwd_seq = template.seq
        rev_seq = template.seq.reverse_complement()
        fwd_AT_map = np.array(array('b',str(fwd_seq)), dtype=complex)
        fwd_AT_map.resize(map_len)
        rev_AT_map = np.array(array('b',str(rev_seq)), dtype=complex)
        rev_AT_map.resize(map_len)
        fwd_GC_map = fwd_AT_map.copy()
        rev_GC_map = rev_AT_map.copy()
        for k,v in cls._T_AT_mapping.iteritems():
            fwd_AT_map[fwd_AT_map == k] = v
            rev_AT_map[rev_AT_map == k] = v
        for k,v in cls._T_GC_mapping.iteritems():
            fwd_GC_map[fwd_GC_map == k] = v
            rev_GC_map[rev_GC_map == k] = v
        return (fwd_seq, rev_seq, 
                (fwd_AT_map[::-1], fwd_GC_map[::-1]),
                (rev_AT_map[::-1], rev_GC_map[::-1]))
    #end def
    
    
    @classmethod
    def _map_pattern(cls, pattern, map_len):
        #this naive algorithm works ~5 times faster 
        #than the one used for template mapping due to short length of patterns
        AT_map = np.zeros(map_len, dtype=complex)
        GC_map = np.zeros(map_len, dtype=complex)
        for i,letter in enumerate(pattern):
            AT_map[i] = cls._map_letter(letter, cls._P_AT_mapping)
            GC_map[i] = cls._map_letter(letter, cls._P_GC_mapping)
        return (AT_map, GC_map)
    #end def

    @classmethod
    def _compile_duplexes(cls, template, primer, matches):
        p_len = len(primer)
        results = []
        for i in matches:
            duplexes = []
            print i #TODO:remove
            for var in primer.sequences:
                duplexes.append(Duplex(var, template[i:i+p_len].reverse_complement()))
                print duplexes[-1] #TODO:remove
            results.append((i+p_len-1, duplexes))
        return results
    #end def
    
    
    @classmethod
    def _find_in_chunk(cls, t_chunk, p_fft, c_stride, map_len):
        t_AT_map = np.array(array('b',str(t_chunk)), dtype=complex)
        t_AT_map.resize(map_len)
        t_GC_map = t_AT_map.copy()
        for k,v in cls._T_AT_mapping.iteritems():
            t_AT_map[t_AT_map == k] = v
        for k,v in cls._T_GC_mapping.iteritems():
            t_GC_map[t_GC_map == k] = v
        AT_score = ifft(fft(t_AT_map[::-1])*p_fft[0])[::-1][:c_stride]
        GC_score = ifft(fft(t_GC_map[::-1])*p_fft[1])[::-1][:c_stride]
        score    = AT_score.real + GC_score.real
        return score
    #end def
    
    
    @classmethod 
    def find_by_chunks(cls, template, primer, mismatches):
        pattern    = primer.master_sequence.seq
        p_len      = len(pattern)
        t_len      = len(template)
        if p_len > t_len: 
            raise ValueError('find: template sequence should be longer or equal'
                             ' to primer sequence.')
        chunk_size   = min(t_len, (cls._chunk_size/p_len+1)*p_len)
        chunk_stride = chunk_size-p_len+1
        map_len    = int(2**(np.ceil(np.log2(chunk_size))))
        p_maps     = cls._map_pattern(pattern, map_len)
        p_AT_fft   = fft(p_maps[0])
        p_GC_fft   = fft(p_maps[1])
        fwd_seq    = template.seq
        rev_seq    = template.seq.reverse_complement()
        fwd_score  = np.ndarray(0,dtype=float)
        rev_score  = np.ndarray(0,dtype=float)
        correction = np.ndarray(chunk_stride)
        correction.fill(p_len/3.0)
        i = 0
        while i < t_len:
            front = min(t_len, i+chunk_size)
            score = cls._find_in_chunk(fwd_seq[i:front], (p_AT_fft,p_GC_fft), 
                                       chunk_stride, map_len)
            score = (score + correction - score/3.0)
            fwd_score = np.concatenate([fwd_score, score])
            score = cls._find_in_chunk(rev_seq[i:front], (p_AT_fft,p_GC_fft), 
                                       chunk_stride, map_len)
            score = (score + correction - score/3.0)
            rev_score = np.concatenate([rev_score, score])
            i += chunk_stride
        #match indeces
        matches     = max(1, p_len - mismatches)-0.5
        fwd_matches = np.arange(t_len-p_len+1)[fwd_score[:t_len-p_len+1] >= matches]
        rev_matches = np.arange(t_len-p_len+1)[rev_score[:t_len-p_len+1] >= matches]
        del fwd_score, rev_score
        #construct duplexes
        fwd_results = cls._compile_duplexes(fwd_seq, primer, fwd_matches)
        rev_results = cls._compile_duplexes(rev_seq, primer, rev_matches)
        return fwd_results,rev_results
    #end def
    
    
    @classmethod
    def find(cls, template, primer, mismatches):
        pattern = primer.master_sequence
        p_len   = len(pattern)
        t_len   = len(template)
        if p_len > t_len: 
            raise ValueError('find: template sequence should be longer or equal'
                             ' to primer sequence.')
        map_len = int(2**(np.ceil(np.log2(t_len))))
        matches = max(1, p_len - mismatches)-0.5
        fwd_seq,rev_seq,fwd_maps,rev_maps = cls._map_template(template, map_len)
        pattern_maps  = cls._map_pattern(pattern, map_len)
        #AT convolution
        P_AT_fft     = fft(pattern_maps[0])
        fwd_AT_score = ifft(fft(fwd_maps[0])*P_AT_fft)[::-1][:t_len-p_len+1]
        rev_AT_score = ifft(fft(rev_maps[0])*P_AT_fft)[::-1][:t_len-p_len+1]
        del P_AT_fft
        #GC convolution
        P_GC_fft     = fft(pattern_maps[1])
        fwd_GC_score = ifft(fft(fwd_maps[1])*P_GC_fft)[::-1][:t_len-p_len+1]
        rev_GC_score = ifft(fft(rev_maps[1])*P_GC_fft)[::-1][:t_len-p_len+1]
        del P_GC_fft
        del fwd_maps, rev_maps
        #total scores
        fwd_score    = fwd_AT_score.real + fwd_GC_score.real
        rev_score    = rev_AT_score.real + rev_GC_score.real
        del fwd_AT_score, fwd_GC_score, rev_AT_score, rev_GC_score
        #correct scores
        correction = np.ndarray(t_len-p_len+1)
        correction.fill(p_len/3.0)
        fwd_score  = (fwd_score + correction - fwd_score/3.0)
        rev_score  = (rev_score + correction - rev_score/3.0)
        del correction
        #match indeces
        fwd_matches = np.arange(t_len-p_len+1)[fwd_score >= matches]
        rev_matches = np.arange(t_len-p_len+1)[rev_score >= matches]
        del fwd_score, rev_score
        #construct duplexes
        fwd_results = cls._compile_duplexes(fwd_seq, primer, fwd_matches)
        rev_results = cls._compile_duplexes(rev_seq, primer, rev_matches)
        return fwd_results,rev_results
    #end def
    
    
    def _init_db(self, filename):
        #if db exists, back it up
        if os.path.isfile(filename):
            shutil.move(filename, filename+'.back')
        #create new db
        self._db = sqlite3.connect(filename, isolation_level='DEFERRED')
        self._cursor = self._db.cursor()
        #init tables
        self._cursor.execute(self._sequences_schema)
#        self._cursor.execute(self._strands_schema)
#        self._cursor.execute(self._maps_schema)
        self._db.commit()
        self.close()
    #end def
    
    
    def create_db(self, filename, sequences):
        #create database tables
        self._init_db(filename)
        #connect to the new database
        self.connect(filename)
        self._cursor.execute('PRAGMA cache_size=500000')
        #populate database with data
        for sequence in sequences:
            fwd_str = str(sequence.seq)
            rev_str = str(sequence.seq.reverse_complement())
            #map strands to the 3d roots of unity
            fwd_AT_map = np.ndarray(len(fwd_str), dtype=complex)
            rev_AT_map = np.ndarray(len(fwd_str), dtype=complex)
            fwd_GC_map = np.ndarray(len(fwd_str), dtype=complex)
            rev_GC_map = np.ndarray(len(fwd_str), dtype=complex)
            for i in range(len(fwd_str)):
                fl = fwd_str[i]
                rl = rev_str[i]
                fwd_AT_map[i] = self._map_letter(fl, self._T_AT_mapping)
                rev_AT_map[i] = self._map_letter(rl, self._T_AT_mapping)
                fwd_GC_map[i] = self._map_letter(fl, self._T_GC_mapping)
                rev_GC_map[i] = self._map_letter(rl, self._T_GC_mapping)
            #fill sequences table
            str_out = StringIO()
            np.save(str_out, fwd_AT_map)
            np.save('fwd_AT_map', fwd_AT_map)
            fwd_AT_map = str_out.getvalue()
            str_out.truncate(0)
            np.save(str_out, rev_AT_map)
            rev_AT_map = str_out.getvalue()
            str_out.truncate(0)
            np.save(str_out, fwd_GC_map)
            fwd_GC_map = str_out.getvalue()
            str_out.truncate(0)
            np.save(str_out, rev_GC_map)
            rev_GC_map = str_out.getvalue()
            str_out.close()
            self._cursor.execute('''
            INSERT INTO sequences (name, sequence, 
                                   fwd_AT_map, rev_AT_map,
                                   fwd_GC_map, rev_GC_map)
            VALUES (?, ?, ?, ?, ?, ?)''', (sequence.id, fwd_str, 
                                           sqlite3.Binary(fwd_AT_map), 
                                           sqlite3.Binary(rev_AT_map),
                                           sqlite3.Binary(fwd_GC_map), 
                                           sqlite3.Binary(rev_GC_map)))
#            #fill in tables in reverse order: maps table first
#            str_out = StringIO()
#            np.save(str_out, fwd_AT_map)
#            self._cursor.execute('''
#            INSERT INTO maps (letters, map)
#            VALUES (?, ?)''', ('AT',  sqlite3.Binary(str_out.getvalue())))
#            fwd_AT_map_id = self._cursor.lastrowid
#            str_out.truncate(0)
#            np.save(str_out, rev_AT_map)
#            self._cursor.execute('''
#            INSERT INTO maps (letters, map)
#            VALUES (?, ?)''', ('AT',  sqlite3.Binary(str_out.getvalue())))
#            rev_AT_map_id = self._cursor.lastrowid
#            str_out.truncate(0)
#            np.save(str_out, fwd_GC_map)
#            self._cursor.execute('''
#            INSERT INTO maps (letters, map)
#            VALUES (?, ?)''', ('GC',  sqlite3.Binary(str_out.getvalue())))
#            fwd_GC_map_id = self._cursor.lastrowid
#            str_out.truncate(0)
#            np.save(str_out, rev_GC_map)
#            self._cursor.execute('''
#            INSERT INTO maps (letters, map)
#            VALUES (?, ?)''', ('GC',  sqlite3.Binary(str_out.getvalue())))
#            rev_GC_map_id = self._cursor.lastrowid
#            str_out.close()
#            #strands table
#            self._cursor.execute('''
#            INSERT INTO strands (sequence, strand, AT_map_id, GC_map_id)
#            VALUES (?, ?, ?, ?)''', (fwd_str.tostring(),  1, fwd_AT_map_id, fwd_GC_map_id))
#            fwd_str_id = self._cursor.lastrowid
#            self._cursor.execute('''
#            INSERT INTO strands (sequence, strand, AT_map_id, GC_map_id)
#            VALUES (?, ?, ?, ?)''', (rev_str,  1, rev_AT_map_id, rev_GC_map_id))
#            rev_str_id = self._cursor.lastrowid
#            #sequences table
#            self._cursor.execute('''
#            INSERT INTO sequences (name, fwd_strand_id, rev_strand_id)
#            VALUES (?, ?, ?)''', (sequence.id, fwd_str_id, rev_str_id))
        self._db.commit()
        self.close()
    #end def
    
    
    def connect(self, filename):
        if not os.path.isfile(filename): return False
        self._db = sqlite3.connect(filename)
        self._cursor = self._db.cursor()
        self._cursor.execute('PRAGMA cache_size=64000')
        self._cursor.execute('PRAGMA synchronous=OFF')
        self._cursor.execute('PRAGMA count_changes=OFF')
        self._cursor.execute('PRAGMA temp_store=MEMORY')
        return True
    #end def
    
    def close(self):
        if not self._db is None: return
        self._db.close()
        self._db = None
        self._cursor = None 
#end class


#tests
if __name__ == '__main__':
    import cProfile
    import sys
    from StringTools import print_exception
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Primer import Primer
    
    os.chdir('../')
    try:
        record_file = open('Ch5_gnm.fa', 'r')
    except IOError, e:
        print 'Unable to open Ch5_gnm.fa'
        print_exception(e)
        sys.exit(1)
    template = SeqIO.read(record_file, 'fasta', IUPAC.unambiguous_dna)
    record_file.close()
    query = Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna)
    
#    'AAGCGTGCTTAGCTAGTCAGTGACGATG'
#    'GACG'
#    '121111101111110211110400'
#    template = SeqRecord(Seq('AAGCGTGCTTAGCTAGTCAGTGACGATG', 
#                             IUPAC.ambiguous_dna), id='test')
#    query = Seq('GACG', IUPAC.ambiguous_dna)

    primer = Primer(SeqRecord(query, id='test'), 0.1e-6)
    seq_db = SeqDB()
    cProfile.run("seq_db.find_by_chunks(template, primer, 3); seq_db.find(template, primer, 3)", 'test.profile')
    #cProfile.run("seq_db.create_db('test.db', (record, ))", 'test.profile')
