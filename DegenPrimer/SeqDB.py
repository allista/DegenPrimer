# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
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
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from StringIO import StringIO


class SeqDB(object):
    '''Create and manage a database of sequences with fast approximate match 
    searching.'''
    
    _ambiguous = set(IUPAC.ambiguous_dna.letters.upper()) - set(IUPAC.unambiguous_dna.letters.upper())
    
    @classmethod
    def is_ambiguous(cls, seqrec):
        sstr = str(seqrec.seq).upper()
        for l in cls._ambiguous:
            if l in sstr: return True
        return False
    
    _sequences_schema = \
    '''CREATE TABLE "sequences" (
           "id"         INTEGER PRIMARY KEY AUTOINCREMENT,
           "name"       TEXT NOT NULL,
           "sequence"   TEXT NOT NULL,
           "length"     INTEGER UNSIGNED NOT NULL)'''
  
    
    def __init__(self):
        self._db_name  = None
        self._db       = None
        self._cursor   = None
    #end def
    
    def __del__(self): self.close()
    
    
    @property
    def name(self): return self._db_name
    
    def __nonzero__(self):
        return not self._cursor is None 
        
    
    def _configure_connection(self):
        self._cursor.execute('PRAGMA cache_size=64000')
        self._cursor.execute('PRAGMA synchronous=OFF')
        self._cursor.execute('PRAGMA count_changes=OFF')
        self._cursor.execute('PRAGMA temp_store=MEMORY')
    #end def
    
        
    def _init_db(self, filename):
        #if db exists, back it up
        if os.path.isfile(filename):
            shutil.move(filename, filename+'.back')
        #create new db
        self._db = sqlite3.connect(filename, isolation_level='DEFERRED')
        self._cursor = self._db.cursor()
        self._configure_connection()
        #init tables
        self._cursor.execute(self._sequences_schema)
        self._db.commit()
        self._db_name = filename
    #end def
    
    
    def _populate_db(self, sequences):
        self._cursor.execute('PRAGMA cache_size=500000')
        #populate database with data
        for sequence in sequences:
#            if self.is_ambiguous(sequence):
#                print '\nSequence %s is excluded because it contains ambiguous symbols\n' % sequence.id
#                continue
            self._cursor.execute('''
            INSERT INTO sequences (name, sequence, length)
            VALUES (?, ?, ?)''', 
            (sequence.description or sequence.name or sequence.id, 
             sequence.format('fasta'), len(sequence)))
        self._db.commit()
        self._cursor.execute('PRAGMA cache_size=64000')
        self._db.commit()
    #end def
    
    
    def create_db(self, filename, sequences):
        '''Create and SQLite database of sequences in a given file.'''
        self._init_db(filename)
        self._populate_db(sequences)
    #end def
    
    
    def _get_sequences(self, files):
        if not files: return None
        sequences = []
        for filename in files:
            try: sequences.extend(SeqIO.parse(filename, 'fasta', IUPAC.unambiguous_dna))
            except Exception, e:
                print '\nError while parsing %s' % filename
                print e
                print
        if not sequences:
            print '\nNo sequences were found in given fasta files.' 
            return None
        return sequences
    #end def
    
    
    def create_db_from_files(self, db_filename, files):
        sequences = self._get_sequences(files)
        if sequences:
            self.create_db(db_filename, sequences)
            return True
        return False
    #end def
    
    
    def add_files_to_db(self, files):
        if not files or not self._db: return False
        sequences = self._get_sequences(files)
        if sequences: 
            self._populate_db(sequences)
            return True
        return False
    #end def
    
    
    def del_from_db(self, ids):
        if self._db is None or not ids: return False
        sql_string = '''
        DELETE FROM sequences
        WHERE id IN (%s)
        ''' % ','.join(['?']*len(ids))
        self._cursor.execute(sql_string, ids)
        self._db.commit()
        return True
    #end def
        
    
    def connect(self, filename):
        '''Connect to an existent file database of sequences.'''
        if self._db is not None:
            print '\nSeqDB.connect: already connected to a database. Call close() first.'
            return False
        if filename != ':memory:' and not os.path.isfile(filename):
            print '\nSeqDB.connect: no such file:\n%s' % filename 
            return False
        self._db = sqlite3.connect(filename, isolation_level='DEFERRED')
        self._cursor = self._db.cursor()
        self._configure_connection()
        self._db_name = filename
        return True
    #end def
    
    
    def close(self):
        if self._db is None: return
        self._db.commit()
        self._db.close()
        self._db      = None
        self._cursor  = None
        self._db_name = None
    #end def
    
    
    def get_names(self, seq_ids=None):
        if self._db is None: return None
        if seq_ids:
            sql_string = '''
            SELECT id, name from sequences
            WHERE id IN (%s)
            ''' % ','.join(['?']*len(seq_ids))
            self._cursor.execute(sql_string, seq_ids)
        else:
            self._cursor.execute('''SELECT id, name from sequences''')
        return dict(self._cursor)
    #end def
    
    
    def get_lengths(self, seq_ids=None):
        if self._db is None: return None
        if seq_ids:
            sql_string = '''
            SELECT id, length from sequences
            WHERE id IN (%s)
            ''' % ','.join(['?']*len(seq_ids))
            self._cursor.execute(sql_string, seq_ids)
        else:
            self._cursor.execute('''SELECT id, length from sequences''')
        return dict(self._cursor)
    #end def
    
    
    def get_seqs(self, seq_ids=None):
        if self._db is None: return None
        if seq_ids:
            sql_string = '''
            SELECT id, length, sequence from sequences
            WHERE id IN (%s)
            ''' % ','.join(['?']*len(seq_ids))
            self._cursor.execute(sql_string, seq_ids)
        else:
            self._cursor.execute('''SELECT id, length, sequence from sequences''')
        return [(sid, length, SeqIO.read(StringIO(seq), 'fasta')) 
                for sid, length, seq in self._cursor]
    #end def
#end class