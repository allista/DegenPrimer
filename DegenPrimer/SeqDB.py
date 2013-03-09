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
import errno
import multiprocessing as mp
from Queue import Empty
from time import time, sleep
from StringIO import StringIO
from StringTools import print_exception
from SearchEngine import SearchEngine


class SeqDB(object):
    '''Create and manage a database of sequences with fast approximate match 
    searching.'''
    
    _sequences_schema = \
    '''CREATE TABLE "sequences" (
           "id"         INTEGER PRIMARY KEY AUTOINCREMENT,
           "name"       TEXT NOT NULL,
           "sequence"   TEXT NOT NULL,
           "length"     INTEGER UNSIGNED NOT NULL)'''
  
    
    def __init__(self):
        self._db       = None
        self._cursor   = None
        self._searcher = SearchEngine() 
    #end def
    
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
    #end def
    
    
    def _populate_db(self, sequences):
        self._cursor.execute('PRAGMA cache_size=500000')
        #populate database with data
        for sequence in sequences:
            self._cursor.execute('''
            INSERT INTO sequences (name, sequence, length)
            VALUES (?, ?, ?)''', 
            (sequence.description, sequence.format('fasta'), len(sequence)))
        self._db.commit()
        self._cursor.execute('PRAGMA cache_size=64000')
    #end def
    
    
    def create_db(self, filename, sequences):
        '''Create and SQLite database of sequences in a given file.'''
        self._init_db(filename)
        self._populate_db(sequences)
    #end def
    
    
    def create_db_from_files(self, db_filename, files):
        if not files: return False
        sequences = []
        for filename in files:
            if os.path.isfile(filename):
                sequences += SeqIO.parse(filename, 'fasta', IUPAC.unambiguous_dna)
            else: print '\nSeqDB: no such file %s' % filename
        if not sequences:
            print '\nNo sequences were found in given fasta files.' 
            return False
        self.create_db(db_filename, sequences)
        return True
    #end def
    
    
    def connect(self, filename):
        '''Connect to an existent file database of sequences.'''
        if filename != ':memory:' and not os.path.isfile(filename):
            print '\nNo such file: %s' % filename 
            return False
        self._db = sqlite3.connect(filename, isolation_level='DEFERRED')
        self._cursor = self._db.cursor()
        self._configure_connection()
        return True
    #end def
    
    
    def close(self):
        if not self._db is None: return
        self._db.close()
        self._db = None
        self._cursor = None
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
    
    
    def find_in_db(self, primer, mismatches, seq_ids=None):
        '''In the connected DB find all occurrences of the possibly degenerate 
        primer with number of mismatches less than or equal to the given number.
        If an optional list of sequence ids is provided, search will be performed 
        only within corresponding subset of sequences.
        Return a dict of following structure: 
        [sequence_id: [fwd_matches, rev_matches], ...]
        where *_matches are lists:
        [(3'-position, [(matching duplex, id), ...]), ...]
        If no DB is connected, return None.'''
        if self._db is None: return None
        seqs = self.get_seqs(seq_ids)
        if not seqs: 
            return None
        elif len(seqs) == 1:
            results = self._searcher.find(seqs[0][2], primer, mismatches)
            if results: return {seqs[0][0]: results}
            else: return None
        else:
            return self._searcher.batch_find(seqs, primer, mismatches)
    #end def
    
    
    def abort_search(self): self._searcher.abort()
#end class


#tests
import signal
import cProfile
import sys, csv, timeit
from functools import partial
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Primer import Primer

searcher = None
ppid     = -1
data1    = []
data2    = []

def sig_handler(signal, frame):
    if ppid != os.getpid():
        return
    print 'Aborting main process %d' % os.getpid()
    if searcher:
        searcher.abort()
    if data1:
        print 'Write out gathered data1...'
        out_file = open('gather_data1-%d.csv' % time(), 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data1)
        out_file.close()
        print 'Done.'
    if data2:
        print 'Write out gathered data2...'
        out_file = open('gather_data2-%d.csv' % time(), 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data2)
        out_file.close()
        print 'Done.'
    sys.exit(1)
#end def


#tests
if __name__ == '__main__':
    #setup signal handler
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)
    
    os.chdir('../')
#    try:
#        record_file = open('Ch5_gnm.fa', 'r')
#    except IOError, e:
#        print 'Unable to open Ch5_gnm.fa'
#        print_exception(e)
#        sys.exit(1)
#    record = SeqIO.read(record_file, 'fasta', IUPAC.unambiguous_dna)
#    record_file.close()
    
    sdb = SeqDB()
    sdb.create_db_from_files(':memory:', ('Ch5_gnm.fa',))
    print sdb.get_names()
    seqs = sdb.get_seqs([1,])
    if seqs:
        print seqs[0][0]
        print seqs[0][1]
    
    query = Primer(SeqRecord(Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna), id='test'), 0.1e-6)
    
    
    def print_out(out, name):
        print name
        for results in out:
            if not results: 
                print 'No results.'
            else:
                for res in results:
                    print res[0]
                    for dup in res[1]:
                        print dup
                    print '\n'
                print '\n'
    #end def
    
#    results = sdb.find_in_db(query, 3)
#    for rid in results:
#        print '\n'
#        print_out(results[rid], rid)
    

    #searcher.find_mp(template[:1000], Primer(SeqRecord(query[:65], id='test'), 0.1e-6), 3)

    #primer = Primer(SeqRecord(query, id='test'), 0.1e-6)
    #out0,out1,out2,out3 = [],[],[],[]
    #cProfile.run("for i in range(10): searcher.find_mp(template, primer, 3);", 'find_mp2.profile')
    #sleep(1)
    #cProfile.run("for i in range(10): searcher.find(template, primer, 3)", 'find.profile')

    
#    ppid = os.getpid()
#    out  = None
#    cProfile.run("for i in range(10):\
#        out = searcher.find_mp(record, Primer(SeqRecord(Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna), id='test'), 0.1e-6), 3)",
#        'find_mp.profile')
#    print_out(out, 'find_mp')
#    cProfile.run("for i in range(10):\
#        out = searcher.find(record, Primer(SeqRecord(Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna), id='test'), 0.1e-6), 3)",
#        'find.profile')
#    print_out(out, 'find')
    
    #print_out(out0, 'find')
    #print_out(out1, 'find_mp')
    print 'Done'
