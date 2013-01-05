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
        self._db     = None
        self._cursor = None
    #end def
    
    def __nonzero__(self):
        return not self._cursor is None 
        
    
    def _init_db(self, filename):
        #if db exists, back it up
        if os.path.isfile(filename):
            shutil.move(filename, filename+'.back')
        #create new db
        self._db = sqlite3.connect(filename, isolation_level='DEFERRED')
        self._cursor = self._db.cursor()
        #init tables
        self._cursor.execute(self._sequences_schema)
        self._db.commit()
        self.close()
    #end def
    
    
    def create_db(self, filename, sequences):
        '''Create and SQLite database of sequences in a given file.'''
        #create database tables
        self._init_db(filename)
        #connect to the new database
        self.connect(filename)
        self._cursor.execute('PRAGMA cache_size=500000')
        #populate database with data
        for sequence in sequences:
            self._cursor.execute('''
            INSERT INTO sequences (name, sequence, length)
            VALUES (?, ?, ?)''', 
            (sequence.id, sequence.format('fasta'), len(sequence)))
        self._db.commit()
        self.close()
    #end def
    
    
    def connect(self, filename):
        '''Connect to an existent file database of sequences.'''
        if not os.path.isfile(filename): return False
        self._db = sqlite3.connect(filename, isolation_level='DEFERRED')
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
    #end def
    
    
    def get_names(self):
        if self._db is None: return []
        self._cursor.execute('''SELECT name, id from sequences''')
        return list(self._cursor)
    #end def
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
    try:
        record_file = open('Ch5_gnm.fa', 'r')
    except IOError, e:
        print 'Unable to open Ch5_gnm.fa'
        print_exception(e)
        sys.exit(1)
    record = SeqIO.read(record_file, 'fasta', IUPAC.unambiguous_dna)
    record_file.close()
    
    sdb = SeqDB()
    sdb.create_db('test.db', (record,))
    if sdb.connect('test.db'):
        print sdb.get_names()
    
    query = Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna)
    
    searcher = SearchEngine()
    
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
