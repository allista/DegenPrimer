# coding=utf-8
#
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
Created on Oct 29, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import re
import csv
import os
import sys
import argparse
import subprocess
import sqlite3
from tempfile import mkdtemp
from shutil import rmtree
from Bio.SeqUtils import GC


def init_sequence_table(db_cursor, ):
    db_cursor.execute('''CREATE TABLE test_sequence (
                            SID      TEXT,
                            sequence TEXT,
                            Na       FLOAT,
                            Mg       FLOAT,
                            dNTP     FLOAT,
                            DNA      FLOAT,
                            Primer   FLOAT,
                            DMSO     FLOAT,
                            Tm       FLOAT)''')
#end def

def populate_sequence_table(db_cursor, csv_file):
    csv_reader = list(csv.reader(csv_file, delimiter='\t', quotechar='"'))
    for row in csv_reader[1:]:  #skip first row
        db_cursor.execute('''INSERT INTO test_sequence (
                                SID, sequence, 
                                Na, Mg, dNTP, DNA, Primer, DMSO, Tm)
                                VALUES (?,?,?,?,?,?,?,?,?)''', row[:9]) #only the first 9 columns are inserted
#end def


if __name__ == '__main__':
    #parse command line arguments
    parser = argparse.ArgumentParser(description='''Tool to test degen_primer 
                                     facilities using table data in csv format 
                                     as input. First 9 columns of the csv should 
                                     contain following information:
                                     sequence name,
                                     sequence (not degenerate),
                                     Na       (mM),
                                     Mg       (mM),
                                     dNTP     (mM),
                                     DNA      (nM),
                                     Primer   (uM),
                                     DMSO     (%v/v),
                                     Tm       (C).
                                     Output is written into another csv file: 
                                     input_file.csv -> input_file-output.csv''')
    parser.add_argument('data_file', 
                        type=str, nargs='+',
                        help='file(s) with test data in csv format.')
    args = parser.parse_args()

    #initialize the database    
    db = sqlite3.connect(':memory:')
    db_cursor = db.cursor()
    init_sequence_table(db_cursor)
    
    #read in test data from input file(s)
    for data_filename in args.data_file:
        data_file = open(data_filename, 'rb')
        populate_sequence_table(db_cursor, data_file)
        data_file.close()
    
    #prepare tmp directory for degen_primer output files
    wdir    = os.getcwdu()
    tmp_dir = mkdtemp(prefix='degen_primer_test_')
    os.chdir(tmp_dir)
    print 'Working directory is now:', tmp_dir
    
    #output columns and regexps for searching corresponding values in report files 
    output = [('SID','sequence','GC%','length','Tm','Tm-predicted', 'delta-Tm', 'Dimers', 'Dimer dG', 'Hairpins', 'Hairpin dG')]
    values = {'Tm'         : [re.compile(' *(Tm|Tm mean): *(\d{1,3}\.\d{1}) C$'), None],
              'Dimers'     : [re.compile(' *(Dimers count): *(\d{1,3})$'), None],
              'Dimer_dG'   : [re.compile(' *(Most stable dimer): *(-\d{1,3}\.\d{2}) kcal/mol$'), None],
              'Hairpins'   : [re.compile(' *(Hairpins count): *(\d{1,3})$'), None],
              'Hairpin_dG' : [re.compile(' *(Most stable hairpin): *(-\d{1,3}\.\d{2}) kcal/mol$'), None], 
              }
    #test each data row
    db_cursor.execute('SELECT * FROM test_sequence')
    for row in db_cursor:
        #remove possible duplicate report
        report_filename = row[0]+'-short-report.txt'
        if os.path.isfile(report_filename):
            os.remove(report_filename)
        #clear previously found values
        for value in values.values(): value[1] = None
        #construct command line
        degen_primer_cli = 'degen_primer ' +\
                            '--fwd-id "%s" ' % row[0] +\
                            '--sense %s '    % row[1] +\
                            '--Na %f '       % row[2] +\
                            '--Mg %f '       % row[3] +\
                            '--dNTP %f '     % row[4] +\
                            '--DNA %f '      % row[5] +\
                            '--Primer %f '   % row[6] +\
                            '--DMSO %f '     % row[7]
        #run subprocess
        child = subprocess.Popen(degen_primer_cli,
                                 shell=(sys.platform!="win32"))
        child.wait()
        #parse results
        report_filename = row[0]+'-full-report.txt'
        #check if there's a report file
        if not os.path.isfile(report_filename):
            print 'Error: no such file:', report_filename
            output.append((row[0], row[1], len(row[1]), row[8], None, None, None, None, None, None, None))
            continue
        try:
            report_file = open(report_filename, 'r')
        except Exception, e:
            print 'Exception occurred while trying to open file:', e._message
            output.append((row[0], row[1], len(row[1]), row[8], None, None, None, None, None, None, None))
            continue
        #parse report file
        for line in report_file:
            all_values_found = True
            for value in values.values():
                if value[1]: continue #if value is already found, skip regexp
                all_values_found = False
                matches = value[0].match(line)
                if matches:
                    value[1] = matches.group(2)
                    break
            if all_values_found: break
        report_file.close()
        #save results
        output.append((row[0],                                  #SID 
                       row[1],                                  #sequence
                       round(GC(row[1]), 2),                    #GC%
                       len(row[1]),                             #sequence length
                       row[8],                                  #measured Tm
                       values['Tm'][1],                         #predicted Tm
                       round(float(values['Tm'][1])-row[8], 2), #delta Tm
                       values['Dimers'][1],                     #N of dimers
                       values['Dimer_dG'][1],                   #dG of the most stable dimer
                       values['Hairpins'][1],                   #N of hairpins
                       values['Hairpin_dG'][1],                 #dG of the most stable hairpin
                       ))
    #end for
    #write output data
    os.chdir(wdir)
    output_filename = args.data_file[0][:args.data_file[0].rfind('.csv')]+'-output.csv'
    output_file = open(output_filename, 'wb')
    csv_writer = csv.writer(output_file, delimiter='\t', quotechar='"')
    csv_writer.writerows(output)
    output_file.close()
    #cleanup
    rmtree(tmp_dir)
#end main