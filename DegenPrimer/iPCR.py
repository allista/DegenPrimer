#!/usr/bin/python
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
Created on Jul 3, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import os
import subprocess

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from StringTools import wrap_text, time_hr, hr, print_exception
from PCR_Results import PCR_Results
from SecStructures import Duplex

class iPCR(object):
    '''Wrapper for ipcress process and parser for it's results'''

    #garbage string which ipcress inserts in the middle of a target sequence ID
    _ipcress_target_garbage = ':filter(unmasked)'
    

    def __init__(self, 
                 job_id, 
                 fwd_primer, 
                 rev_primer, 
                 min_amplicon, 
                 max_amplicon, 
                 polymerase, 
                 with_exonuclease,
                 num_cycles,
                 side_reactions=None,
                 side_concentrations=None):
        self._ipcress_subprocess  = None
        #files
        self._program_filename    = job_id+'.ipcr'
        self._raw_report_filename = job_id+'-ipcr-raw-report.txt'
        self._PCR_report_filename = job_id+'-ipcr-PCR-report.txt'
        #parameters
        self._fwd_primer         = fwd_primer
        self._rev_primer         = rev_primer
        self._max_mismatches     = None
        self._min_amplicon       = min_amplicon
        self._max_amplicon       = max_amplicon
        self._polymerase         = polymerase
        self._with_exonuclease   = with_exonuclease
        self._num_cycles         = num_cycles
        #results
        self._PCR_products = PCR_Results([fwd_primer, rev_primer],
                                         self._min_amplicon,
                                         self._max_amplicon,
                                         self._polymerase,
                                         self._with_exonuclease,
                                         self._num_cycles)
        if side_reactions:
            self._PCR_products.add_side_reactions(side_reactions)
        if side_concentrations:
            self._PCR_products.add_side_concentrations(side_concentrations)
        self._results      = None
        self._have_results = False
        self._have_report  = False
    #end def
    
    
    def __del__(self):
        if self._ipcress_subprocess != None:
            if self._ipcress_subprocess.poll() == None:
                self._ipcress_subprocess.terminate()
                self._ipcress_subprocess.wait()
    #end def
    
    
    @classmethod
    def _ipcress_executable(cls):
        path = os.environ.get('PATH', None)
        if not path: return None
        for p in os.environ.get('PATH', '').split(os.pathsep):
            p = os.path.join(p, 'ipcress')
            if os.access(p, os.X_OK):
                return p
        return None
    #end def
    
    
    def ipcress_pid(self):
        if self._ipcress_subprocess != None:
            if self._ipcress_subprocess.poll() == None:
                return self._ipcress_subprocess.pid
        return None
    #end def
    
    
    #property functions
    def have_results(self): return self._have_results
    
    
    def write_program(self):
        try:
            ipcr_program = open(self._program_filename, 'w')
        except IOError, e:
            print '\nUnable to open %s' % self._program_filename
            print e.message
            return False
        for fwd_primer in self._fwd_primer.seq_records:
            for rev_primer in self._rev_primer.seq_records:
                pcr_string  = fwd_primer.id + '-' + rev_primer.id
                pcr_string += ' '+str(fwd_primer.seq)
                pcr_string += ' '+str(rev_primer.seq)
                pcr_string += ' '+str(self._min_amplicon)
                pcr_string += ' '+str(self._max_amplicon)
                pcr_string += '\n'
                ipcr_program.write(pcr_string)
        ipcr_program.close()
        print '\nThe ipcress program was written to:\n   ', self._program_filename
        return True
    #end def
    
    
    def execute_program(self, fasta_files, max_mismatches):
        #check if program file exists
        if not os.path.isfile(self._program_filename):
            print '\nFile with ipcress program was not found: %s' % self._program_filename 
            return False
        #check if any of the fasta-files exists
        existing_fasta_files = []
        for ffile in fasta_files:
            if os.path.isfile(ffile):
                existing_fasta_files.append(ffile)
        if not existing_fasta_files: 
            print '\nFailed to execute ipcress: no existing fasta files were provided.'
            return False
        #find ipcress in the path
        self._max_mismatches = max_mismatches
        ipcress_executable = self._ipcress_executable()
        if ipcress_executable is None:
            print '"ipcress" executable is not found in the PATH.'
            print 'NOTE: it is provided by the "exonerate" package in debian-like distributions.'
            return False
        #if found, construct command line and execute
        ipcr_cli = [ipcress_executable]
        ipcr_cli.append('--input')
        ipcr_cli.append(self._program_filename)
        ipcr_cli.append('--mismatch')
        ipcr_cli.append(str(max_mismatches))
        ipcr_cli.append('--sequence')
        for fasta_file in existing_fasta_files:
            ipcr_cli.append(fasta_file)
        try:
            self._ipcress_subprocess = subprocess.Popen(ipcr_cli,
                                                        stdin=subprocess.PIPE,
                                                        stdout=subprocess.PIPE,
                                                        stderr=subprocess.PIPE)
            ipcr_report = open(self._raw_report_filename, 'w')
            output = self._ipcress_subprocess.communicate()[0]
            if self._ipcress_subprocess.returncode != 0:
                print '\nipcress exited with error:\n'
                print output
                return False
            ipcr_report.write(output)
            ipcr_report.close()
            print '\nRaw ipcress report was written to:\n   ',self._raw_report_filename
            self._ipcress_subprocess = None
        except OSError, e:
            print '\nFaild to execute ipcress.'
            print_exception(e)
            return False
        except Exception, e:
            print '\nFaild to execute ipcress.'
            print_exception(e)
            return False
        return True
    #end def
    
    def execute_and_analyze(self, fasta_files, max_mismatches):
        return self.execute_program(fasta_files, max_mismatches) \
        and    self.load_results() \
        and    self.calculate_quantities()
    #end def


    def _clear_target_string(self, target):
        clean_target = ''
        while target:
            pos = target.find(self._ipcress_target_garbage)
            if pos == -1:
                clean_target += target
                target = None
            else:
                clean_target += target[:pos]
                target = target[pos+len(self._ipcress_target_garbage):]
        return clean_target
    #end def
    
    
    def load_results(self):
        #check for file existence
        if not os.path.isfile(self._raw_report_filename): return False
        #open results file
        try:
            ipcr_report = open(self._raw_report_filename, 'r')
        except IOError, e:
            print '\nFailed to open raw iPCR report file:\n   %s' % self._raw_report_filename
            print e.message
            return False
        #parse results file
        self._results = []
        hits = []
        for line in ipcr_report:
            line = line.rstrip('\n')
            if not line: continue
            if line == 'Ipcress result':
                self._results.append(dict())
            else:
                words = line.split()
                if   words[0] == 'Experiment:':
                    self._results[-1]['experiment']   = words[1]
                elif words[0] == 'Target:':
                    target_name = (''.join(w+' ' for w in words[1:])).strip()
                    target_name = self._clear_target_string(target_name)
                    self._results[-1]['hit'] = target_name
                    if not target_name in hits:
                        hits.append(target_name)
                elif words[0] == 'Product:':
                    self._results[-1]['length']  = int(words[1])
                elif words[0] == 'ipcress:':
                    self._results[-1]['start']   = int(words[5])
                    self._results[-1]['end']     = int(words[8])
                elif words[-1] == 'forward':
                    seq = Seq(words[0].strip('.')[::-1], IUPAC.unambiguous_dna).complement()  #forward template 5'->3'
                    self._results[-1]['fwd_seq'] = seq
                elif words[-1] == 'primers':
                    self._results[-1]['fwd_primer'] = Seq(words[0][3:-3], IUPAC.unambiguous_dna) #forward primer 5'->3'
                    self._results[-1]['rev_primer'] = Seq(words[1][3:-3][::-1], IUPAC.unambiguous_dna) #reverse primer 5'->3'
                elif words[-1] == 'revcomp':
                    seq = Seq(words[0].strip('.'), IUPAC.unambiguous_dna).complement() #forward template 5'->3'
                    self._results[-1]['rev_seq'] = seq
        ipcr_report.close()
        if not self._results:
            print '\nNo results found in raw iPCR report:\n   %s' % self._raw_report_filename
            self._results = None
            return False
        return True
    #end def
    
    
    def calculate_quantities(self):
        if not self._results: return False
        #add results to PCR_Results object
        for result in self._results:
            self._PCR_products.add_product(result['hit'], 
                                           result['start']+len(result['fwd_primer']),#ipcress defines start position as 5'-end of the forward primer 
                                           result['end'],
                                           Duplex(result['fwd_primer'], result['fwd_seq']), 
                                           Duplex(result['rev_primer'], result['rev_seq']))
        #compute PCR products quantities
        self._PCR_products.calculate_quantities()
        if not self._PCR_products:
            print '\nAll results found in raw iPCR report were filtered out'
            return False
        self._have_results = True
        return True
    #end def
    
    
    def write_PCR_report(self):
        if not self._have_results: return
        #open report file
        try:
            ipcr_report = open(self._PCR_report_filename, 'w')
        except IOError, e:
            print '\nFailed to open iPCR report file for writing:\n   %s' % self._PCR_report_filename
            print e.message
            return
        #format report
        ipcr_report.write(time_hr())
        #header
        ipcr_report.write(wrap_text('All possible PCR products are ' 
                                     'filtered by amplicon size.\n'
                                     'If --no-exonuclease option was ' 
                                     'provided, products formed by primers with ' 
                                     "mismatches on 3'-end are ignored.\n"
                                     'Relative quantities of the remaining products ' 
                                     'are estimated using equilibrium equations and '
                                     'current PCR parameters.\n'))
        ipcr_report.write('\n')
        #filter parameters
        ipcr_report.write(self._PCR_products.format_report_header())
        if self._max_mismatches != None:
            ipcr_report.write('Number of mismatches allowed: %d\n' % self._max_mismatches)
        #if no PCR products have been found
        if not self._PCR_products:
            ipcr_report.write(hr(' No PCR products have been found ', symbol='!'))
            ipcr_report.close()
            return
        #else...
        ipcr_report.write(self._PCR_products.format_quantity_explanation())
        #PCR products by hit 
        if len(self._PCR_products.hits()) == 1:
            #all products histogram
            hit = self._PCR_products.hits()[0]
            ipcr_report.write(hr(' histogram of all possible PCR products ', symbol='='))
            ipcr_report.write(self._PCR_products.per_hit_header(hit))
            ipcr_report.write(self._PCR_products.per_hit_histogram(hit))
            ipcr_report.write('\n')
            ipcr_report.write(hr(' electrophorogram of all possible PCR products ', symbol='='))
            ipcr_report.write(self._PCR_products.per_hit_electrophoresis(hit))
        else:
            #all products histogram
            ipcr_report.write(hr(' histogram of all possible PCR products ', symbol='='))
            ipcr_report.write(self._PCR_products.all_products_histogram())
            ipcr_report.write('\n\n\n')
            #per hit histogram and phoresis
            ipcr_report.write(hr(' histograms electrophorograms of PCR products of each hit ', symbol='='))
            ipcr_report.write(self._PCR_products.all_graphs_grouped_by_hit())            
        ipcr_report.close()
        print 'iPCR report was written to:\n   ',self._PCR_report_filename
        self._have_report = True
    #end def
    
    
    def reports(self):
        if self._have_report:
            return ({'report_name': 'iPCR report', 'report_file': self._PCR_report_filename},)
        else: return None
    #end def
#end class
        