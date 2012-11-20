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

import sys, os
import subprocess
from StringTools import print_exception, wrap_text, time_hr, hr
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from PCR_Results import PCR_Results
from TD_Functions import format_PCR_conditions

class iPCR(object):
    '''Wrapper for ipcress process and parser for it's results'''

    #garbage string which ipcress inserts in the middle of a target sequence ID
    _ipcress_target_garbage = ':filter(unmasked)'
    

    def __init__(self, 
                 job_id, 
                 fwd_primers, 
                 rev_primers, 
                 min_amplicon, 
                 max_amplicon, 
                 no_exonuclease):
        #files
        self._program_filename    = job_id+'.ipcr'
        self._raw_report_filename = job_id+'-ipcr-raw-report.txt'
        self._PCR_report_filename = job_id+'-ipcr-PCR-report.txt'
        #parameters
        self._fwd_primers        = fwd_primers
        self._rev_primers        = rev_primers
        self._max_mismatches     = None
        self._min_amplicon       = min_amplicon
        self._max_amplicon       = max_amplicon
        self._no_exonuclease     = no_exonuclease
        #results
        self._PCR_products = PCR_Results(self._min_amplicon,
                                        self._max_amplicon,
                                        self._no_exonuclease)
        self._results     = None
        self._hits        = None
        self.have_results = False
    #end def
    
    
    def writeProgram(self):
        ipcr_program = open(self._program_filename, 'w')
        for fwd_primer in self._fwd_primers:
            for rev_primer in self._rev_primers:
                pcr_string  = fwd_primer.id + '-' + rev_primer.id
                pcr_string += ' '+str(fwd_primer.seq)
                pcr_string += ' '+str(rev_primer.seq)
                pcr_string += ' '+str(self._min_amplicon)
                pcr_string += ' '+str(self._max_amplicon)
                pcr_string += '\n'
                ipcr_program.write(pcr_string)
        ipcr_program.close()
        print '\nThe ipcress program was written to:\n   ', self._program_filename
    #end def
    
    
    def executeProgram(self, fasta_files, max_mismatches):
        #check if any of the fasta-files exists
        have_fasta = False
        for ffile in fasta_files:
            if os.path.isfile(ffile):
                have_fasta = True
                break
        if not have_fasta: 
            print '\nFailed to execute ipcress: no existing fasta files were provided.'
            return False
        #if so, execute a program
        self._max_mismatches = max_mismatches
        ipcr_cli = 'ipcress %s -m %d -s' % (self._program_filename, max_mismatches)
        for fasta_file in fasta_files:
            ipcr_cli += ' "'+fasta_file+'"'
        try:
            print '\nExecuting iPCR program. This may take awhile...'
            child = subprocess.Popen(ipcr_cli,
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=(sys.platform!="win32"))
            ipcr_report = open(self._raw_report_filename, 'w')
            ipcr_report.write(child.stdout.read())
            ipcr_report.close()
            print '\nRaw ipcress report was written to:\n   ',self._raw_report_filename
        except OSError, e:
            print '\nFaild to execute ipcress'
            print_exception(e)
            print 'It seems that "ipcress" executable is not found in the PATH.'
            print 'NOTE: it is provided by the "exonerate" package in debian-like distributions.'
            return False
        except Exception, e:
            print '\nFaild to execute ipcress'
            print_exception(e)
            return False
        return self.load_results()
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
        #open results file
        try:
            ipcr_report = open(self._raw_report_filename, 'r')
        except Exception:
            print('\nFailed to open raw iPCR report file:\n   %s' % self._raw_report_filename)
            return False
        #parse results file
        self._results = []
        self._hits = []
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
                    if not target_name in self._hits:
                        self._hits.append(target_name)
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
            print('\nNo results found in raw iPCR report:\n   %s' % self._raw_report_filename)
            self._results = None
            return False
        #add results to PCR_Results object
        for result in self._results:
            self._PCR_products.add_product(result['hit'], 
                                          result['start'], 
                                          result['end']+len(result['rev_primer']), #ipcress defines end position as 3'-start of the reverse primer
                                          result['length'], 
                                          result['fwd_seq'], 
                                          result['rev_seq'], 
                                          result['fwd_primer'], 
                                          result['rev_primer'])
        #compute PCR products quantities
        self._PCR_products.compute_quantities()
        if not self._PCR_products:
            print('\nAll results found in raw iPCR report were filtered out\n')
            return False
        self.have_results = True
        return True
    #end def
    
    
    def write_PCR_report(self):
        #open report file
        try:
            ipcr_report = open(self._PCR_report_filename, 'w')
        except Exception:
            print('\nFailed to open iPCR report file for writing:\n   %s' % self._PCR_report_filename)
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
        ipcr_report.write(hr(' filtration parameters '))
        if self._no_exonuclease:
            ipcr_report.write("DNA polymerase DOES NOT have 3'-5'-exonuclease activity\n")
        else: ipcr_report.write("DNA polymerase has 3'-5'-exonuclease activity\n")
        ipcr_report.write('Minimum amplicon size:      %d\n' % self._min_amplicon)
        ipcr_report.write('Maximum amplicon size:      %d\n' % self._max_amplicon)
        if self._max_mismatches:
            ipcr_report.write('Number of mismatches allowed: %d\n' % self._max_mismatches)
        ipcr_report.write('\n')
        ipcr_report.write(hr(' PCR conditions '))
        ipcr_report.write(format_PCR_conditions()+'\n')
        ipcr_report.write('\n\n\n')
        #if no PCR products have been found
        if not self._PCR_products:
            ipcr_report.write(hr(' No PCR products have been found ', symbol='!'))
            ipcr_report.close()
            return
        #else...
        ipcr_report.write(self._PCR_products.format_quantity_explanation())
        #all products histogram
        ipcr_report.write(hr(' histogram of all possible PCR products ', symbol='='))
        ipcr_report.write(self._PCR_products.all_products_histogram())
        ipcr_report.write('\n\n\n')
        #all products electrophoresis
        ipcr_report.write(hr(' electrophorogram of all possible PCR products ', symbol='='))
        ipcr_report.write(self._PCR_products.all_products_electrophoresis())
        ipcr_report.write('\n\n')
        #PCR products by hit, if there more than one hit 
        if len(self._PCR_products.hits()) > 1:
            ipcr_report.write(hr(' the same grouped by target sequence ', symbol='='))
            ipcr_report.write(self._PCR_products.all_graphs_grouped_by_hit())
        ipcr_report.close()
        print '\niPCRess report was written to:\n   ',self._PCR_report_filename
    #end def
    
    
    def register_reports(self, args):
        args.register_report('iPCR report', self._PCR_report_filename)
    #end def
#end class
        