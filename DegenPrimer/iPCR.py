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

import sys
import subprocess
import StringTools
from StringTools import print_exception, format_histogram, wrap_text, time_hr, hr
from Electrophoresis import print_electrophoresis

class iPCR(object):
    '''Wrapper for ipcress process and parser for it's results'''

    #garbage string which ipcress inserts in the middle of a target sequence ID
    _ipcress_target_garbage = ':filter(unmasked)'
    
    #titles of a histogram columns
    _col_titles = ('PCR-product', 'relative concentration')


    def __init__(self, job_id, fwd_primers, rev_primers):
        '''Constructor'''
        self._program_filename  = job_id+'.ipcr'
        self._report_filename   = job_id+'-ipcr-full-report.txt'
        self._summary_filename  = job_id+'-ipcr-short-report.txt'
        self._fwd_primers = fwd_primers
        self._rev_primers = rev_primers
        self._results = None
    #end def
    
    
    def writeProgram(self, min_amplicon, max_amplicon):
        ipcr_program = open(self._program_filename, 'w')
        for fwd_primer in self._fwd_primers:
            for rev_primer in self._rev_primers:
                pcr_string  = fwd_primer.id + '-' + rev_primer.id
                pcr_string += ' '+str(fwd_primer.seq)
                pcr_string += ' '+str(rev_primer.seq)
                pcr_string += ' '+str(min_amplicon)
                pcr_string += ' '+str(max_amplicon)
                pcr_string += '\n'
                ipcr_program.write(pcr_string)
        ipcr_program.close()
        print '\nThe ipcress program was written to:\n   ', self._program_filename
    #end def
    
    
    def executeProgram(self, fasta_files, max_mismatches):
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
            ipcr_report = open(self._report_filename, 'w')
            ipcr_report.write(child.stdout.read())
            ipcr_report.close()
            print '\nFull ipcress report was written to:\n   ',self._report_filename
            #parse the results and write a summary
        except OSError, e:
            print 'Faild to execute ipcress'
            print_exception(e)
            print 'It seems that "ipcress" executable is not found in the PATH.'
            print 'NOTE: it is provided by the "exonerate" package in debian-like distributions.'
            return
        except Exception, e:
            print 'Faild to execute ipcress'
            print_exception(e)
            return
        #parse results
        self._parse_results()
        self._write_summary()
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
    
    
    def _parse_results(self):
        self._results = []
        ipcr_report = open(self._report_filename, 'r')
        for line in ipcr_report:
            line = line.rstrip('\n')
            if not line: continue
            if line == 'Ipcress result':
                self._results.append(dict())
                self._results[-1]['quantity'] = 1
            else:
                words = line.split()
                if   words[0] == 'Experiment:':
                    self._results[-1]['name']   = words[1]
                elif words[0] == 'Target:':
                    target_name = (''.join(w+' ' for w in words[1:])).strip()
                    self._results[-1]['target'] = self._clear_target_string(target_name)
                elif words[0] == 'Product:':
                    self._results[-1]['len']    = int(words[1])
                elif words[0] == 'ipcress:':
                    self._results[-1]['left']   = int(words[5])
                    self._results[-1]['right']  = int(words[8])
        ipcr_report.close()
    #end def
    
    
    def _results_histogram(self):
        text_width = StringTools.text_width
        hist_width = len(self._col_titles[1])
        #construct histogram
        histogram  = dict()
        max_target = max(len(r['target']) for r in self._results)
        max_len    = max(len(str(r['len'])) for r in self._results)
        for result in self._results:
            target_spec   = ' %d%s b [%d-%d]' % (result['len'],
                                                 ' '*(max_len-len(str(result['len']))),
                                                 result['left'],
                                                 result['right'])
            target_limit  = text_width-hist_width-2 - len(target_spec)
            target_spacer = max_target - len(result['target'])
            if target_limit < 0: target_limit = 0
            colname = (result['target']+' '*target_spacer)[:target_limit] + target_spec
            if colname in histogram:
                histogram[colname][0] += 1
            else: histogram[colname] = [1,result['len']]
        #normalize histogram
        norm = float(max(histogram.values(), key=lambda x: x[0])[0])
        for colname in histogram:
            histogram[colname][0] /= norm
        #sort histogram by product len
        histogram = sorted(zip(histogram, histogram.values()), key=lambda x: x[1][1])
        histogram = tuple((line[0], line[1][0]) for line in histogram)
        #format histogram
        return format_histogram(histogram, self._col_titles, hist_width)
    #end def
    
    
    def _write_summary(self):
        if not self._results: return
        ipcr_summary = open(self._summary_filename, 'w')
        #format report
        summary_text  = ''
        summary_text += time_hr()
        summary_text += 'A summary report of the in silica PCR simulation.\n'
        summary_text += 'Number of mismatches allowed: %d\n' % self._max_mismatches
        summary_text += '\n'
        summary_text += hr(' iPCR products histogram ')
        summary_text += self._results_histogram()
        summary_text += wrap_text('"relative concentration" of a product is a normalized '
                                  'number of primer pairs that yield this particular product.\n')
        summary_text += '\n'
        summary_text += hr(' iPCR products electrophoresis ')
        summary_text += print_electrophoresis(self._results)
        #write report
        ipcr_summary.write(summary_text)
        ipcr_summary.close()
        print '\nShort ipcress report was written to:\n   ',self._summary_filename
    #end def
#end class
        