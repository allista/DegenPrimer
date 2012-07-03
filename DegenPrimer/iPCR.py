'''
Created on Jul 3, 2012

@author: allis
'''

import sys
import subprocess

class iPCR(object):
    '''Wrapper for ipcress process and parser for it's results'''

    def __init__(self, job_id, fwd_primers, rev_primers):
        '''Constructor'''
        self._program_filename = job_id+'.ipcr'
        self._report_filename  = job_id+'-ipcr-report.txt'
        self._fwd_primers = fwd_primers
        self._rev_primers = rev_primers
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
        print '\nAn ipcress program was written to:\n   ', self._program_filename
    #end def
    
    
    def executeProgram(self, fasta_files, max_mismatches):
        ipcr_cli = 'ipcress %s -m %d -s' % (self._program_filename, max_mismatches)
        for fasta_file in fasta_files:
            ipcr_cli += ' "'+fasta_file+'"'
        try:
            child = subprocess.Popen(ipcr_cli,
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=(sys.platform!="win32"))
            ipcr_report = open(self._report_filename, 'w')
            print '\nWriting ipcress report to:\n   ',self._report_filename
            ipcr_report.write(child.stdout.read())
            ipcr_report.close()
        except OSError, e:
            print 'Faild to execute ipcress'
            print_exception(e)
            print 'It seems that "ipcress" executable is not found in the PATH.'
            print 'NOTE: it is provided by the "exonerate" package in debian-like distributions.'
        except Exception, e:
            print 'Faild to execute ipcress'
            print_exception(e)
    #end def
#end class
        