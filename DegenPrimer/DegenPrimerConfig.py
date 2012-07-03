#!/usr/bin/python
# coding=utf-8
#
# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# indicator_gddccontrol is free software: you can redistribute it and/or modify it
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
Created on Jul 1, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import argparse
from ConfigParser import SafeConfigParser
from StringTools import random_text, print_exception

class DegenPrimerConfig(object):
    """
    Parse command line and/or configuration file and store all job parameters
    """

    def __init__(self):
        """Constructor"""
        #setup a list of options and a parser for command line arguments
        self._options = []
        self._parser = argparse.ArgumentParser(description='This is a tool to compute degenerate \
        primer parameters. At least one primer (sense or antisense) should be provided.')
        #configuration file
        conf_group = self._parser.add_argument_group('Preset configuration')
        conf_group.add_argument('config_file', metavar='file.cfg', 
                            type=str, nargs='?',
                            help='Path to a configuration file containing some '
                                 'or all of the options listed below. If any '
                                 'option is present in the configuration file '
                                 'and on the command line, the latter overrides '
                                 'the former. '
                                 '(NOTE, that the "--do-blast" option '
                                 'always must be set explicitly on the command line.)')
        #primers with ids
        prim_group = self._parser.add_argument_group('Primers with IDs')
        #self._options.append(('SECTION','OPTION', 'VALUE', 'TYPE'))
        self._options.append(('primers','sense_primer', None, 's'))
        prim_group.add_argument('-s', '--sense', '-f', '--forwad', dest='sense_primer', metavar='SEQUENCE', 
                            required=False, type=str,
                            help='A sense primer sequence (5\'->3\') composed of '
                                 'capital letters of extended IUPAC DNA alphabet')
        self._options.append(('primers','sense_primer_id', None, 's'))
        prim_group.add_argument('--sense-id', '--fwd-id', dest='sense_primer_id', metavar='ID', 
                            required=False, type=str,
                            help='A sense primer identification string')
        self._options.append(('primers','antisense_primer', None, 's'))
        prim_group.add_argument('-a', '--antisense', '-r', '--reverse', dest='antisense_primer', metavar='SEQUENCE',
                            required=False, type=str,  
                            help='An antisense primer sequence (5\'->3\') composed of \
                            capital letters of extended IUPAC DNA alphabet')
        self._options.append(('primers','antisense_primer_id', None, 's'))
        prim_group.add_argument('--antisense-id', '--rev-id', dest='antisense_primer_id', metavar='ID', 
                            required=False, type=str,
                            help='An antisense primer identification string')
        #Tm and dG calculations
        TmdG_group = self._parser.add_argument_group('PCR conditions for Tm and dG calculation')
        self._options.append(('PCR','Na', 50, 'f'))
        TmdG_group.add_argument('--Na', metavar='C(Na) mM', 
                            required=False, type=float,
                            help='Concentration of monovalent ions in mM for Tm and dG correction (def=50)')
        self._options.append(('PCR','Mg', 1.5, 'f'))
        TmdG_group.add_argument('--Mg', metavar='C(Mg) mM', 
                            required=False, type=float,
                            help='Concentration of divalent ions in mM for Tm and dG correction (def=1.5)')
        self._options.append(('PCR','dNTP', 0, 'f'))
        TmdG_group.add_argument('--dNTP', metavar='C(Mg) mM', 
                            required=False, type=float,
                            help='Concentration of dNTP in mM for Tm and dG correction (def=0)')
        self._options.append(('PCR','DNA', 50, 'f'))
        TmdG_group.add_argument('--DNA', metavar='C(DNA) nM', 
                            required=False, type=float,
                            help='Concentration of target DNA in nM for Tm and dG correction (def=50)')
        self._options.append(('PCR','Primer', 0.25, 'f'))
        TmdG_group.add_argument('--Primer', metavar='C(Primer) uM', 
                            required=False, type=float,
                            help='Concentration of primer (assume C(sense)=C(antisense)) \
                            in uM for Tm and dG correction (def=0.25)')
        #in silica PCR
        iPCR_group = self._parser.add_argument_group('In silica PCR simulation parameters')
        self._options.append(('iPCR','min_amplicon', 50, 'd'))
        iPCR_group.add_argument('--min-amplicon', metavar='bp', 
                            required=False, type=int,
                            help='Minimum amplicon size for ipcress simulation (default 50)')
        self._options.append(('iPCR','max_amplicon', 3000, 'd'))
        iPCR_group.add_argument('--max-amplicon', metavar='bp', 
                            required=False, type=int,
                            help='Maximum amplicon size for ipcress simulation (default 3000)')
        self._options.append(('iPCR','fasta_files', None, 's'))
        iPCR_group.add_argument('--fasta-files', metavar='path', 
                            required=False, nargs='+', type=str,
                            help='Path(s) to fasta files containing target sequences for ipcress simulation. '
                            'If fasta files are provided, ipcress simulation will be launched automatically.')
        self._options.append(('iPCR','max_mismatches', 1, 'd'))
        iPCR_group.add_argument('--max-mismatches', metavar='b', 
                            required=False, type=int,
                            help='Maximum mismatches for ipcress simulation (default 1)')
        #BLAST
        BLAST_group = self._parser.add_argument_group('BLAST parameters')
        self._options.append(('BLAST','do_blast', False, 'd'))
        BLAST_group.add_argument('--do-blast', default=False, 
                            required=False, action='store_true', 
                            help='Do blast search for specificity of primers/primer-pairs. '
                            'This option must be always set explicitly.')
        self._options.append(('BLAST','organisms', None, 's'))
        BLAST_group.add_argument('--organisms', metavar='name', 
                            required=False, type=str, nargs='+',
                            help='List of organisms or higher taxons to be used in Entrez \
                            query in blast searches (e.g. bacteria)')
        self._options.append(('BLAST','top_hits', 10, 'd'))
        BLAST_group.add_argument('--top-hits', dest='top_hits', metavar='N',
                            required=False, type=int,
                            help='Number of top hits to include in the report')
        self._options.append(('BLAST','top_hsps', 4, 'd'))
        BLAST_group.add_argument('--top-hsps', dest='top_hsps', metavar='N',
                            required=False, type=int,
                            help='Number of top HSPs of each hit to include in the report')
        #output
        output_group= self._parser.add_argument_group('Output parameters')
        self._options.append(('output','dG_threshold', -5.0, 'f'))
        output_group.add_argument('--dG-threshold', metavar='kcal/mol', 
                            required=False, type=float,
                            help='Dimers with free energy ABOVE this threshold will not '
                            'be reported in SHORT report (default is -5 kcal/mol. '
                            'For hairpins the threshold is grater by 2 kcal/mol. '
                            'For 3\' structures corresponding thresholds are grater by '
                            'another 1 kcal/mol' )
    #end def
    
    
    def _fill_option(self, option):
        option_dict = {'section':option[0], 
                       'option' :option[1], 
                       'value'  :option[2], 
                       'type'   :option[3]}
        #setup class member for the option
        option_line = 'self.%(option)s = %(value)'+('%(type)s\n' % option_dict)
        exec (option_line % option_dict)
        #try to read in command line
        exec_line   = ('if self._args.%(option)s:\n'
                       '    self.%(option)s = self._args.%(option)s\n')
        exec (exec_line % option_dict)
        #if not, try to read in config file
        exec_line   = ('if not self._args.%(option)s and self._config '
                       'and self._config.has_option("%(section)s","%(option)s") '
                       'and self._config.get("%(section)s","%(option)s") != "None":\n')
        if   option[3] == 's':
            exec_line += '    self.%(option)s = self._config.get("%(section)s","%(option)s")\n'
        elif option[3] == 'f':
            exec_line += '    self.%(option)s = self._config.getfloat("%(section)s","%(option)s")\n'
        elif option[3] == 'd':
            exec_line += '    self.%(option)s = self._config.getint("%(section)s","%(option)s")\n'
        try: 
            exec (exec_line % option_dict)
        except Exception, e:
            print 'DegenPrimerConfig._fill_option:'
            print_exception(e) 
            pass
    #end def
    
    
    def parse_configuration(self):
        #parse command line arguments
        self._args = self._parser.parse_args()
        
        #read in configuration file if it is provided
        self._config = None
        if self._args.config_file:
            self._config = SafeConfigParser()
            if not self._config.read(self._args.config_file):
                self._config = None
        
        #fill in the configuration
        for option in self._options:
            self._fill_option(option)
            
        #set 'do_blast' explicitly from the command line
        self.do_blast = self._args.do_blast
        
        #if fasta_files was read in from the config file, convert it to the list
        if type(self.fasta_files) == str:
            self.fasta_files = eval(self.fasta_files)
            
        #set job_id
        self.job_id = ''
        if self.sense_primer_id:
            self.job_id  = self.sense_primer_id
        if self.antisense_primer_id:
            if self.job_id != '': self.job_id += '-'
            self.job_id += self.antisense_primer_id
        #if there's no ids provided, generate a random one
        if not self.job_id: 
            self.job_id = random_text(6)+'-'+random_text(6)
    #end def
    
    
    def save_configuration(self):
        config = SafeConfigParser()
        config.optionxform = str
        for option in self._options:
            #do not save 'do_blast' option
            if option[1] == 'do_blast': continue
            #save all other
            if not config.has_section(option[0]):
                config.add_section(option[0])
            exec 'config.set("%(section)s", "%(option)s", ' \
                 'str(self.%(option)s))' % {'section': option[0],
                                            'option' : option[1]}
        #write output
        config_filename = self.job_id + '.cfg'
        config_file = open(config_filename, 'wb')
        config.write(config_file)
        config_file.close()
        print 'Configuration was written to:\n   ', config_filename
        print 'NOTE: you may always re-run current analysis with this file \n' \
              '      or use it as a template to run the analysis with modified \n' \
              '      parameters.'
    #end def
    
    
    def print_options(self):
        for option in self._options:
            exec 'print "%(option)s:", self.%(option)s' % {'option':option[1]}
#end class


#tests
if __name__ == '__main__':
    conf = DegenPrimerConfig()
    conf.parse_configuration()
    conf.print_options()
    print 'job_id:', conf.job_id
    conf.save_configuration()