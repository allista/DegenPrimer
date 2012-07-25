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
Created on Jul 1, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from ConfigParser import SafeConfigParser
from StringTools import random_text, print_exception
from OligoFunctions import load_sequence

class DegenPrimerConfig(object):
    """
    Base class for degen_primer configuration parsers
    """
    
    #program description
    _description = 'This is a tool to compute degenerate ' \
                   'primer parameters. At least one primer (sense or ' \
                   'antisense) should be provided.'
    
    #directory of option goups
    #          #name
    _groups  = {'primers':'Primers with IDs',
               'PCR'    :'PCR conditions for Tm and dG calculation',
               'iPCR'   :'In silica PCR simulation parameters',
               'BLAST'  :'BLAST parameters'}

    #dictionary of options
    #           name            group
    _options = [                #primers with ids
               ('sense_primer','primers',
                               #command-line arguments 
                               ('-s', '--sense', '-f', '--forward'),
                               #number of arguments
                               1,
                               #metavar
                               'SEQUENCE',
                               #help string
                               'A sense primer sequence (5\'->3\'). '
                               'It may be a fasta or genbank file '
                               'or simply a raw sequence string composed of '
                               'letters of extended IUPAC DNA alphabet',
                               #type
                               str, #for argparse
                               's', #for string formatting
                               ('string','file'), #for gui
                               #default value
                               None),
               ('antisense_primer','primers',
                               #command-line arguments 
                               ('-a', '--antisense', '-r', '--reverse'),
                               #number of arguments
                               1,
                               #metavar
                               'SEQUENCE',
                               #help string
                               'An antisense primer sequence (5\'->3\'). '
                               'It may be a fasta or genbank file '
                               'or simply a raw sequence string composed of '
                               'letters of extended IUPAC DNA alphabet',
                               #type
                               str, #for argparse
                               's', #for string formatting
                               ('string','file'), #for gui
                               #default value
                               None),
               ('sense_primer_id','primers',
                               #command-line arguments 
                               ('--sense-id', '--fwd-id'),
                               #number of arguments
                               1,
                               #metavar
                               'ID',
                               #help string
                               'A sense primer identifier',
                               #type
                               str, #for argparse
                               's', #for string formatting
                               ('string',), #for gui
                               #default value
                               None),
               ('antisense_primer_id','primers',
                               #command-line arguments 
                               ('--antisense-id', '--rev-id'),
                               #number of arguments
                               1,
                               #metavar
                               'ID',
                               #help string
                               'An antisense primer identifier',
                               #type
                               str, #for argparse
                               's', #for string formatting
                               ('string',), #for gui
                               #default value
                               None),
                #Tm and dG calculations
               ('Na','PCR',
                               #command-line arguments 
                               ('--Na',),
                               #number of arguments
                               1,
                               #metavar
                               'C(Na) mM',
                               #help string
                               'Concentration of monovalent ions in mM for '
                               'Tm and dG correction (def=50)',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               50.0),
               ('Mg','PCR',
                               #command-line arguments 
                               ('--Mg',),
                               #number of arguments
                               1,
                               #metavar
                               'C(Na) mM',
                               #help string
                               'Concentration of divalent ions in mM for '
                               'Tm and dG correction (def=1.5)',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               1.5),
               ('dNTP','PCR',
                               #command-line arguments 
                               ('--dNTP',),
                               #number of arguments
                               1,
                               #metavar
                               'C(Na) mM',
                               #help string
                               'Concentration of dNTP in mM for '
                               'Tm and dG correction (def=0)',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               0),
               ('DNA','PCR',
                               #command-line arguments 
                               ('--DNA',),
                               #number of arguments
                               1,
                               #metavar
                               'C(Na) nM',
                               #help string
                               'Concentration of target DNA in nM for '
                               'Tm and dG correction (def=50)',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               0.25),
               ('Primer','PCR',
                               #command-line arguments 
                               ('--Primer',),
                               #number of arguments
                               1,
                               #metavar
                               'C(Na) uM',
                               #help string
                               'Concentration of primer (assume C(sense)=C(antisense)) '
                               'in uM for Tm and dG correction (def=0.25)',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               50),
               ('sec_dG','PCR',
                               #command-line arguments 
                               ('--sec-dG',),
                               #number of arguments
                               1,
                               #metavar
                               'kcal/mol',
                               #help string
                               'Dimers with free energy ABOVE this threshold will not '
                               'be reported in SHORT report (default is -5 kcal/mol. '
                               'For hairpins the threshold is grater by 2 kcal/mol. '
                               'For 3\' structures corresponding thresholds are grater by '
                               'another 1 kcal/mol',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               -5.0),
               #in silica PCR
               ('min_amplicon','iPCR',
                               #command-line arguments 
                               ('--min-amplicon',),
                               #number of arguments
                               1,
                               #metavar
                               'bp',
                               #help string
                               'Minimum amplicon size (default 50). Applies to '
                               'both ipcress simulation and blast results parsing.',
                               #type
                               int, #for argparse
                               'd', #for string formatting
                               ('integer',), #for gui
                               #default value
                               50),
               ('max_amplicon','iPCR',
                               #command-line arguments 
                               ('--max-amplicon',),
                               #number of arguments
                               1,
                               #metavar
                               'bp',
                               #help string
                               'Maximum amplicon size (default 3000). Applies to '
                               'both ipcress simulation and blast results parsing.',
                               #type
                               int, #for argparse
                               'd', #for string formatting
                               ('integer',), #for gui
                               #default value
                               3000),
               ('max_mismatches','iPCR',
                               #command-line arguments 
                               ('--max-mismatches',),
                               #number of arguments
                               1,
                               #metavar
                               'bp',
                               #help string
                               'Maximum number of mismatches between a primer '
                               'and a target sequence '
                               '(default 20%% of the biggest primer length). '
                               'Applies only to ipcress simulation.',
                               #type
                               int, #for argparse
                               'd', #for string formatting
                               ('integer',), #for gui
                               #default value
                               None),
               ('hsp_dG','iPCR',
                               #command-line arguments 
                               ('--hsp-dG',),
                               #number of arguments
                               1,
                               #metavar
                               'kcal/mol',
                               #help string
                               'HSPs with free energy ABOVE this threshold '
                               'will be considered unstable, not yielding PCR products '
                               '(default -10 kcal/mol). Applies only to blast results parsing.',
                               #type
                               float, #for argparse
                               'f', #for string formatting
                               ('float',), #for gui
                               #default value
                               -10.0),
               ('no_exonuclease','iPCR',
                               #command-line arguments 
                               ('--no-exonuclease',),
                               #number of arguments
                               1,
                               #metavar
                               None,
                               #help string
                               "Set this flag if you are planning to use a "
                               "DNA-polymerase without 3'-5'-exonuclease activity. "
                               "Applies only to blast results parsing.",
                               #type
                               bool, #for argparse
                               'd', #for string formatting
                               ('boolean',), #for gui
                               #default value
                               -10.0),
               ('fasta_files','iPCR',
                               #command-line arguments 
                               ('--fasta-files',),
                               #number of arguments
                               '+',
                               #metavar
                               'path',
                               #help string
                               'Path(s) to fasta file(s) containing target sequences. '
                               'If fasta files are provided, ipcress simulation will be '
                               'launched automatically. Applies only to ipcress simulation.',
                               #type
                               str, #for argparse
                               's', #for string formatting
                               ('string','file'), #for gui
                               #default value
                               None),
               #BLAST
               ('do_blast','BLAST',
                               #command-line arguments 
                               ('--do-blast',),
                               #number of arguments
                               1,
                               #metavar
                               None,
                               #help string
                               'Do blast search for specificity of primers/primer-pairs. '
                               'This option must be always set explicitly.',
                               #type
                               bool, #for argparse
                               'd', #for string formatting
                               ('boolean',), #for gui
                               #default value
                               False),
               ('organisms','BLAST',
                               #command-line arguments 
                               ('--organisms',),
                               #number of arguments
                               '+',
                               #metavar
                               None,
                               #help string
                               'List of organisms or higher taxons to be used in Entrez '
                               'query in blast searches (e.g. bacteria)',
                               #type
                               str, #for argparse
                               's', #for string formatting
                               ('string',), #for gui
                               #default value
                               None),
              ]


    #base class constructor, there's nothing to initialize yet
    def __init__(self):
        pass
    
    
    #virtual
    def _override_option(self, option_name):
        return None
    
    
    def _fill_option(self, option):
        option_dict = {'section':option[1], 
                       'option' :option[0], 
                       'type'   :option[7], 
                       'value'  :option[9]}
        #setup class member for the option
        if option_dict['value'] == None: option_line = 'self.%(option)s = None'
        else: option_line = 'self.%(option)s = %(value)'+('%(type)s\n' % option_dict) 
        exec (option_line % option_dict)
        #try to override default value
        value_override = self._override_option(option[0])
        if value_override:
            exec_line   = ('self.%(option)s = value_override\n')
            exec (exec_line % option_dict)
            return
        #if failed, try to read in config file
        exec_line   = ('if self._config '
                       'and self._config.has_option("%(section)s","%(option)s") '
                       'and self._config.get("%(section)s","%(option)s") != "None":\n')
        if   option_dict['type'] == 's':
            exec_line += '    self.%(option)s = self._config.get("%(section)s","%(option)s")\n'
        elif option_dict['type'] == 'f':
            exec_line += '    self.%(option)s = self._config.getfloat("%(section)s","%(option)s")\n'
        elif option_dict['type'] == 'd':
            exec_line += '    self.%(option)s = self._config.getint("%(section)s","%(option)s")\n'
        try: 
            exec (exec_line % option_dict)
        except Exception, e:
            print 'DegenPrimerConfig._fill_option:'
            print_exception(e) 
            pass
    #end def
    
    
    def parse_configuration(self, config_file=None):
        #read in configuration file if it is provided
        self._config_file = config_file
        self._config      = None
        if self._config_file:
            self._config = SafeConfigParser()
            if not self._config.read(self._config_file):
                self._config = None
        
        #fill in the configuration
        for option in self._options:
            self._fill_option(option)
        
        #customize some options    
        #set max_mismatches to be the 20% of the length of the smallest primer
        if not self.max_mismatches:
            if self.sense_primer:
                self.max_mismatches = int(0.2*len(self.sense_primer))
            if self.antisense_primer:
                self.max_mismatches = max(int(0.2*len(self.antisense_primer)), 
                                          self.max_mismatches)
            if not self.max_mismatches:
                self.max_mismatches = 1

        #if 'fasta_files' was read in from the config file, convert it to the list
        if type(self.fasta_files) == str:
            self.fasta_files = eval(self.fasta_files)
                
        #load primers
        self.primers = [[],[]]
        if self.sense_primer:
            seq_record = load_sequence(self.sense_primer, self.sense_primer_id, 'sense')
            self.primers[0] = [seq_record, [], dict()]
        if self.antisense_primer:
            seq_record = load_sequence(self.antisense_primer, self.antisense_primer_id, 'antisense')
            self.primers[1] = [seq_record, [], dict()]
            
        #set job_id
        #if there's no ids provided, generate a random one
        self.job_id = ''
        for primer in self.primers:
            if not primer: continue
            if self.job_id != '': self.job_id += '-'
            if primer[0].id:
                self.job_id += primer[0].id
        if not self.job_id: 
            self.job_id = random_text(6)+'-'+random_text(6)
    #end def
    
    
    def save_configuration(self):
        config = SafeConfigParser()
        config.optionxform = str
        for option in self._options:
            #do not save 'do_blast' option
            if option[0] == 'do_blast': continue
            #save all other
            if not config.has_section(option[1]):
                config.add_section(option[1])
            exec 'config.set("%(section)s", "%(option)s", ' \
                 'str(self.%(option)s))' % {'section': option[1],
                                            'option' : option[0]}
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
            exec 'print "%(option)s:", self.%(option)s' % {'option':option[0]}
#end class


#tests
if __name__ == '__main__':
    conf = DegenPrimerConfig()
    conf.parse_configuration()
    conf.print_options()
    print 'job_id:', conf.job_id
    conf.save_configuration()