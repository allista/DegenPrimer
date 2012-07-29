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
               {'option':'sense_primer',
                               'section'   :'primers',
                               #command-line arguments 
                               'args'      :('-s', '--sense', '-f', '--forward'),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'SEQUENCE',
                               #help string
                               'help'      :'A sense primer sequence (5\'->3\'). '
                               'It may be a fasta or genbank file '
                               'or simply a raw sequence string composed of '
                               'letters of extended IUPAC DNA alphabet',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None)},
               {'option':'antisense_primer',
                               'section'   :'primers',
                               #command-line arguments 
                               'args'      :('-a', '--antisense', '-r', '--reverse'),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'SEQUENCE',
                               #help string
                               'help'      :'An antisense primer sequence (5\'->3\'). '
                               'It may be a fasta or genbank file '
                               'or simply a raw sequence string composed of '
                               'letters of extended IUPAC DNA alphabet',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None)},
               {'option':'sense_primer_id',
                               'section'   :'primers',
                               #command-line arguments 
                               'args'      :('--sense-id', '--fwd-id'),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'ID',
                               #help string
                               'help'      :'A sense primer identifier',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'string', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None)},
               {'option':'antisense_primer_id',
                               'section'   :'primers',
                               #command-line arguments 
                               'args'      :('--antisense-id', '--rev-id'),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'ID',
                               #help string
                               'help'      :'An antisense primer identifier',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'string', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None)},
                #Tm and dG calculations
               {'option':'Na',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--Na',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'mM',
                               #help string
                               'help'      :'Concentration of monovalent ions in mM for '
                               'Tm and dG correction (def=50)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :50.0,
                               'limits'    :(50, 1100)},
               {'option':'Mg',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--Mg',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'mM',
                               #help string
                               'help'      :'Concentration of divalent ions in mM for '
                               'Tm and dG correction (def=1.5)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :1.5,
                               'limits'    :(0, 1000)},
               {'option':'dNTP',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--dNTP',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'mM',
                               #help string
                               'help'      :'Concentration of dNTP in mM for '
                               'Tm and dG correction (def=0)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0,
                               'limits'    :(0, 1000)},
               {'option':'DNA',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--DNA',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'nM',
                               #help string
                               'help'      :'Concentration of target DNA in nM for '
                               'Tm and dG correction (def=50)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :50,
                               'limits'    :(0, 1000)},
               {'option':'Primer',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--Primer',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'uM',
                               #help string
                               'help'      :'Concentration of primer (assume C(sense)=C(antisense)) '
                               'in uM for Tm and dG correction (def=0.25)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0.25,
                               'limits'    :(0, 1000)},
               {'option':'sec_dG',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--sec-dG',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'kcal/mol',
                               #help string
                               'help'      :'Dimers with free energy ABOVE this threshold will not '
                               'be reported in SHORT report (default is -5 kcal/mol. '
                               'For hairpins the threshold is grater by 2 kcal/mol. '
                               'For 3\' structures corresponding thresholds are grater by '
                               'another 1 kcal/mol',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :-5.0,
                               'limits'    :(-1000, 0)},
               #in silica PCR
               {'option':'min_amplicon',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--min-amplicon',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'bp',
                               #help string
                               'help'      :'Minimum amplicon size (default 50). Applies to '
                               'both ipcress simulation and blast results parsing.',
                               #type
                               'py_type'   :int, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'integer', #for gui
                               #default value
                               'default'   :50,
                               'limits'    :(1, 1000000000)},
               {'option':'max_amplicon',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--max-amplicon',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'bp',
                               #help string
                               'help'      :'Maximum amplicon size (default 3000). Applies to '
                               'both ipcress simulation and blast results parsing.',
                               #type
                               'py_type'   :int, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'integer', #for gui
                               #default value
                               'default'   :3000,
                               'limits'    :(1, 1000000000)},
               {'option':'max_mismatches',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--max-mismatches',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'bp',
                               #help string
                               'help'      :'Maximum number of mismatches between a primer '
                               'and a target sequence '
                               '(default 20 percent of the biggest primer length). '
                               'Applies only to ipcress simulation.',
                               #type
                               'py_type'   :int, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'integer', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(0, 1000000000)},
               {'option':'hsp_dG',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--hsp-dG',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'kcal/mol',
                               #help string
                               'help'      :'HSPs with free energy ABOVE this threshold '
                               'will be considered unstable, not yielding PCR products '
                               '(default -10 kcal/mol). Applies only to blast results parsing.',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :-10.0,
                               'limits'    :(-1000, 0)},
               {'option':'no_exonuclease',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--no-exonuclease',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :None,
                               #help string
                               'help'      :"Set this flag if you are planning to use a "
                               "DNA-polymerase without 3'-5'-exonuclease activity. "
                               "Applies only to blast results parsing.",
                               #type
                               'py_type'   :bool, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'boolean', #for gui
                               #default value
                               'default'   :False,
                               'limits'    :(None, None)},
               {'option':'fasta_files',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--fasta-files',),
                               #number of arguments
                               'nargs'     :'+',
                               #metavar
                               'metavar'   :'path',
                               #help string
                               'help'      :'Path(s) to fasta file(s) containing target sequences. '
                               'If fasta files are provided, ipcress simulation will be '
                               'launched automatically. Applies only to ipcress simulation.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None)},
               #BLAST
               {'option':'do_blast',
                               'section'   :'BLAST',
                               #command-line arguments 
                               'args'      :('--do-blast',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :None,
                               #help string
                               'help'      :'Do blast search for specificity of primers/primer-pairs. '
                               'This option must be always set explicitly.',
                               #type
                               'py_type'   :bool, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'boolean', #for gui
                               #default value
                               'default'   :False,
                               'limits'    :(None, None)},
               {'option':'organisms',
                               'section'   :'BLAST',
                               #command-line arguments 
                               'args'      :('--organisms',),
                               #number of arguments
                               'nargs'     :'+',
                               #metavar
                               'metavar'   :None,
                               #help string
                               'help'      :'List of organisms or higher taxons to be used in Entrez '
                               'query in blast searches (e.g. bacteria)',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'string', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None)},
              ]


    #base class constructor, there's nothing to initialize yet
    def __init__(self):
        self._reports = list()
        return
    
    
    def _multiple_args(self, option):
        nargs = option['nargs']
        return nargs == '*' or nargs == '+'
    #end def
    
    
    #virtual
    def _override_option(self, option):
        return None
    
    
    def _fill_option(self, option):
        #setup class member for the option
        if option['default'] == None: option_line = 'self.%(option)s = None'
        else: option_line = 'self.%(option)s = %(default)'+('%(str_type)s\n' % option) 
        exec (option_line % option)
        #try to override default value
        value_override = self._override_option(option)
        if value_override:
            exec_line   = ('self.%(option)s = value_override\n')
            exec (exec_line % option)
            return
        #if failed, try to read in config file
        exec_line   = ('if self._config '
                       'and self._config.has_option("%(section)s","%(option)s") '
                       'and self._config.get("%(section)s","%(option)s") != "None":\n')
        if   option['py_type'] == str:
            exec_line += '    self.%(option)s = self._config.get("%(section)s","%(option)s")\n'
        elif option['py_type'] == float:
            exec_line += '    self.%(option)s = self._config.getfloat("%(section)s","%(option)s")\n'
        elif option['py_type'] == int:
            exec_line += '    self.%(option)s = self._config.getint("%(section)s","%(option)s")\n'
        elif option['py_type'] == bool:
            exec_line += '    self.%(option)s = self._config.getboolean("%(section)s","%(option)s")\n'
        try:
            #read value from configuration 
            exec (exec_line % option)
            #if it is a list, evaluate it
            if self._multiple_args(option):
                exec_line   = ('if self.%(option)s:\n'
                               '    self.%(option)s = eval(self.%(option)s)\n')
                exec (exec_line % option)
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
                print 'Unable to load configuration from', self._config_file
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
            if option['option'] == 'do_blast': continue
            #save all other
            if not config.has_section(option['section']):
                config.add_section(option['section'])
            exec 'config.set("%(section)s", "%(option)s", ' \
                 'str(self.%(option)s))' % {'section': option['section'],
                                            'option' : option['option']}
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
            exec 'print "%(option)s:", self.%(option)s' % {'option':option['option']}
            
    
    def register_report(self, report_name, report_file):
        self._reports.append((report_name, report_file))
#end class


#tests
if __name__ == '__main__':
    conf = DegenPrimerConfig()
    conf.parse_configuration()
    conf.print_options()
    print 'job_id:', conf.job_id
    conf.save_configuration()