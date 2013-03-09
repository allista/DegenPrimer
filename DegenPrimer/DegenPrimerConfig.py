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
from Primer import Primer, load_sequence

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
                'seq_db' :'Management of a sequence database',
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
                               'limits'    :(None, None),
                               #scale of a numerical argument
                               'scale'     :None,
                               'save'      :True},
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
                               'limits'    :(None, None),
                               #scale of a numerical argument
                               'scale'     :None,
                               'save'      :True},
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
                               'limits'    :(None, None),
                               #scale of a numerical argument
                               'scale'     :None,
                               'save'      :True},
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
                               'limits'    :(None, None),
                               #scale of a numerical argument
                               'scale'     :None,
                               'save'      :True},
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
                               'limits'    :(1, 5000), #5M max
                               #scale of a numerical argument
                               'scale'     :1e-3,
                               'save'      :True},
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
                               'limits'    :(0, 5000), #5M max
                               #scale of a numerical argument
                               'scale'     :1e-3,
                               'save'      :True},
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
                               'Tm and dG correction (def=0.1)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0.1,
                               'limits'    :(1e-3, 1000), #1M max
                               'scale'     :1e-3,
                               'save'      :True},
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
                               'Tm and dG correction (def=1)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :1,
                               'limits'    :(1e-6, 1e9),  #1M max
                               'scale'     :1e-9,
                               'save'      :True},
               {'option':'S_primer',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--sense-primer-concentration --sense-C',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'uM',
                               #help string
                               'help'      :'Concentration of the sense primer '
                               'in uM for Tm and dG correction (def=0.25)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0.25,
                               'limits'    :(0.001, 1e6), #1M max
                               'scale'     :1e-6,
                               'save'      :True},
                {'option':'A_primer',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--antisense-primer-concentration --anti-C',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'uM',
                               #help string
                               'help'      :'Concentration of the antisense primer '
                               'in uM for Tm and dG correction (def=0.25)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0.25,
                               'limits'    :(0.001, 1e6), #1M max
                               'scale'     :1e-6,
                               'save'      :True},
               {'option':'DMSO',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--DMSO',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'%',
                               #help string
                               'help'      :'Concentration of DMSO (%% v/v) '
                               'for Tm and dG correction (def=0)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0,
                               'limits'    :(0, 100),
                               'scale'     :1,
                               'save'      :True},
               {'option':'PCR_T',
                               'section'   :'PCR',
                               #command-line arguments 
                               'args'      :('--PCR-T',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'C',
                               #help string
                               'help'      :'Temperature at which primer annealing '
                               'will be performed during PCR (def. 60 C).',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :60.0,
                               'limits'    :(20, 80),
                               'scale'     :1,
                               'save'      :True},
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
                               'limits'    :(1, 1000000000),
                               'scale'     :1,
                               'save'      :True},
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
                               'limits'    :(1, 1000000000),
                               'scale'     :1,
                               'save'      :True},
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
                               'It is recommended to gradually increase this parameter '
                               'up to 6-7, for some annealings with 6 or even 7 mismatches may be '
                               'more stable than some with 5 mismatches. '
                               'Note, that PCR simulation with 6 and more mismatches '
                               'may take a considerable time due to an increasing number of possible products. '
                               'Applies only to iPCR simulation.',
                               #type
                               'py_type'   :int, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'integer', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(0, 1000000000),
                               'scale'     :1,
                               'save'      :True},
                {'option':'polymerase',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--polymerase_activity',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'u/ul',
                               #help string
                               'help'      :'Concentration of polymerase in '
                               'Units per microliter. A Unit of polimerase is '
                               'the amount of enzyme which adds 10 nmol of dNTP '
                               'in 30 minutes to the elongated strand ' 
                               '(def. 0.05, i.e. 1u in 20ul)',
                               #type
                               'py_type'   :float, #for argparse
                               'str_type'  :'f', #for string formatting
                               'field_type':'float', #for gui
                               #default value
                               'default'   :0.05,
                               'limits'    :(1e-3, 1),
                               'scale'     :1e6,
                               'save'      :True},
               {'option':'with_exonuclease',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--with-exonuclease',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :None,
                               #help string
                               'help'      :"Set this flag to True if you are planning to use a "
                               "DNA-polymerase with 3'-5'-exonuclease activity.",
                               #type
                               'py_type'   :bool, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'boolean', #for gui
                               #default value
                               'default'   :False,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :True},
               {'option':'cycles',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--cycles',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :None,
                               #help string
                               'help'      :'Number of PCR cycles in simulation (def. 30).',
                               #type
                               'py_type'   :int, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'integer', #for gui
                               #default value
                               'default'   :30,
                               'limits'    :(5, 1000),
                               'scale'     :1,
                               'save'      :True},
                {'option':'analyze_all_annealings',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--analyze-all-annealings',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :None,
                               #help string
                               'help'      :'When simulating PCR, include all possible '
                               'annealing sites of all primers into equilibrium system '
                               'analysis. This gives more precise estimation of final '
                               'concentrations of PCR products and sometimes even changes '
                               'their distribution. But using this option, especially combined '
                               'with --max-mismatches grater than 5 or 6, usually '
                               'leads to analysis time of several minutes '
                               'to tens of minutes.',
                               #type
                               'py_type'   :bool, #for argparse
                               'str_type'  :'d', #for string formatting
                               'field_type':'boolean', #for gui
                               #default value
                               'default'   :False,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :True},
               {'option':'fasta_files',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--fasta-files',),
                               #number of arguments
                               'nargs'     :'+',
                               #metavar
                               'metavar'   :'path',
                               #help string
                               'help'      :'Path(s) to fasta file(s) containing '
                               'target sequences for iPCR simulation.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :True},
                {'option':'sequence_db',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--sequence-db',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'path',
                               #help string
                               'help'      :'Path to an sqlite database containing '
                               'target sequences for iPCR simulation. If provided, '
                               '--fasta-files option will be ignored.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :True},
                {'option':'use_sequences',
                               'section'   :'iPCR',
                               #command-line arguments 
                               'args'      :('--use-sequences',),
                               #number of arguments
                               'nargs'     :'+',
                               #metavar
                               'metavar'   :'ids',
                               #help string
                               'help'      :'ID(s) of the sequence(s) from the database '
                               'specified by --sequence-db option to be used as target '
                               'sequences for iPCR cimulation. If empty, all the sequences '
                               'in the database will be used.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'string', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :True},
               #seq_db
               {'option':'list_db',
                               'section'   :'seq_db',
                               #command-line arguments 
                               'args'      :('--list-db',),
                               #number of arguments
                               'nargs'     :1,
                               #metavar
                               'metavar'   :'path',
                               #help string
                               'help'      :'List contents of a sequence database '
                               'specified.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :False},
               {'option':'add_sequences',
                               'section'   :'seq_db',
                               #command-line arguments 
                               'args'      :('--add-sequences',),
                               #number of arguments
                               'nargs'     :'+',
                               #metavar
                               'metavar'   :'path',
                               #help string
                               'help'      :'Path(s) to fasta file(s) containing '
                               'sequences to be included into the database specified '
                               'by --sequence-db option. The database itself is '
                               'created if nesessary.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'file', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :False},
                {'option':'del_sequences',
                               'section'   :'seq_db',
                               #command-line arguments 
                               'args'      :('--del-sequences',),
                               #number of arguments
                               'nargs'     :'+',
                               #metavar
                               'metavar'   :'ids',
                               #help string
                               'help'      :'ID(s) of the sequence(s) to be excluded '
                               'from the database specified by --sequence-db option.',
                               #type
                               'py_type'   :str, #for argparse
                               'str_type'  :'s', #for string formatting
                               'field_type':'string', #for gui
                               #default value
                               'default'   :None,
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :False},
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
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :False},
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
                               'limits'    :(None, None),
                               'scale'     :None,
                               'save'      :True},
              ]


    #base class constructor, there's nothing to initialize yet
    def __init__(self):
        self._config_file   = None
        self._config        = None
        self.primers        = [None, None]
        self.job_id         = ''
        self._reports       = list()
        self.max_mismatches = None
    #end def

    @property
    def reports(self): return self._reports
    
    @property
    def config_file(self): return self._config_file
    

    def __str__(self):
        ret = ''
        for option in self._options:
            if hasattr(self, option['option']):
                ret += '%s: %s\n' % (option['option'], self._get_option(option))
        return ret
    #end def
    
    def __repr__(self): return str(self)    
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __ne__(self, other):
        return self.__hash__() != other.__hash__()
    
    def __hash__(self):
        values = []
        for option in self._options:
            if hasattr(self, option['option']):
                value = getattr(self, option['option'])
                if type(value) == list: value = tuple(value)
                values.append(value)
        return hash(tuple(values)) & 0xFFFFFFF
    #end def
    
    
    def _get_option(self, option):
        if hasattr(self, option['option']):
            value = getattr(self, option['option'])
        else: value = None
        return value
    #end def
    
    def _set_option(self, option, value):
        setattr(self, option['option'], value)
    
    
    @classmethod
    def _multiple_args(cls, option):
        nargs = option['nargs']
        return nargs == '*' or nargs == '+'
    #end def
    
    
    #virtual
    def _override_option(self, option): pass
    
    
    def _check_limits(self, option):
        value = self._get_option(option)
        if value is None: return
        if option['limits'][0] and value < option['limits'][0] \
        or option['limits'][1] and value > option['limits'][1]:
            raise ValueError(('The value of "%(option)s" parameter should be within' % option) + \
                             (' [' + str(option['limits'][0]) + \
                              ', ' + str(option['limits'][1]) + \
                              '] ' + option['metavar']))
    #end def
    
    
    def _scale_value(self, option):
        value = self._get_option(option)
        if option['scale'] is None or value is None: return
        value *= option['scale']
        self._set_option(option, value)
    #end def
        
    
    def _unscale_value(self, option):
        value = self._get_option(option)
        if option['scale'] is None or value is None: return value
        return value/option['scale']
    #end def
    
    
    def _fill_option(self, option):
        #setup class member for the option
        self._set_option(option, option['default'])
        #try to override default value
        value_override = self._override_option(option)
        if value_override != None:
            self._set_option(option, value_override)
        #if no value given, try to read in from config file
        elif self._config:
            if  self._config.has_option(option['section'], option['option']) \
            and self._config.get(option['section'], option['option']) != None:
                #read value from configuration
                if   option['py_type'] == str:
                    self._set_option(option, self._config.get(option['section'], option['option']))
                elif option['py_type'] == float:
                    self._set_option(option, self._config.getfloat(option['section'], option['option']))
                elif option['py_type'] == int:
                    self._set_option(option, self._config.getint(option['section'], option['option']))
                elif option['py_type'] == bool:
                    self._set_option(option, self._config.getboolean(option['section'], option['option']))
                #if it is a list, evaluate it
                if self._multiple_args(option):
                    self._set_option(option, eval(self._get_option(option)))
        #check if option value is within limits
        self._check_limits(option)
        #scale option
        self._scale_value(option)
    #end def
    
    
    def check_configuration(self):
        #check if at least one primer is provided
        has_primers = False
        for primer in self.primers: has_primers |= primer is not None
        if not has_primers:
            print 'At least one primer (sense or antisense) should be provided.'
            return False
        #test for self-complementarity
        for primer in self.primers:
            if primer.self_complement:
                print 'Error: %s primer "%s" [%s] is self-complementary.' \
                    % (primer.master_sequence.description, 
                       primer.master_sequence.id, 
                       primer.master_sequence.seq)
                print 'You should not use self-complementary oligonucleotides as primers.\n'
                return False
        return True
    #end def
    
    
    def parse_configuration(self, config_file=None):
        #read in configuration file if it is provided
        self._config_file = config_file
        self._config      = None
        if self._config_file:
            self._config = SafeConfigParser()
            if not self._config.read(self._config_file):
                print '\nUnable to load configuration from', self._config_file
                self._config = None
        
        #fill in the configuration
        for option in self._options:
            self._fill_option(option)
        
        #load primers
        #if no primer ids are provided, generate a random ones
        if self.sense_primer:
            seq_record = load_sequence(self.sense_primer, 
                                       self.sense_primer_id, 'sense')
            if not seq_record.id:
                seq_record.id = random_text(6)
            self.primers[0] = Primer(seq_record, self.S_primer)
        if self.antisense_primer:
            seq_record = load_sequence(self.antisense_primer, 
                                       self.antisense_primer_id, 'antisense')
            if not seq_record.id:
                seq_record.id = random_text(6)
            self.primers[1] = Primer(seq_record, self.A_primer)
            
        #set max_mismatches to be the 20% of the length of the biggest primer
        if self.max_mismatches == None:
            self.max_mismatches = 1
            for primer in self.primers:
                if primer is None: continue
                self.max_mismatches = max(int(0.2*len(self.primer)), 
                                          self._max_mismatches)

        #set job_id
        self.job_id = ''
        for primer in self.primers:
            if not primer: continue
            if self.job_id: self.job_id += '-'
            if primer.id:
                self.job_id += primer.id
    #end def
    
    
    def save_configuration(self, silent=False):
        config = SafeConfigParser()
        config.optionxform = str
        for option in self._options:
            #do not save options not marked to be saved
            if not option['save']: continue
            #save all other
            if not config.has_section(option['section']):
                config.add_section(option['section'])
            value = self._unscale_value(option)
            config.set(option['section'], option['option'], str(value))
        #write output
        self._config_file = self.job_id + '.cfg'
        try:
            config_file = open(self._config_file, 'wb')
            config.write(config_file)
            config_file.close()
        except IOError, e:
            print '\nFailed to write configuration file:\n   %s' % self._config_file
            print e.message
            return
        if not silent:
            print '\nConfiguration was written to:\n      ', self._config_file
            print ('NOTE: you may always re-run current analysis with this file\n'
                   '      or use it as a template to run the analysis with modified\n'
                   '      parameters.')
    #end def
    
    
    def register_report(self, report_name, report_file):
        self._reports.append((report_name, report_file))
#end class


#tests
if __name__ == '__main__':
    conf = DegenPrimerConfig()
    conf.parse_configuration()
    print conf
    print 'job_id:', conf.job_id
    conf.save_configuration()