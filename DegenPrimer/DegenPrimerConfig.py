# coding=utf-8
#
# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
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
from StringTools import random_text
from Primer import Primer, load_sequence
from Option import Option, OptionGroup
from itertools import chain
from copy import deepcopy
                               


class DegenPrimerConfig(object):
    """
    Base class for degen_primer configuration parsers
    """
    
    #program description
    _description = ('This is a tool to compute degenerate primer parameters. '
                    'At least one primer should be provided.')
    
    
    _obsolete_options = \
    [
     OptionGroup('primers',
                 'Primers with IDs and concentrations',
                 options=(Option(name='sense_primer', 
                                 desc='Sense primer',
                                 nargs=1,
                                 py_type=str,
                                 ),
                          Option(name='antisense_primer', 
                                 desc='AntiSense primer',
                                 nargs=1,
                                 py_type=str,
                                 ),
                          Option(name='sense_primer_id', 
                                 desc='Sense primer ID',
                                 nargs=1,
                                 py_type=str,
                                 ),
                          Option(name='antisense_primer_id', 
                                 desc='AntiSense primer ID',
                                 nargs=1,
                                 py_type=str,
                                 ),
                          )),
     OptionGroup('PCR',
                 'PCR conditions',
                 options=(Option(name='S_primer',
                                 desc='Primer concentration',
                                 nargs=1,
                                 py_type=float,
                                 units='uM',
                                 default=0.25,
                                 limits=(0.001, 1e6), #1M max
                                 scale=1e-6,
                                 ),
                          Option(name='A_primer',
                                 desc='Primer concentration',
                                 nargs=1,
                                 py_type=float,
                                 units='uM',
                                 default=0.25,
                                 limits=(0.001, 1e6), #1M max
                                 scale=1e-6,
                                 ),
                          )),
    ]
    
    
    _option_groups = \
    [
        OptionGroup('primers',
                    'Primers with IDs and concentrations',
                    options=(Option(name='primer', 
                                    desc='Primer.',
                                    nargs='+',
                                    py_type=None,
                                    options=(Option(name='sequence',
                                                    desc='Primer sequence (5\'->3\'). '
                                                    'It may be a fasta or genbank file '
                                                    'or simply a raw sequence string composed of '
                                                    'letters of extended IUPAC DNA alphabet.',
                                                    nargs=1,
                                                    py_type=str,
                                                    field_type='file',
                                                    ),
                                             Option(name='id',
                                                    desc='Primer ID string',
                                                    nargs='?',
                                                    py_type=str,
                                                    ),
                                             Option(name='concentration',
                                                    desc='Primer concentration',
                                                    nargs=1,
                                                    py_type=float,
                                                    units='uM',
                                                    default=0.25,
                                                    limits=(0.001, 1e6), #1M max
                                                    scale=1e-6,
                                                    ),
                                             )),
#                             Option(name='only_tm', 
#                                    desc='Only perform Tm calculations',
#                                    nargs='?',
#                                    py_type=bool,
#                                    default=False,
#                                    save=False,
#                                    ),
                             )),
        OptionGroup('BLAST',
                    'BLAST parameters',
                    options=(Option(name='do_blast', 
                                    desc='Do blast search for specificity of primers/primer-pairs. '
                                    'This option must always be set explicitly.',
                                    nargs='?',
                                    py_type=bool,
                                    default=False,
                                    save=False,
                                    ),
                             Option(name='organisms', 
                                    desc='List of organisms or higher taxons to be used in Entrez '
                                    'query in blast searches (e.g. bacteria)',
                                    nargs='+',
                                    py_type=str,
                                    ),
                             )),
        OptionGroup('PCR',
                    'PCR conditions',
                    options=(Option(name='Na', 
                                    desc='Concentration of monovalent ions in PCR buffer',
                                    nargs=1,
                                    py_type=float,
                                    units='mM',
                                    default=50.0,
                                    limits=(1, 5000), #5M max
                                    scale=1e-3,
                                    ),
                             Option(name='Mg', 
                                    desc='Concentration of divalent ions in PCR buffer',
                                    nargs=1,
                                    py_type=float,
                                    units='mM',
                                    default=1.5,
                                    limits=(1, 5000), #5M max
                                    scale=1e-3,
                                    ),
                             Option(name='dNTP', 
                                    desc='Concentration of dNTP in PCR buffer',
                                    nargs=1,
                                    py_type=float,
                                    units='mM',
                                    default=0.1,
                                    limits=(1e-3, 1000), #1M max
                                    scale=1e-3,
                                    ),
                             Option(name='DNA', 
                                    desc='Concentration of target DNA in PCR buffer',
                                    nargs=1,
                                    py_type=float,
                                    units='nM',
                                    default=0.1,
                                    limits=(1e-6, 1e9), #1M max
                                    scale=1e-9,
                                    ),
                             Option(name='DMSO', 
                                    desc='Concentration of DMSO in PCR buffer',
                                    nargs=1,
                                    py_type=float,
                                    units='%% v/v',
                                    default=0.0,
                                    limits=(0, 100), #1M max
                                    scale=1,
                                    ),
                             Option(name='PCR_T', 
                                    desc='Temperature at which primer annealing '
                                    'will be performed during PCR',
                                    nargs=1,
                                    py_type=float,
                                    units='C',
                                    default=60.0,
                                    limits=(20, 80), #1M max
                                    scale=1,
                                    ),
                             )),
        OptionGroup('iPCR',
                    'In silica PCR simulation parameters',
                    options=(Option(name='min_amplicon', 
                                    desc='Minimum amplicon size',
                                    nargs=1,
                                    py_type=int,
                                    units='bp',
                                    default=50,
                                    limits=(1, 1000000000), #5M max
                                    scale=1,
                                    ),
                             Option(name='max_amplicon', 
                                    desc='Maximum amplicon size',
                                    nargs=1,
                                    py_type=int,
                                    units='bp',
                                    default=3000,
                                    limits=(1, 1000000000), #5M max
                                    scale=1,
                                    ),
                             Option(name='max_mismatches', 
                                    desc='Maximum number of mismatches between a primer '
                                    'and a target sequence '
                                    '(default 20 percent of the biggest primer length). '
                                    'It is recommended to gradually increase this parameter '
                                    'up to 6-7, for some annealings with 6 or even 7 mismatches may be '
                                    'more stable than some with 5 mismatches. '
                                    'Note, that PCR simulation with 6 and more mismatches '
                                    'may take a considerable time due to an increasing number of possible products. ',
                                    nargs=1,
                                    py_type=int,
                                    units='bp',
                                    default=None,
                                    limits=(0, 1000000000), #5M max
                                    scale=1,
                                    ),
                             Option(name='polymerase', 
                                    desc='Concentration of polymerase in '
                                    'Units per microliter. A Unit of polimerase is '
                                    'the amount of enzyme which adds 10 nmol of dNTP '
                                    'in 30 minutes to the elongated strand',
                                    nargs=1,
                                    py_type=float,
                                    units='u/ul',
                                    default=0.05,
                                    limits=(1e-3, 1), #1M max
                                    scale=1e6,
                                    ),
                             Option(name='with_exonuclease', 
                                    desc="Set this flag if you are planning to use a "
                                    "DNA-polymerase with 3'-5'-exonuclease activity.",
                                    nargs='?',
                                    py_type=bool,
                                    default=False,
                                    ),
                             Option(name='cycles', 
                                    desc='Number of PCR cycles in simulation',
                                    nargs=1,
                                    py_type=int,
                                    default=30,
                                    limits=(5, 1000), #5M max
                                    scale=1,
                                    ),
                             Option(name='analyse_all_annealings', 
                                    desc='When simulating PCR, include all possible '
                                   'annealing sites of all primers into equilibrium system '
                                   'analysis. This gives more precise estimation of final '
                                   'concentrations of PCR products and sometimes even changes '
                                   'their distribution. But using this option, especially combined '
                                   'with --max-mismatches grater than 5 or 6, usually '
                                   'leads to analysis time of several minutes '
                                   'to tens of minutes.',
                                    nargs='?',
                                    py_type=bool,
                                    default=False,
                                    ),
                             Option(name='fasta_files', 
                                    desc='Path(s) to fasta file(s) containing '
                                    'target sequences for iPCR simulation.',
                                    nargs='+',
                                    py_type=str,
                                    field_type='file',
                                    ),
                             Option(name='sequence_db', 
                                    desc='Path to an sqlite database containing '
                                    'target sequences for iPCR simulation. If provided, '
                                    '--fasta-files option will be ignored.',
                                    nargs=1,
                                    py_type=str,
                                    field_type='file',
                                    ),
                             Option(name='use_sequences', 
                                    desc='ID(s) of the sequence(s) from the database '
                                    'specified by --sequence-db option to be used as target '
                                    'sequences for iPCR cimulation. If empty, all the sequences '
                                    'in the database will be used.',
                                    nargs='+',
                                    py_type=str,
                                    ),
                             )),
        OptionGroup('seq_db',
                    'Management of a sequence database',
                    options=(Option(name='list_db', 
                                    desc='List contents of a sequence database '
                                    'specified.',
                                    nargs=1,
                                    py_type=str,
                                    field_type='file',
                                    save=False,
                                    ),
                             Option(name='add_sequences', 
                                    desc='Path(s) to fasta file(s) containing '
                                    'sequences to be included into the database specified '
                                    'by --sequence-db option. The database itself is '
                                    'created if nesessary.',
                                    nargs='+',
                                    py_type=str,
                                    field_type='file',
                                    save=False,
                                    ),
                             Option(name='del_sequences', 
                                    desc='ID(s) of the sequence(s) to be excluded '
                                    'from the database specified by --sequence-db option.',
                                    nargs='+',
                                    py_type=str,
                                    save=False,
                                    ),
                             )),
        OptionGroup('optimize',
                    'PCR parameters optimization',
                    options=(Option(name='optimize', 
                                    desc='If set, run PCR parameters optimization '
                                    'instead of regular analysis. Warning: this is '
                                    'usually a slow operation.',
                                    nargs='?',
                                    py_type=bool,
                                    default=False,
                                    save=False,
                                    ),
                             Option(name='max_simulations', 
                                    desc='Maximum number of PCR simulations to '
                                    'run during optimization.',
                                    nargs=1,
                                    py_type=int,
                                    default=100,
                                    limits=(10,1000)
                                    ),
                             Option(name='product_purity', 
                                    desc='The grater this value, the grater weight has '
                                    'the purity of the target product(s) during optimization. '
                                    'Zero purity means only product concentration is considered.',
                                    nargs=1,
                                    py_type=int,
                                    default=3,
                                    limits=(0,5)
                                    ),
                             Option(name='target_product', 
                                    desc='PCR product which concentration and purity will be '
                                    'maximized during an optimization. If several products are '
                                    'given, then their sum is maximized.',
                                    nargs='*',
                                    py_type=None,
                                    options=(Option(name='start',
                                                       desc='Product start as defined in PCR simulation report.',
                                                       nargs=1,
                                                       py_type=int,
                                                       limits=(1,None),
                                                       ),
                                                Option(name='end',
                                                       desc='Product end as defined in PCR simulation report.',
                                                       nargs=1,
                                                       py_type=int,
                                                       limits=(1,None),
                                                       ),
                                                )),
                             Option(name='optimization_parameter', 
                                    desc='PCR parameter to optimize. Optimization is performed '
                                    'in terms of maximization of an objective function which is: '
                                    'sum_over_target_products(product_concentration*product_fraction^product_purity).',
                                    nargs='*',
                                    py_type=None,
                                    options=(Option(name='name',
                                                    desc='Name of the parameter. '
                                                    'Any of PCR conditions as well as '
                                                    'primer and polymerase concentrations may be used. '
                                                    'To optimize primer concentration, '
                                                    'set this to the ID of the primer.',
                                                    nargs=1,
                                                    py_type=str,
                                                    ),
                                             Option(name='min',
                                                    desc='Lower boundary of optimization interval.',
                                                    nargs=1,
                                                    py_type=str,
                                                    ),
                                             Option(name='ini',
                                                    desc='Initial value of the parameter. Optimization '
                                                    'starts from this value and tries to find the nearest '
                                                    'optimum. Thus changing it may change the outcome of '
                                                    'the optimization.',
                                                    nargs=1,
                                                    py_type=str,
                                                    ),
                                             Option(name='max',
                                                    desc='Upper boundary of optimization interval.',
                                                    nargs=1,
                                                    py_type=str,
                                                    ),
                                             Option(name='enabled', 
                                                    desc='If unset, current parameter is not used '
                                                    'in the optimization.',
                                                    nargs=1,
                                                    py_type=bool,
                                                    default=True,
                                                    ),
                                                )),
                             )),
    ]


    def __init__(self):        
        self.job_id         = ''
        self.primers        = []
        
        self._config_file   = None
        self._config        = None

        self._reports       = list()
        self._groups_dict   = dict()
        for i, group in enumerate(self._option_groups):
            self._groups_dict[group.name] = i
        self._options       = list(chain(*(grp.options for grp in self._option_groups)))
        self._options_dict  = dict((opt.dest, opt) for opt in self._options)
    #end def

    @property
    def reports(self): return self._reports
    
    @property
    def config_file(self): return self._config_file
    
    @property
    def options(self):
        opt_names = [opt.dest for opt in self._options]
        return dict((key, value) for key, value in vars(self).items() 
                    if key in opt_names)
    #end def
    
    @classmethod
    def from_options(cls, options):
        new_config = cls()
        for option, value in options.items():
            setattr(new_config, option, value)
        new_config._process_options()
        return new_config
    #end def
    
    
    def find_option(self, full_name, options=None):
        if options is None: options = self._option_groups
        names  = full_name.split('/')
        option = None 
        for opt in options:
            if opt.name == names[0]:
                option = opt
                break
        if len(names) == 1: return option
        if option is None or option.options is None: return None
        return self.find_option('/'.join(names[1:]), option.options)
    #end def
            
        
    def __str__(self):
        ret = ''
        for option in self._options:
            if hasattr(self, option.dest):
                ret += '%s: %s\n' % (option.dest, self._get_option(option))
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
            if hasattr(self, option.dest):
                value = getattr(self, option.dest)
                if type(value) == list: value = tuple(value)
                values.append(value)
        return hash(tuple(values)) & 0xFFFFFFF
    #end def
    
    
    def _get_option(self, option):
        if hasattr(self, option.dest):
            return getattr(self, option.dest)
        else: return None
    #end def
    
    def _copy_option(self, option):
        return deepcopy(self._get_option(option))
    
    def _set_option(self, option, value):
        setattr(self, option.dest, value)
    

    #virtual
    def _override_option(self, option): pass
    
    
    def _apply_to_option(self, value, option, func):
        if value is None: return value
        if option.is_compound:
            values = value if option.is_poly else [value]
            for sub in option.options:
                for val in values:
                    if sub.name in val:
                        val[sub.name] = self._apply_to_option(val[sub.name], sub, func)
            return value
        return func(value, option)
    #end def
    
    
    def _check_limits(self, option):
        def check_limits(value, option):
            if option.limits is None: return value
            if option.limits[0] and value < option.limits[0] \
            or option.limits[1] and value > option.limits[1]:
                raise ValueError(('The value of "%s" parameter should be '
                                  'within %s %s') % (option.name, 
                                                     option.limits_str,
                                                     option.units))
            return value
        value = self._get_option(option)
        self._apply_to_option(value, option, check_limits)
    #end def
    
    
    def _scale_value(self, option):
        def scale_value(value, option):
            if option.scale is None: return value
            return value*option.scale
        value = self._get_option(option)
        self._set_option(option, self._apply_to_option(value, option, scale_value))
    #end def
        
    
    def _unscaled_value(self, option):
        def unscale_value(value, option):
            if option.scale is None: return value
            return value/option.scale
        value = self._copy_option(option)
        return self._apply_to_option(value, option, unscale_value)
    #end def
    
    
    def _fill_option(self, option):
        #setup class member for the option
        self._set_option(option, option.default)
        #try to override default value
        value_override = self._override_option(option)
        if value_override is not None:
            self._set_option(option, value_override)
        #if no value given, try to read in from config file
        elif self._config:
            if  self._config.has_option(option.group, option.dest) \
            and self._config.get(option.group, option.dest) != 'None':
                #read value from configuration
                if   option.py_type == str:
                    self._set_option(option, self._config.get(option.group, option.dest).decode('UTF-8'))
                elif option.py_type == float:
                    self._set_option(option, self._config.getfloat(option.group, option.dest))
                elif option.py_type == int:
                    self._set_option(option, self._config.getint(option.group, option.dest))
                elif option.py_type == bool:
                    self._set_option(option, self._config.getboolean(option.group, option.dest))
                else: self._set_option(option, self._config.get(option.group, option.dest).decode('UTF-8'))
                #if it is a list, evaluate it
                if option.is_compound or option.is_multi:
                    self._set_option(option, eval(self._get_option(option)))
        #check if option value is within limits
        self._check_limits(option)
        #scale option
        self._scale_value(option)
    #end def
    
    
    def _convert_obsolete_options(self):
        for option in chain(*(grp.options for grp in self._obsolete_options)):
            self._fill_option(option)
        if not self.primer_list:
            try:
                if self.sense_primer:
                    primer = dict()
                    primer['sequence'] = self.sense_primer
                    primer['id'] = self.sense_primer_id
                    primer['concentration'] = self.S_primer
                    self.primer_list.append(primer)
                if self.antisense_primer:
                    primer = dict()
                    primer['sequence'] = self.antisense_primer
                    primer['id'] = self.antisense_primer_id
                    primer['concentration'] = self.A_primer
                    self.primer_list.append(primer)
            except AttributeError: pass
    #end def    
    
    
    def _process_options(self):
        '''Compile all attributes that use information in options'''
        #load primers
        try: 
            self.primers = []
            if self.primer_list:
                for primer in self.primer_list:
                    seq_record = load_sequence(primer['sequence'], 
                                               primer['id'])
                    if not seq_record.id:
                        seq_record.id = random_text(6)
                    self.primers.append(Primer(seq_record, primer['concentration']))
                    primer['sequence'] = str(seq_record.seq)
                    primer['id']       = seq_record.id
        except AttributeError: pass
        #set max_mismatches to be the 20% of the length of the biggest primer
        try:
            if self.max_mismatches is None:
                self.max_mismatches = 1
                for primer in self.primers:
                    if primer is None: continue
                    self.max_mismatches = max(int(0.2*len(primer)), 
                                              self.max_mismatches)
        except AttributeError: pass
        #set job_id
        self.job_id = ''
        for primer in self.primers:
            if not primer: continue
            if self.job_id: self.job_id += '-'
            if primer.id:
                self.job_id += primer.id
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
        #process filled options
        self._convert_obsolete_options()
        self._process_options()
    #end def
    
    
    def reset_options(self): self.parse_configuration(None)
    
    def reset_temporary_options(self):
        '''Reset only options which are not saved into a config file'''
        for option in self._options:
            if not option.save:
                self._set_option(option, option.default)
        #process filled options
        self._process_options()
    #end def
    
    
    def save_configuration(self, silent=False):
        config = SafeConfigParser()
        config.optionxform = str
        for option in self._options:
            #do not save options not marked to be saved
            if not option.save: continue
            #save all other
            if not config.has_section(option.group):
                config.add_section(option.group)
            value = self._unscaled_value(option)
            config.set(option.group, option.dest, unicode(value).encode('UTF-8'))
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
            print '\nConfiguration was written to:\n      %s' % self._config_file
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
    conf.parse_configuration('../A9F-A915R.cfg')
    opts = conf.options
    print conf.find_option('primers/primer/sequence').full_name
    conf.from_options(opts)
    conf.job_id = 'test'
    conf.save_configuration()
