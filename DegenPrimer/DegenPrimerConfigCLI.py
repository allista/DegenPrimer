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

import argparse
from DegenPrimerConfig import DegenPrimerConfig

class DegenPrimerConfigCLI(DegenPrimerConfig):
    '''Command line configuration parser for dege_primer CLI'''

    def __init__(self):
        DegenPrimerConfig.__init__(self)
        #setup a list of options and a parser for command line arguments
        self._parser = argparse.ArgumentParser(description=self._description)
        #configuration file
        self.config_files = None
        conf_group = self._parser.add_argument_group('Preset configuration')
        conf_group.add_argument('config_files', metavar='file.cfg', 
                                type=str, nargs='*',
                                help='Path to the configuration file(s) containing some '
                                'or all of the options listed below. If more than '
                                'one file is provided, all other options are discarded '
                                'and the program enters the batch mode. Otherwise, if any '
                                'option is present in the configuration file '
                                'and on the command line, the latter overrides '
                                'the former. '
                                '(NOTE, that the "--do-blast" option '
                                'always must be set explicitly on the command line.)')
        #all other options
        arg_groups = dict()
        for option in self._options:
            if option['section'] not in arg_groups:
                arg_groups[option['section']] = self._parser.add_argument_group(self._groups[option['section']])
            if option['py_type'] == bool:
                arg_groups[option['section']].add_argument(*option['args'], dest=option['option'], 
                                                           metavar='True/False', choices=('True', 'False'), help=option['help'])
            else:
                arg_groups[option['section']].add_argument(*option['args'], dest=option['option'], nargs=option['nargs'], 
                                                           type=option['py_type'], metavar=option['metavar'], help=option['help'])
    #end def
    
    
    def _override_option(self, option):
        value_override = None
        exec_line = 'if self._args.%(option)s:\n'
        if self._multiple_args(option):
            exec_line += \
                     '    value_override = self._args.%(option)s\n'
        else: #if option represents a single argument, pass this argument, not a list with it
            exec_line += \
                     '    if type(self._args.%(option)s) == list:\n'
            exec_line += \
                     '        value_override = self._args.%(option)s[0]\n'
            exec_line += \
                     '    else: value_override = self._args.%(option)s\n'
        exec (exec_line % option)
        if value_override != None \
        and option['py_type'] == bool:
            value_override = True if value_override == 'True' else False
        return value_override
    #end def
    
    
    def parse_configuration(self):
        #parse command line arguments
        self._args = self._parser.parse_args()
        self.config_files = self._args.config_files
        #call parent function
        #if a single configuration file was provided, use it
        if self.config_files and len(self.config_files) == 1:
            DegenPrimerConfig.parse_configuration(self, config_file=self.config_files[0])
        #otherwise, only defaults and command-line arguments are used
        else: DegenPrimerConfig.parse_configuration(self)
    #end def
    
    
    def batch(self):
        return (self.config_files and len(self.config_files) > 1)
    #end def
#end class


#tests
if __name__ == '__main__':
    conf = DegenPrimerConfigCLI()
    conf.parse_configuration()
    conf.print_options()
    print 'job_id:', conf.job_id
    conf.save_configuration()