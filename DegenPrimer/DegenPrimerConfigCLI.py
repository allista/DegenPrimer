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
        """Constructor"""
        #setup a list of options and a parser for command line arguments
        self._parser = argparse.ArgumentParser(description=self._description)
        #configuration file
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
            if option[1] not in arg_groups:
                arg_groups[option[1]] = self._parser.add_argument_group(self._groups[option[1]])
            if option[6] == bool:
                arg_groups[option[1]].add_argument(*option[2], dest=option[0], 
                                               help=option[5], default=option[9],
                                               action='store_true')
            else:
                arg_groups[option[1]].add_argument(*option[2], dest=option[0], nargs=option[3], 
                                               metavar=option[4], help=option[5], type=option[6])
    #end def
    
    
    def _override_option(self, option_name):
        value_override = None
        exec_line   = ('if self._args.%(option)s:\n'
                       '    value_override = self._args.%(option)s\n')
        exec (exec_line % {'option':option_name})
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