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

import argparse
from DegenPrimerConfig import DegenPrimerConfig


class BoolAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is None: 
            setattr(namespace, self.dest, True)
            return
        setattr(namespace, self.dest, True if values == 'True' else False)
#end class

                
class DegenPrimerConfigCLI(DegenPrimerConfig):
    '''Command line configuration parser for dege_primer CLI'''

    def __init__(self):
        DegenPrimerConfig.__init__(self)
        #setup a list of options and a parser for command line arguments
        self._parser = argparse.ArgumentParser(description=self._description)
        self._args   = None
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
        for group in self._option_groups:
            if group.is_mutex:
                arg_group = self._parser.add_mutually_exclusive_group()
            else: arg_group = self._parser.add_argument_group(group.desc)
            for option in group.options:
                if option.is_compound:
                    arg_group.add_argument(*option.args, dest=option.dest, nargs=option.nargs, action=option.action, 
                                           metavar=option.metavar, help=option.desc, required=option.required)
                elif option.py_type == bool:
                    arg_group.add_argument(*option.args, dest=option.dest, nargs=option.nargs, action=BoolAction,
                                           metavar=option.metavar, choices=('True', 'False'), help=option.desc, required=option.required)
                else:
                    arg_group.add_argument(*option.args, dest=option.dest, nargs=option.nargs,
                                           type=option.py_type, metavar=option.metavar, help=option.desc, required=option.required)
    #end def
    
    
    def _override_option(self, option):
        value_override = None
        if hasattr(self._args, option.dest):
            value_override = getattr(self._args, option.dest)
            #single values should not be lists
            if not option.is_multi and not option.is_poly \
            and type(value_override) == list:
                value_override = value_override[0]
        return value_override
    #end def
    
    
    def parse_configuration(self, config_file = None):
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