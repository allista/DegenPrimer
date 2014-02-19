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
Created on Mar 10, 2013

@author: Allis Tauri <allista@gmail.com>
'''

from PipelineTaskBase import PipelineTaskBase
from SeqDB import SeqDB
from StringTools import print_dict


class DBManagmentTask(PipelineTaskBase):
    '''
    Given a DegenPrimerConfig object, perform sequence database management 
    operations defined by options in the configuration.
    '''

    @staticmethod
    def check_options(args):
        if args.list_db: return True
        if (args.sequence_db and 
            (args.add_sequences and (len(args.add_sequences) > 1 or args.add_sequences[0] != '') 
             or
             args.del_sequences and (len(args.del_sequences) > 1 or args.del_sequences[0] != ''))):
            return True
        return False
    #end def
    

    def __init__(self, abort_event):
        PipelineTaskBase.__init__(self, abort_event)
        self._seq_db = SeqDB(self._abort_event)
    #end def
    
    
    def _print_db_content(self):
        if not self._seq_db: return
        seq_names = self._seq_db.get_names()
        print '\nContent of the database %s:' % self._seq_db.name
        print '\n', print_dict(seq_names, header=('ID', 'Sequence name'))
    #end def
    
    
    def run(self, args):
        '''Run database managment tasks as specified in the args.'''
        if args.list_db:
            if not self._seq_db.connect(args.list_db): 
                print 'DBManagmentTask: unable to connect to a database %s' % args.list_db
                return -1
            self._print_db_content()
            self._seq_db.close()
            return 1
        elif args.add_sequences:
            if not self._seq_db.connect(args.sequence_db):
                print '\nCreating new database: %s' % args.sequence_db
                self._seq_db.create_db_from_files(args.sequence_db, args.add_sequences) 
            else: self._seq_db.add_files_to_db(args.add_sequences)
            self._print_db_content()
            self._seq_db.close()
            return 1
        elif args.del_sequences:
            if not self._seq_db.connect(args.sequence_db): 
                print 'DBManagmentTask: unable to connect to a database %s' % args.sequence_db
                return -1
            print '\nDeleting sequences %s from database %s' % (', '.join(args.del_sequences), args.sequence_db)
            self._seq_db.del_from_db(args.del_sequences)
            self._print_db_content()
            self._seq_db.close()
            return 1
        return 0
    #end def
#end class