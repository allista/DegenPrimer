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
Created on Mar 28, 2013

@author: Allis Tauri <allista@gmail.com>
'''


from PipelineTaskBase import PipelineTaskBase


class PrimerTaskBase(PipelineTaskBase):
    
    @staticmethod
    def check_options(args):
        #check if at least one primer is provided
        has_primers = False
        for primer in args.primers:
            if primer is None: continue
            if primer.self_complement:
                print 'Error: %s primer "%s" [%s] is self-complementary.' \
                    % (primer.master_sequence.description, 
                       primer.master_sequence.id, 
                       primer.master_sequence.seq)
                print 'You should not use self-complementary oligonucleotides as primers.\n'
                return False 
            has_primers = True
        if not has_primers:
            print 'At least one primer should be provided.'
            return False
        return True
    #end def
#end class