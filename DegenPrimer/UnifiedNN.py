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
Created on Jun 24, 2012

@author: Allis Tauri <allista@gmail.com>

All calculations and data are based on:
[1] ﻿SantaLucia, J., & Hicks, D. (2004). 
The thermodynamics of DNA structural motifs. Annual review of biophysics and 
biomolecular structure, 33, 415-40. doi:10.1146/annurev.biophys.32.110601.141800
'''


import os
import csv
from math import log

from SeqUtils import unambiguous_sequences
###############################################################################


def _load_csv(filenames):
    '''Load a csv (tab delimited, quoting character ") file into a dictionary of 
    dicts, assuming that the first row and first column define names, and each cell 
    contains float value. The file is loaded from the first existing file in the 
    given list of paths'''
    #check filepaths provided
    csv_file = None
    csv_filename = None
    for filename in filenames:
        if os.path.isfile(filename):
            csv_file = open(filename, 'rb')
            csv_filename = filename
            break
    if not csv_file:
        raise ImportError('UnifiedNN._load_csv: file %s not found in the system.' \
                          % os.path.basename(filenames[0]))
    #try to load csv
    try:
        csv_reader = list(csv.reader(csv_file, delimiter='\t', quotechar='"'))
    except:
        print 'UnifiedNN._load_csv: exception occured while trying to load %s' \
        % csv_filename
        raise
    #read csv and fill the dict
    table_dict = dict()
    #dictionary of column names
    row_dict   = dict()
    for ci in xrange(1,len(csv_reader[0])):
        row_dict[ci] = csv_reader[0][ci]
    #read row by row
    for row in csv_reader[1:]: #skip the first row
        table_dict[row[0]] = dict()
        for ci in row_dict:
            table_dict[row[0]][row_dict[ci]] = float(row[ci])
    return table_dict
#end def


def _install_paths(possible_paths, filename):
    return tuple(path+name for path,name in zip(possible_paths, (filename,)*len(possible_paths)))


class UnifiedNN(object):
    '''
    Load thermodynamic data of Unified Nearest Neighbour model
    and provide methods for convenient access to it.
    '''
    #class initialization flag
    _inited = False
    
    #paths to the thermodynamic tables
    _possible_paths = ('./', '/usr/local/share/degen_primer/', '/usr/share/degen_primer/')
     
    internal_NN_filename             = 'internal-NN.csv'
    _internal_NN_paths               = _install_paths(_possible_paths, internal_NN_filename)
    
    loops_filename                   = 'loops.csv'
    _loops_paths                     = _install_paths(_possible_paths, loops_filename)
    
    tri_tetra_hairpin_loops_filename = '3-4-hairpin-loops.csv'
    _tri_tetra_hairpin_loops_paths   = _install_paths(_possible_paths, tri_tetra_hairpin_loops_filename)
    
    dangling_ends_filename           = 'dangling-ends.csv'
    _dangling_ends_paths             = _install_paths(_possible_paths, dangling_ends_filename)
    ###############################################################################
    
    
    #constants
    R                       =  1.9872  #Universal gas constant cal/(K*mol)
    K0                      = -273.15  #Absolute temperature zero in degree Celsius
    K37                     =  37 - K0 #37C
    dG_Na_coefficient_oligo = -0.114   #kcal/mol for oligomers with length =< 16
    dG_Na_coefficient_poly  = -0.175   #kcal/mol for longer polymers
    dS_Na_coefficient       =  0.368   #e.u.
    T_DMSO_coefficient      =  0.75
    Loop_coefficient        =  2.44
    #The NN stabilities at 37◦ C range from −1.23 to −0.21 kcal/mol 
    #for CG/GA and AC/TC, respectively
    Terminal_mismatch_mean  = (-1.23 + -0.21)/2
    
    
    #thermodynamic tables
    internal_NN   = None #table of thermodynamic values for internal nearest neighbour nucleotide pairs
    dangling_ends = None
    loops         = None
    tri_tetra_hairpin_loops = None
    
    
    #constructor
    def __init__(self):
        if self._inited: return
        self.load_tables()
    #end def
    
    
    #bool value of an instance
    def __nonzero__(self):
        return self._inited
    
    
    @classmethod    
    def load_tables(cls):
        '''load thermodynamic tables'''
        cls.tri_tetra_hairpin_loops = _load_csv(cls._tri_tetra_hairpin_loops_paths)
        cls.loops         = dict((int(k), v) for k, v 
                                 in _load_csv(cls._loops_paths).iteritems())
        cls.internal_NN   = _load_csv(cls._internal_NN_paths)
        for nn in cls.internal_NN.keys():
            cls.internal_NN[nn[::-1]] = cls.internal_NN[nn]
        cls.dangling_ends = _load_csv(cls._dangling_ends_paths)
        cls._inited       = True
    #end def
    
    #'standard' enthalpy, enthropy and Gibbs energy
    @classmethod
    def pair_dPar_37(cls, seq, rev):
        '''return 'standard' parameter (enthalpy, enthropy and Gibbs energy) 
        for the given dinucleotide duplex'''
        key = '%s/%s' % (seq, rev)
        if '-' in key: pars = cls.dangling_ends.get(key, None)
        else: pars = cls.internal_NN.get(key, None)
        if pars: return pars
        #if not found anywhere, try to unambiguate sequences
        seqs = unambiguous_sequences(seq)
        revs = unambiguous_sequences(rev)
        pars = dict(dH=0, dS=0, dG=0)
        if len(seqs) == 1 and len(revs) == 1:
            print 'UnifiedNN.dPar: sequence is not in the database: %s' % key
            return pars
        for s in seqs:
            for r in revs:
                p = cls.pair_dPar_37(s, r)
                for par in pars: pars[par] += p[par]
        numpars = float(len(seqs)*len(revs))
        for par in pars: pars[par] /= numpars
        #cache results
        cls.internal_NN[key] = pars
        cls.internal_NN[key[::-1]] = pars
        return pars
    #end def

    @classmethod
    def _extrapolate_loop_dG_37(cls, length, loop_type):
        if length < 30:
            exp_len = length
            while exp_len not in cls.loops: exp_len += 1
        else: exp_len = 30
        return cls.loops[exp_len][loop_type] + cls.Loop_coefficient * cls.R * cls.K37/1000 * log(float(length)/exp_len)
    #end def
    
    
    @classmethod
    def internal_loop_dG_37(cls, length):
        #check for loop length
        if length > 30: 
            print 'Warning: thermodynamic parameters for loops of length 30 and \
            more are extrapolated from experimental data and may be erroneous.'
        #if loop length is in the table
        try: loop_dG = cls.loops[length]['internal'] 
        except KeyError: #interpolate loop dG
            loop_dG = cls._extrapolate_loop_dG_37(length, 'internal')
        loop_dG += 2*cls.Terminal_mismatch_mean
        return loop_dG
    #end def
    
    
    @classmethod
    def hairpin_loop_dG_37(cls, seq):
        length = len(seq) - 2 #two boundary nucleotides are not included
        #check for loop length
        if length > 30: 
            print 'Warning: thermodynamic parameters for loops of length 30 and \
            more are extrapolated from experimental data and may be erroneous.'
        #if loop length is in the table
        if length in cls.loops:
            loop_dG = cls.loops[length]['hairpin']
            if   length == 3:
                if seq in cls.tri_tetra_hairpin_loops:
                    #special tri-loop correction
                    loop_dG += cls.tri_tetra_hairpin_loops[seq]['dG']
                if seq[0] in 'AT':
                    loop_dG += 0.5 #kcal/mol; AT-closing penalty
            elif length == 4:
                if seq in cls.tri_tetra_hairpin_loops:
                    #special tri-loop correction
                    loop_dG += cls.tri_tetra_hairpin_loops[seq]['dG']
                loop_dG += cls.Terminal_mismatch_mean
            else: loop_dG += cls.Terminal_mismatch_mean
            return loop_dG
        else: #interpolate loop dG
            loop_dG  = cls._extrapolate_loop_dG_37(length, 'hairpin')
            loop_dG += cls.Terminal_mismatch_mean
            return loop_dG
    #end def
    
    
    @classmethod
    def hairpin_loop_dH_37(cls, seq):
        if len(seq)-2 > 4: #two boundary nucleotides are not included 
            return 0 # All loop dH◦ parameters are assumed to equal zero. [1]
        try: return cls.tri_tetra_hairpin_loops[seq]['dH']
        except KeyError: return 0
    #end def
    
    @classmethod
    def internal_loop_dS_37(cls, seq):
        return cls.internal_loop_dG_37(len(seq)-2)*-1000.0/cls.K37
    
    @classmethod
    def hairpin_loop_dS_37(cls, seq):
        loop_dG = cls.hairpin_loop_dG_37(seq)
        loop_dH = cls.hairpin_loop_dH_37(seq)
        return (loop_dG-loop_dH)*-1000.0/cls.K37
    #end def
    
    
    #Gibbs energy at specified temperature
    @classmethod
    def dG_T(cls, dH, dS, temperature):
        return dH - (temperature-cls.K0)*dS/1000.0
#end class
    

#tests
if __name__ == '__main__':
    UnifiedNN.load_tables()
    print UnifiedNN.loops