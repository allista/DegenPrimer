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
Created on Jun 25, 2012

@author: Allis Tauri <allista@gmail.com>

All calculations and data are based on:
﻿SantaLucia, J., & Hicks, D. (2004). 
The thermodynamics of DNA structural motifs. Annual review of biophysics and 
biomolecular structure, 33, 415-40. doi:10.1146/annurev.biophys.32.110601.141800
'''

Tri_Tetra_Loops = {
                   #dG(37C) kcal/mol for tri- and tetra-loops
                   #triloops
                   'AGAAT': -1.5,
                   'AGCAT': -1.5,
                   'AGGAT': -1.5,
                   'AGTAT': -1.5,
                   'CGAAG': -2,
                   'CGCAG': -2,
                   'CGGAG': -2,
                   'CGTAG': -2,
                   'GGAAC': -2,
                   'GGCAC': -2,
                   'GGGAC': -2,
                   'GGTAC': -2,
                   'TGAAA': -1.5,
                   'TGCAA': -1.5,
                   'TGGAA': -1.5,
                   'TGTAA': -1.5,

                   #tetraloops:
                   'AAAAAT':  0.7,
                   'AAAACT':  0.2,
                   'AAACAT':  0.5,
                   'ACTTGT': -1.3,
                   'AGAAAT': -1.6,
                   'AGAGAT': -1.6,
                   'AGATAT': -2,
                   'AGCAAT': -2.1,
                   'AGCGAT': -1.6,
                   'AGCTTT': -0.3,
                   'AGGAAT': -1.6,
                   'AGGGAT': -1.6,
                   'AGGGGT':  0.3,
                   'AGTAAT': -2.1,
                   'AGTGAT': -1.6,
                   'AGTTCT':  0.3,
                   'ATTCGT': -0.7,
                   'ATTTGT': -0.5,
                   'ATTTTT': -1,
                   'CAAAAG':  0.9,
                   'CAAACG':  0.7,
                   'CAACAG':  1,
                   'CAACCG':  0,
                   'CCTTGG': -0.8,
                   'CGAAAG': -1.1,
                   'CGAGAG': -1.1,
                   'CGATAG': -1.5,
                   'CGCAAG': -1.6,
                   'CGCGAG': -1.1,
                   'CGCTTG':  0.2,
                   'CGGAAG': -1.1,
                   'CGGGAG': -1,
                   'CGGGGG':  0.8,
                   'CGTAAG': -1.6,
                   'CGTGAG': -1.1,
                   'CGTTCG':  0.8,
                   'CTTCGG': -0.2,
                   'CTTTGG':  0,
                   'CTTTTG': -0.5,
                   'GAAAAC':  1.5,
                   'GAAACC':  0.7,
                   'GAACAC':  1,
                   'GCTTGC': -0.8,
                   'GGAAAC': -1.1,
                   'GGAGAC': -1.1,
                   'GGATAC': -1.6,
                   'GGCAAC': -1.6,
                   'GGCGAC': -1.1,
                   'GGCTTC':  0.2,
                   'GGGAAC': -1.1,
                   'GGGGAC': -1.1,
                   'GGGGGC':  0.8,
                   'GGTAAC': -1.6,
                   'GGTGAC': -1.1,
                   'GGTTCC':  0.8,
                   'GTTCGC': -0.2,
                   'GTTTGC':  0,
                   'GTTTTC': -0.5,
                   'GAAAAT':  1.5,
                   'GAAACT':  1,
                   'GAACAT':  1,
                   'GCTTGT': -0.5,
                   'GGAAAT': -1.1,
                   'GGAGAT': -1.1,
                   'GGATAT': -1.6,
                   'GGCAAT': -1.6,
                   'GGCGAT': -1.1,
                   'GGCTTT': -0.1,
                   'GGGAAT': -1.1,
                   'GGGGAT': -1.1,
                   'GGGGGT':  0.8,
                   'GGTAAT': -1.6,
                   'GGTGAT': -1.1,
                   'GTATAT': -0.5,
                   'GTTCGT': -0.4,
                   'GTTTGT': -0.4,
                   'GTTTTT': -0.5,
                   'TAAAAA':  0.4,
                   'TAAACA':  0.2,
                   'TAACAA':  0.5,
                   'TCTTGA': -1.3,
                   'TGAAAA': -1.6,
                   'TGAGAA': -1.6,
                   'TGATAA': -2.1,
                   'TGCAAA': -2.1,
                   'TGCGAA': -1.6,
                   'TGCTTA': -0.3,
                   'TGGAAA': -1.6,
                   'TGGGAA': -1.6,
                   'TGGGGA':  0.3,
                   'TGTAAA': -2.1,
                   'TGTGAA': -1.6,
                   'TGTTCA':  0.3,
                   'TTTCGA': -0.7,
                   'TTTTGA': -0.5,
                   'TTTTTA': -1,
                   'TAAAAG':  1,
                   'TAAACG':  0.5,
                   'TAACAG':  0.5,
                   'TCTTGG': -1,
                   'TGAAAG': -1.5,
                   'TGAGAG': -1.5,
                   'TGATAG': -2,
                   'TGCAAG': -2,
                   'TGCGAG': -1.5,
                   'TGCTTG': -0.6,
                   'TGGAAG': -1.5,
                   'TGGGAG': -1.5,
                   'TGGGGG':  0.3,
                   'TGTAAG': -2,
                   'TGTGAG': -1.5,
                   'TTTCGG': -0.9,
                   'TTTTAG': -1.5,
                   'TTTTGG': -0.9,
                   'TTTTTG': -1,
}