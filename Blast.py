#!/usr/bin/python
# coding=utf-8

# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# indicator_gddccontrol is free software: you can redistribute it and/or modify it
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
Created on Jun 26, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from time import ctime
from Bio.Blast import NCBIWWW, NCBIXML
from StringTools import hr, text_width, wrap_text

def blast_short(query, blast_id='blast_short', entrez_query=''):
    #values adjusted for short sequences
    e_val  =  1000 #E-value 
    w_size =  7    #word size
    n_pen  = -2    #Reward and penalty for matching and mismatching bases
    n_rew  =  1    #Reward and penalty for matching and mismatching bases
    fltr   = 'none' #Turn off filtering of the results
    database = 'nt' #database to search: see  NCBI's Program Selection Guide for details
    #make blast query
    blast_results = NCBIWWW.qblast('blastn', database, query, 
                               expect=e_val, 
                               word_size=w_size,
                               nucl_penalty=n_pen,
                               nucl_reward=n_rew,
                               filter=fltr,
                               entrez_query=entrez_query)
    #save results to a file
    results_filename = blast_id+'-blast.xml'
    results_file = open(results_filename, 'w')
    results_file.write(blast_results.read())
    results_file.close()
    blast_results.close()
    print 'Blast output was written to: '+results_filename
    #parse results
    results_file  = open(results_filename, 'r')
    blast_results = list(NCBIXML.parse(results_file))
    results_file.close()
    return blast_results


def hits_histogram(blast_results, top_hits=10):
    histogram = dict()
    w = text_width
    for blast_record in blast_results:
        top_hits = top_hits if len(blast_record.alignments) > top_hits else len(blast_record.alignments)
        for a in range(top_hits):
            hit = blast_record.alignments[a].title.split('|')[-1]
            hit = hit[:w-16]+'...' if len(hit) > w-16 else hit
            if hit in histogram:
                histogram[hit] += 1
            else: histogram[hit] = 1
    num_records = len(blast_results)
    for hit in histogram:
        histogram[hit] *= 10.0/num_records
    return histogram
#end def


def format_histogram(histogram):
    histogram_string = ''
    max_col      = len(max(histogram, key=lambda x: len(x)))
    col_names    = ['hit description', ' % of records']
    names_spacer = max_col-len(col_names[0])
    histogram_string += '-'*(names_spacer/2)+col_names[0] + \
                        '-'*(names_spacer-names_spacer/2) + col_names[1] + '\n'
    for col in sorted(zip(histogram, histogram.values()), reverse=True, key=lambda x: x[1]):
        name_spacer = max_col - len(col[0])
        col_spacer  = 10-int(col[1])
        histogram_string += col[0]  + ' '*name_spacer + ' '
        histogram_string += ':' + '#'*int(col[1]) + ' '*col_spacer + ':' + '\n'
    histogram_string += '\n'
    histogram_string += wrap_text('The value of a histogram column is a percent of RECORDS '
    'in which there is an alignment of a query with the named sequence with a score '
    'within the "top-hits".\n')
    return histogram_string
#end def


def write_blast_report(blast_results, blast_report_filename, top_hits=10, top_hsps=4):
    blast_report = open(blast_report_filename, 'w')
    blast_report.write(hr(' %s ' % ctime(), symbol='#'))
    #header
    blast_report.write(wrap_text('If both sense and antisense primers are provided, '
                       'each primer from the sense group is blasted with every primer '
                       'from the antisense group using a compound query: senseNNNNNNNNNNantisense.\n'
                       'In that case each RECORD of these blast results corresponds '
                       'to a particular pair of primers. '
                       'If only one primer is provided, each RECORD corresponds to a '
                       'single unambiguous primer.'
                       '\n\n'))
    #hits histogram
    blast_report.write(hr(' hits histogram ', symbol='#'))
    histogram = hits_histogram(blast_results)
    blast_report.write(format_histogram(histogram))
    blast_report.write('\n')
    for r in range(len(blast_results)):
        blast_record = blast_results[r]
        blast_report.write('\n')
        blast_report.write(hr(' RECORD #%d: %s ' % (r+1, blast_record.query), symbol='#'))
        blast_report.write('Query length: %s\n'   % blast_record.query_length)
        blast_report.write('Hits:         %s\n\n' % len(blast_record.alignments))
        #print top 10? hits
        top_hits = top_hits if len(blast_record.alignments) > top_hits else len(blast_record.alignments) 
        blast_report.write(hr(' top %d hits ' % top_hits,  symbol='#'))
        for a in range(top_hits):
            alignment = blast_record.alignments[a]
            blast_report.write(hr(' Hit #%d ' % (a+1), symbol='='))
            blast_report.write(alignment.title+'\n')
            blast_report.write('Length:     %d\n' % alignment.length)
            blast_report.write('Max bits:   %d\n' % blast_record.descriptions[a].bits)
            blast_report.write('Alignments: %d\n' % len(alignment.hsps))
            #print to 4? hsps
            top_hsps = top_hsps if len(alignment.hsps) > top_hsps else len(alignment.hsps)  
            blast_report.write(hr(' top %01d alignments ' % top_hsps))
            for h in range(top_hsps):
                hsp = alignment.hsps[h]
                blast_report.write('\n'+str(hsp)+'\n')
                blast_report.write('Strand: ')
                if hsp.frame[0] == 1:
                    blast_report.write('Plus/')
                elif hsp.frame[0] == -1:
                    blast_report.write('Minus/')
                if hsp.frame[1] == 1:
                    blast_report.write('Plus\n')
                elif hsp.frame[1] == -1:
                    blast_report.write('Minus\n')
                blast_report.write('\n')
        blast_report.write(hr('', symbol='#'))
    blast_report.close()
    print 'Blast report was written to: '+blast_report_filename