#!/usr/bin/python
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
Created on Jun 26, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import os
from ConfigParser import SafeConfigParser
from StringTools import hr, wrap_text, time_hr, print_exception
from SecStructures import Dimer, SecStructures
from TD_Functions import dimer_dG, format_PCR_conditions
from PCR_Results import PCR_Results
try:
    from Bio.Blast import NCBIWWW, NCBIXML 
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
except Exception, e:
    print_exception(e)
    raise ImportError('The BioPython must be installed in your system.')


class Blast(object):
    '''Wrapper for Biopython NCBI-BLAST and parser for it's results'''

    #values adjusted for short sequences
    e_val    =  1000  #E-value 
    w_size   =  7     #word size
    n_pen    = -2     #Reward and penalty for matching and mismatching bases
    n_rew    =  1     #Reward and penalty for matching and mismatching bases
    fltr     = 'none' #Turn off filtering of the results
    database = 'nt'   #database to search: see  NCBI's Program Selection Guide for details
    spacer   = 'N'*20 #single query separator for concatenation

    #histogram column titles
    _col_titles = ('hit description', '% of records')
    _products_col_titles = ('<target name>', 'relative concentration')
    _all_products_col_titles = ('PCR-product', 'relative concentration')
    
    
    def __init__(self, job_id, min_amplicon, max_amplicon, with_exonuclease):
        self._job_id           = job_id
        self._blast_results    = None
        self._boundaries       = None
        self._query            = None
        self._min_amplicon     = min_amplicon
        self._max_amplicon     = max_amplicon
        self._with_exonuclease = with_exonuclease
        self._PCR_products     = None
        self._have_hits_report = False
        self._have_PCR_report  = False
        self._have_results     = False
    #end def

    #property functions
    def have_results(self): return self._have_results
    
    
    def load_results(self):
        results_filename    = self._job_id+'-blast.xml'
        query_filename      = self._job_id+'-blast.cfg'
        #if there's no saved results, quit silently
        if not os.path.isfile(results_filename): return
        #load blast results
        print '\nLoading previously saved BLAST results:\n   ', results_filename
        try:
            results_file        = open(results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            #load query configuration
            query_config        = SafeConfigParser()
            query_config.read(query_filename)
            #parse query configuration
            self._query      = SeqRecord(Seq(query_config.get('query', 'query'), 
                                             IUPAC.unambiguous_dna), self._job_id)
            self._boundaries = eval(query_config.get('query', 'boundaries'))
            self._have_results = True
        except Exception, e:
            print '\nFailed to load blast results.'
    #end def
    

    def blast_short(self, query, entrez_query=''):
        if not query: return
        print '\nStarting a BLAST search. This may take awhile...'
        try:
            blast_results = NCBIWWW.qblast('blastn', 
                                           self.database, 
                                           self._query.format('fasta'), 
                                           expect       = self.e_val, 
                                           word_size    = self.w_size,
                                           nucl_penalty = self.n_pen,
                                           nucl_reward  = self.n_rew,
                                           filter       = self.fltr,
                                           entrez_query = entrez_query)
            #save results to a file
            results_filename = self._job_id+'-blast.xml'
            results_file = open(results_filename, 'w')
            results_file.write(blast_results.read())
            results_file.close()
            blast_results.close()
            print '\nBlast output was written to:\n   '+results_filename
            #parse results
            results_file  = open(results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
        except Exception, e:
            print '\nFailed to obtain BLAST search results from NCBI.'
            print_exception(e)
            return
        self._have_results = True
    #end def


    def blast_primers(self, all_primers, entrez_query=''):
        #construct a concatenated query: primer1NNNNNNNprimer2NNNNNNprimer3...
        #also save primer boundaries positions in the concatenate
        if not all_primers: return
        self._boundaries = [[1,1]]
        self._query      = ''
        for p in range(len(all_primers)):
            primer       = all_primers[p]
            #construct query
            self._query += str(primer.seq) 
            #calculate boundaries
            self._boundaries[-1][1] += len(primer)-1
            #if not the last primer...
            if p < len(all_primers)-1:
                self._query += self.spacer
                next_start   = self._boundaries[-1][1] + len(self.spacer)+1
                self._boundaries.append([next_start, next_start])
        self._query = SeqRecord(Seq(self._query, IUPAC.unambiguous_dna), self._job_id)
        #save query config
        query_config   = SafeConfigParser()
        query_config.optionxform = str
        query_config.add_section('query')
        query_config.set('query', 'query', str(self._query.seq))
        query_config.set('query', 'boundaries', str(self._boundaries))
        query_filename = self._job_id+'-blast.cfg'
        query_file     = open(query_filename, 'wb')
        query_config.write(query_file)
        query_file.close()
        print '\nBlast query configuration was written to:\n   '+query_filename
        #do the blast
        self.blast_short(self._query.format('fasta'), entrez_query)
    #end def
    
    
    def _check_hsp(self, hsp):
        #check dG constraint
        dimer  = Dimer()
        query  = Seq(hsp.query, IUPAC.unambiguous_dna) #5'-3'
        target = Seq(hsp.sbjct, IUPAC.unambiguous_dna).reverse_complement() #revcomp(3'-5')
        for i in range(len(hsp.match)):
            if hsp.match[i] == '|':
                dimer.add(i,i)
        dG    = dimer_dG(dimer, query, target)
        if dG > SecStructures.max_dimer_dG: return False
        #check 3' constraint
        query_start  = hsp.query_start
        query_end    = hsp.query_end
        for primer in self._boundaries:
            if query_start < primer[0] \
            or query_end  > primer[1]:
                continue
            if not self._with_exonuclease and query_end < primer[1]:
                return False
        return True
    #end def 
        
    
    def find_PCR_products(self, side_reactions=None, side_concentrations=None):
        self._PCR_products     = dict()
        sorted_hits = []
        for record in self._blast_results:
            self._PCR_products[record.query] = PCR_Results(self._min_amplicon, self._max_amplicon, self._with_exonuclease) 
            if side_reactions: self._PCR_products[record.query].add_reactions(side_reactions, side_concentrations)
            sorted_hits.append(dict())
            record_hits = sorted_hits[-1]
            for alignment in record.alignments:
                hit_title = (alignment.title.split('|')[-1]).strip()
                record_hits[hit_title] = {'fwd':list(), 'rev':list()}
                hits = record_hits[hit_title]
                #check and sort all hsps 
                for hsp in alignment.hsps:
                    #add hsp to corresponding dictionary
                    target_end   = hsp.sbjct_end
                    primer_dir    = hsp.frame[0]
                    target_dir   = hsp.frame[1]
                    primer       = Seq(hsp.query, IUPAC.unambiguous_dna) #5'->3'
                    target       = Seq(hsp.sbjct, IUPAC.unambiguous_dna).reverse_complement() #5'->3'
                    #fwd hsp
                    if primer_dir == 1 and target_dir == 1:
                        hits['fwd'].append({'start':target_end,
                                            'fwd_seq':target,
                                            'fwd_primer':primer})
                    #rev hsp
                    elif primer_dir == 1 and target_dir == -1:
                        hits['rev'].append({'end':target_end,
                                            'rev_seq':target,
                                            'rev_primer':primer})
                #search for possible PCR products
                for fwd_hit in hits['fwd']:
                    for rev_hit in hits['rev']:
                        if fwd_hit['start'] < rev_hit['end']:
                            self._PCR_products[record.query].add_product(
                                                    hit_title,
                                                    fwd_hit['start'], 
                                                    rev_hit['end'], 
                                                    rev_hit['end'] - fwd_hit['start'],
                                                    fwd_hit['fwd_seq'],
                                                    rev_hit['rev_seq'],
                                                    fwd_hit['fwd_primer'],
                                                    rev_hit['rev_primer']) 
            #compute PCR products quantities
            self._PCR_products[record.query].compute_quantities()
            #remove empty query products dict
            if not self._PCR_products[record.query]: 
                del self._PCR_products[record.query]
        if not self._PCR_products:
            print '\nNo possible PCR products have been found using BLAST hits.'
    #end def
    

    def write_PCR_report(self):
        if not self._blast_results: return
        if not self._PCR_products:  return
        self._blast_PCR_report_filename = self._job_id + '-blast-PCR-report.txt'
        blast_report = open(self._blast_PCR_report_filename, 'w')
        blast_report.write(time_hr())
        #header
        blast_report.write(wrap_text('All hits are sorted into "forward" and ' 
                                     '"reverse" groups. Pairs of forward and reverse ' 
                                     'hits comprise possible PCR products which are ' 
                                     'filtered by amplicon size.\n'
                                     'If --no-exonuclease option was '
                                     'provided, hits with mismatches on ' 
                                     "3'-end are ignored.\n"
                                     'Relative quantities of the remaining products '
                                     'are estimated using equilibrium equations and ' 
                                     'current PCR parameters.\n'))
        blast_report.write('\n')
        #filter parameters
        blast_report.write(hr(' filtration parameters '))
        if self._with_exonuclease:
            blast_report.write("DNA polymerase HAS 3'-5'-exonuclease activity\n")
        else: blast_report.write("DNA polymerase doesn't have 3'-5'-exonuclease activity\n")
        blast_report.write('Minimum amplicon size:      %d\n' % self._min_amplicon)
        blast_report.write('Maximum amplicon size:      %d\n' % self._max_amplicon)
        blast_report.write('\n')
        blast_report.write(hr(' PCR conditions '))
        blast_report.write(format_PCR_conditions()+'\n')
        blast_report.write('\n\n\n')
        for record_name in self._PCR_products:
            blast_report.write(hr(' query ID: %s ' % record_name, symbol='#'))
            blast_report.write(self._PCR_products[record_name].format_quantity_explanation())
            #all products histogram
            blast_report.write(hr(' histogram of all possible PCR products ', symbol='='))
            blast_report.write(self._PCR_products[record_name].all_products_histogram())
            blast_report.write('\n\n\n')
            #all products electrophoresis
            blast_report.write(hr(' electrophorogram of all possible PCR products ', symbol='='))
            blast_report.write(self._PCR_products[record_name].all_products_electrophoresis())
            blast_report.write('\n\n')
            #PCR products by hit, if there more than one hit 
            if len(self._PCR_products[record_name].hits()) > 1:
                blast_report.write(hr(' the same grouped by target sequence ', symbol='='))
                blast_report.write(self._PCR_products[record_name].all_graphs_grouped_by_hit())
        blast_report.close()
        print '\nPossible PCR products defined by hits from BLAST search were written to:\n   ' + \
            self._blast_PCR_report_filename
        self._have_PCR_report = True
    #end def

    
    def write_hits_report(self):
        if not self._blast_results: return
        self._blast_report_filename = self._job_id + '-blast-hits-report.txt'
        blast_report = open(self._blast_report_filename, 'w')
        blast_report.write(time_hr())
        #header
        blast_report.write(wrap_text('All hits are filtered by number of mismatches '
                                     'and, if --no-exonuclease option was '
                                     'provided, hits with mismatches on '
                                     "3'-end are also filtered.\n"))
        blast_report.write('\n')
        #filter parameters
        blast_report.write(hr(' filtration parameters '))
        if self._with_exonuclease:
            blast_report.write("DNA polymerase HAS 3'-5'-exonuclease activity\n")
        else: blast_report.write("DNA polymerase doesn't have 3'-5'-exonuclease activity\n")
        blast_report.write('Maximum dG of an alignment: %.2f\n' % SecStructures.max_dimer_dG)
        blast_report.write('\n')
        blast_report.write(hr(''))
        blast_report.write('\n\n')
        #print records
        for r in range(len(self._blast_results)):
            blast_record = self._blast_results[r]
            num_hits = len(blast_record.alignments)
            blast_report.write(hr(' query ID: %s '  % blast_record.query, symbol='#'))
            blast_report.write(hr(' %d hits ' % num_hits,  symbol='#'))
            for a in range(num_hits):
                alignment = blast_record.alignments[a]
                blast_report.write(hr(' Hit #%d ' % (a+1), symbol='='))
                blast_report.write(alignment.title+'\n')
                blast_report.write('Length:     %d\n' % alignment.length)
                blast_report.write('Max bits:   %d\n' % blast_record.descriptions[a].bits)
                blast_report.write('Alignments: %d\n' % len(alignment.hsps))
                blast_report.write('\n')
                #print hsps
                hsps_text    = ''
                printed_hsps = 0
                for hsp in alignment.hsps:
                    if not self._check_hsp(hsp):
                        continue
                    printed_hsps += 1
                    hsps_text += str(hsp)+'\n'
                    hsps_text += 'Strand: '
                    if hsp.frame[0] == 1:
                        hsps_text += 'Plus/'
                    elif hsp.frame[0] == -1:
                        hsps_text += 'Minus/'
                    if hsp.frame[1] == 1:
                        hsps_text += 'Plus\n'
                    elif hsp.frame[1] == -1:
                        hsps_text += 'Minus\n'
                    hsps_text += '\n\n'
                hsps_text = hr(' %d alignments after filtration ' % printed_hsps) + hsps_text[:-1]
                blast_report.write(hsps_text)
                blast_report.write(hr('', symbol='='))
                if a < num_hits-1: blast_report.write('\n\n')
            blast_report.write(hr('', symbol='#'))
            if r < len(self._blast_results)-1: blast_report.write('\n\n')
        blast_report.close()
        print '\nTop hits with top HSPs from BLAST results were written to:\n   ' + \
            self._blast_report_filename
        self._have_hits_report = True
    #end def
    
    
    def reports(self):
        reports = []
        if self._have_hits_report: 
            reports.append({'report_name': 'BLAST hits', 'report_file': self._blast_report_filename})
        if self._have_PCR_report:
            reports.append({'report_name': 'BLAST PCR', 'report_file': self._blast_PCR_report_filename})
        if reports: return reports
        else: return None
#end class