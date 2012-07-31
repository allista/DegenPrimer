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

from itertools import chain
from ConfigParser import SafeConfigParser
import StringTools
from StringTools import hr, wrap_text, format_histogram, time_hr, print_exception
from Electrophoresis import print_electrophoresis
from TD_Functions import dimer_dG
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
    
    
    def __init__(self, job_id):
        self._job_id        = job_id
        self._blast_results = None
        self._boundaries    = None
        self._query         = None
        self._PCR_products  = None
        self.have_results   = False
    #end def

    
    def load_results(self):
        results_filename    = self._job_id+'-blast.xml'
        query_filename      = self._job_id+'-blast.cfg'
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
            self.have_results = True
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
        self.have_results = True
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
    
    
    def _check_hsp(self, hsp, max_dG, no_exonuclease):
        #check dG constraint
        dimer  = (set(), set(), [])
        query  = Seq(hsp.query, IUPAC.unambiguous_dna) #5'-3'
        target = Seq(hsp.sbjct, IUPAC.unambiguous_dna).reverse_complement() #revcomp(3'-5')
        for i in range(len(hsp.match)):
            if hsp.match[i] == '|':
                dimer[0].add(i)
                dimer[1].add(i)
        dG    = dimer_dG(dimer, query, target)
        if dG > max_dG: return False
        #check 3' constraint
        query_start  = hsp.query_start
        query_end    = hsp.query_end
        for primer in self._boundaries:
            if query_start < primer[0] \
            or query_end  > primer[1]:
                continue
            if no_exonuclease and query_end < primer[1]:
                return False
        return True
    #end def 
        
    
    def _find_PCR_products(self, min_amplicon, max_amplicon, max_dG, no_exonuclease):
        self._PCR_products = dict()
        sorted_hits = []
        for record in self._blast_results:
            self._PCR_products[record.query] = dict()
            sorted_hits.append(dict())
            record_hits = sorted_hits[-1]
            for alignment in record.alignments:
                hit_title = (alignment.title.split('|')[-1]).strip()
                record_hits[hit_title] = {'fwd':dict(), 'rev':dict()}
                hits = record_hits[hit_title]
                #check and sort all hsps 
                for hsp in alignment.hsps:
                    #check for 3' matching and mismatch number
                    if not self._check_hsp(hsp, max_dG, no_exonuclease):
                        continue
                    #add hsp to corresponding dictionary
                    target_end   = hsp.sbjct_end
                    query_dir    = hsp.frame[0]
                    target_dir   = hsp.frame[1]
                    #fwd hsp
                    if query_dir == 1 and target_dir == 1:
                        if target_end not in hits['fwd']:
                            hits['fwd'][target_end] = 1
                        else: hits['fwd'][target_end] += 1
                    #rev hsp
                    elif query_dir == 1 and target_dir == -1:
                        if target_end not in hits['rev']:
                            hits['rev'][target_end] = 1
                        else: hits['rev'][target_end] += 1
                #search for possible PCR products
                self._PCR_products[record.query][hit_title] = []
                products = self._PCR_products[record.query][hit_title]
                for fwd_hit in hits['fwd']:
                    for rev_hit in hits['rev']:
                        if fwd_hit < rev_hit \
                        and rev_hit-fwd_hit < max_amplicon \
                        and rev_hit-fwd_hit < max_amplicon:
                            products.append({'target'  :hit_title,
                                             'start'   :fwd_hit, 
                                             'end'     :rev_hit, 
                                             'len'     :(rev_hit-fwd_hit), 
                                             'quantity':min(hits['fwd'][fwd_hit],
                                                            hits['rev'][rev_hit])})
                #remove empty target products dict
                if not products: del self._PCR_products[record.query][hit_title]
            #remove empty query products dict
            if not self._PCR_products[record.query]: del self._PCR_products[record.query]
    #end def

    
    def _products_histogram(self, name, products):
        text_width = StringTools.text_width
        hist_width = len(self._products_col_titles[1])
        #construct histogram
        histogram  = dict()
        max_start  = max(len(str(p['start'])) for p in products)
        max_end    = max(len(str(p['end'])) for p in products)
        max_len    = max(len(str(p['len'])) for p in products)
        prod_limit = text_width-hist_width-2
        for product in products:
            product_spec   = '%d%s bp [%d%s-%s%d]' % (product['len'],
                                                      ' '*(max_len-len(str(product['len']))),
                                                      product['start'],
                                                      ' '*(max_start-len(str(product['start']))),
                                                      ' '*(max_end-len(str(product['end']))),
                                                      product['end'])
            
            colname = product_spec+' '*(prod_limit - len(product_spec))
            if colname not in histogram:
                histogram[colname] = [product['quantity'],product['start']]
            else: histogram[colname][0] += product['quantity'] 
        #normalize histogram
        norm = float(max(histogram.values(), key=lambda x: x[0])[0])
        for colname in histogram:
            histogram[colname][0] /= norm
        #sort histogram by product len
        histogram = sorted(zip(histogram, histogram.values()), key=lambda x: x[1][1])
        histogram = tuple((line[0], line[1][0]) for line in histogram)
        #format histogram
        return format_histogram(histogram, (name, self._products_col_titles[1]), hist_width)
    #end def

    
    def _all_products_histogram(self, products):
        all_products = list(chain(*(products.values())))
        text_width = StringTools.text_width
        hist_width = len(self._all_products_col_titles[1])
        #construct histogram
        histogram  = dict()
        max_target = max(len(p['target']) for p in all_products)
        max_len    = max(len(str(p['len'])) for p in all_products)
        for product in all_products:
            target_spec   = ' %d%s bp [%d-%d]' % (product['len'],
                                                 ' '*(max_len-len(str(product['len']))),
                                                 product['start'],
                                                 product['end'])
            target_limit  = text_width-hist_width-2 - len(target_spec)
            target_spacer = max_target - len(product['target'])
            if target_limit < 0: target_limit = 0
            colname = (product['target']+' '*target_spacer)[:target_limit] + target_spec
            if colname not in histogram:
                histogram[colname] = [product['quantity'],product['len']]
            else: histogram[colname][0] += product['quantity'] 
        #normalize histogram
        norm = float(max(histogram.values(), key=lambda x: x[0])[0])
        for colname in histogram:
            histogram[colname][0] /= norm
        #sort histogram by product len
        histogram = sorted(zip(histogram, histogram.values()), key=lambda x: x[1][1])
        histogram = tuple((line[0], line[1][0]) for line in histogram)
        #format histogram
        return format_histogram(histogram, self._all_products_col_titles, hist_width)
    #end def
    

    def write_PCR_report(self, min_amplicon, max_amplicon, max_dG, no_exonuclease):
        if not self._blast_results: return
        #parse results
        self._find_PCR_products(min_amplicon, max_amplicon, max_dG, no_exonuclease)
        #
        self._blast_PCR_report_filename = self._job_id + '-blast-PCR-report.txt'
        blast_report = open(self._blast_PCR_report_filename, 'w')
        blast_report.write(time_hr())
        #header
        blast_report.write(wrap_text('All hits are filtered by number of mismatches '
                                     'and, if --no-exonuclease option was '
                                     'provided, hits with mismatches on '
                                     "3'-end are also filtered.\n"
                                     'Then hits are sorted into "forward" and '
                                     '"reverse" groups. Pairs of forward and reverse '
                                     'hits comprise possible PCR products which are '
                                     'in their turn filtered by amplicon size.\n'))
        blast_report.write('\n')
        #filter parameters
        blast_report.write(hr(' filtration parameters '))
        if no_exonuclease:
            blast_report.write("DNA polymerase DOES NOT have 3'-5'-exonuclease activity\n")
        else: blast_report.write("DNA polymerase has 3'-5'-exonuclease activity\n")
        blast_report.write('Minimum amplicon size:      %d\n' % min_amplicon)
        blast_report.write('Maximum amplicon size:      %d\n' % max_amplicon)
        blast_report.write('Maximum dG of an alignment: %.2f\n' % max_dG)
        blast_report.write('\n\n\n')
        #if no PCR products have been found
        if not self._PCR_products:
            blast_report.write(hr(' No PCR products have been found ', symbol='!'))
            blast_report.close()
            return
        #else...
        for record_name in self._PCR_products:
            blast_report.write(hr(' query ID: %s ' % record_name, symbol='#'))
            #all products histogram
            blast_report.write(hr(' histogram of all possible PCR products ', symbol='='))
            blast_report.write(self._all_products_histogram(self._PCR_products[record_name]))
            blast_report.write(wrap_text('"relative concentration" of a product is a normalized '
                                         'number of primer pairs that yield this particular product.\n'))
            blast_report.write('\n\n\n')
            #all products electrophoresis
            blast_report.write(hr(' electrophorogram of all possible PCR products ', symbol='='))
            blast_report.write(print_electrophoresis(list(chain(*(self._PCR_products[record_name].values())))))
            blast_report.write('\n\n')
            #PCR products by target
            blast_report.write(hr(' the same grouped by target sequence ', symbol='='))
            for target_name in self._PCR_products[record_name]:
                blast_report.write(self._products_histogram(target_name, self._PCR_products[record_name][target_name]))
                blast_report.write(print_electrophoresis(self._PCR_products[record_name][target_name]))
                blast_report.write('\n')
                blast_report.write(hr(''))
                blast_report.write('\n')
        blast_report.close()
        print '\nPossible PCR products defined by hits from BLAST search were written to:\n   ' + \
            self._blast_PCR_report_filename
    #end def

    
    def write_hits_report(self, max_dG, no_exonuclease):
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
        if no_exonuclease:
            blast_report.write("DNA polymerase DOES NOT have 3'-5'-exonuclease activity\n")
        else: blast_report.write("DNA polymerase has 3'-5'-exonuclease activity\n")
        blast_report.write('Maximum dG of an alignment: %.2f\n' % max_dG)
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
                    if not self._check_hsp(hsp, max_dG, no_exonuclease):
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
            if a < len(self._blast_results)-1: blast_report.write('\n\n')
        blast_report.close()
        print '\nTop hits with top HSPs from BLAST results were written to:\n   ' + \
            self._blast_report_filename
    #end def
    
    
    def register_reports(self, args):
        args.register_report('BLAST hits', self._blast_report_filename)
        args.register_report('BLAST PCR', self._blast_PCR_report_filename)
#end class