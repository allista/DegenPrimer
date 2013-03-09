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
try:
    from Bio.Blast import NCBIWWW, NCBIXML 
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
except ImportError:
    print'The BioPython must be installed in your system.'
    raise
from iPCR_Interface import iPCR_Interface
from ConfigParser import SafeConfigParser
from StringTools import hr, wrap_text, time_hr, print_exception
from SecStructures import SecStructures, Dimer, Duplex
import TD_Functions


class BlastPrimers(iPCR_Interface):
    '''Check specificity of primers using BLAST query and iPCR_Simulation'''

    #values adjusted for short sequences
    e_val    =  100000 #E-value grater than 1e5 results in a veeeery long search 
    w_size   =  7      #word size
    n_pen    = -4      #Reward and penalty for matching and mismatching bases
    n_rew    =  3      #Reward and penalty for matching and mismatching bases
    fltr     = 'none'  #Turn off filtering of the results
    database = 'nt'    #database to search: see  NCBI's Program Selection Guide for details
    no_gaps  =  True   #no gaps in primers!
    spacer   = 'N'*20  #single query separator for concatenation
    
    
    def __init__(self, job_id, *args, **kwargs):
        iPCR_Interface.__init__(self, *args, **kwargs) 
        self._job_id           = job_id
        self._blast_results    = None
        self._bounds           = None
        self._query            = None
        self._format_query()
        #PCR parameters
        self._PCR_Simulations  = dict()
        #results
        self._results_filename = self._job_id+'-blast.xml'
        self._query_filename   = self._job_id+'-blast.cfg'
        #reports
        self._PCR_report_filename  = self._job_id + '-blast-PCR-report.txt'
        self._hits_report_filename = self._job_id + '-blast-hits-report.txt'
        #flags
        self._have_blast_results = False
    #end def
   
    
    def _format_query(self):
        '''Construct a concatenated query: primer1NNNNNNNprimer2NNNNNNprimer3...
        Also save primer boundaries positions in the concatenate.'''
        all_primers = []
        for primer in self._primers: all_primers.extend(primer.sequences)
        self._bounds = [[[1,1],]]
        query = ''
        for p in range(len(all_primers)):
            #construct query
            query += str(all_primers[p])
            #calculate boundaries
            self._bounds[-1][0][1] += len(all_primers[p])-1
            self._bounds[-1].append(all_primers[p]) 
            #if not the last primer, add polyN spacer
            if p < len(all_primers)-1:
                query += self.spacer
                next_start = self._bounds[-1][0][1] + len(self.spacer)+1
                self._bounds.append([[next_start, next_start],])
        self._query = SeqRecord(Seq(query, IUPAC.ambiguous_dna), self._job_id)
    #end def
    
    
    def _save_query_config(self):
        #save query config
        query_config   = SafeConfigParser()
        query_config.optionxform = str
        query_config.add_section('query')
        query_config.set('query', 'query', str(self._query.seq))
        query_config.set('query', 'boundaries', str(self._bounds))
        query_filename = self._job_id+'-blast.cfg'
        try:
            query_file     = open(query_filename, 'wb')
            query_config.write(query_file)
            query_file.close()
        except IOError, e:
            print '\nUnable to write blast configuration file:\n   %s' % query_filename
            print e.message
        print '\nBlast query configuration was written to:\n   %s' % query_filename
    #end def
    
    
    def _duplex_from_hsp(self, hsp):
        query_start = hsp.query_start
        query_end   = hsp.query_end
        for bounds,primer in self._bounds:
            if query_start < bounds[0] \
            or query_end   > bounds[1]:
                continue
            dimer = Dimer()
            for i in range(len(hsp.match)):
                if hsp.match[i] == '|':
                    dimer.add(query_start-bounds[0]+i,i)
            template = Seq(hsp.sbjct, IUPAC.unambiguous_dna).reverse_complement()
            return Duplex(primer, template, dimer)
        return None
    #end def 
    
    
    def blast_query(self, entrez_query=''):
        self._save_query_config()
        try:
            blast_results = NCBIWWW.qblast('blastn', 
                                           self.database, 
                                           self._query.format('fasta'), 
                                           expect       = self.e_val, 
                                           word_size    = self.w_size,
                                           nucl_penalty = self.n_pen,
                                           nucl_reward  = self.n_rew,
                                           filter       = self.fltr,
                                           entrez_query = entrez_query,
                                           ungapped_alignment = self.no_gaps,)
            #save results to a file
            results_file = open(self._results_filename, 'w')
            results_file.write(blast_results.read())
            results_file.close()
            blast_results.close()
            print '\nBlast output was written to:\n   %s' % self._results_filename
            #parse results
            results_file  = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
        except Exception, e:
            print '\nFailed to obtain BLAST query results from NCBI.'
            print_exception(e)
            return False
        self._have_blast_results = (len(self._blast_results) > 0 
                              and max(len(record.alignments) 
                                      for record in self._blast_results) > 0)
        if not self._have_blast_results:
            print '\nNCBI provided no results for the BLAST query.'
        return self._have_blast_results
    #end def


    def load_results(self):
        #check for file existence
        if  not os.path.isfile(self._results_filename) \
        or  not os.path.isfile(self._query_filename): 
            return False
        #load blast results
        try:
            #load query configuration
            query_config = SafeConfigParser()
            query_config.read(self._query_filename)
            #parse query configuration
            query = SeqRecord(Seq(query_config.get('query', 'query'), 
                                  IUPAC.ambiguous_dna), self._job_id)
            #if saved query represents different set of primers, discard saved results
            if str(self._query.seq) != str(query.seq): 
                raise AssertionError('BlastPrimers.load_results: saved query '
                                     'differs from current query.')
            #load blast query results
            results_file = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            self._have_blast_results = True
        except Exception, e:
            print '\nFailed to load blast results:\n   %s\n   %s' \
                % (self._results_filename, self._query_filename)
            print_exception(e)
            return False
        print '\nPreviously saved BLAST results were loaded:\n   %s\n   %s' \
              % (self._results_filename, self._query_filename)
        return True
    #end def
    
    
    def simulate_PCR(self):
        if not self._have_blast_results: return False
        for record in self._blast_results:
            _PCR_Sim = self._PCR_Simulation_factory()
            self._PCR_Simulations[record.query] = _PCR_Sim
            for alignment in record.alignments:
                hit_title = (alignment.title.split('|')[-1]).strip()
                fwd_annealings = []
                rev_annealings = []
                #check and sort all hsps
                for hsp in alignment.hsps:
                    #construct annealing duplex
                    hsp_duplex = self._duplex_from_hsp(hsp)
                    #find id of the primer
                    hsp_id = ''
                    for primer in self._primers:
                        hsp_id = primer.find_id(hsp_duplex.fwd_seq)
                        if hsp_id: break
                    #add hsp_duplex to corresponding dictionary
                    primer_dir = hsp.frame[0]
                    target_dir = hsp.frame[1]
                    #fwd hsp
                    if primer_dir == 1 and target_dir == 1:
                        product_bound = hsp.sbjct_end + hsp_duplex.fwd_3_overhang
                        fwd_annealings.append((product_bound, [(hsp_duplex, hsp_id),]))
                    #rev hsp
                    elif primer_dir == 1 and target_dir == -1:
                        product_bound = hsp.sbjct_end - hsp_duplex.fwd_3_overhang
                        rev_annealings.append((product_bound, [(hsp_duplex, hsp_id),]))
                #add primers' annealings to PCR simulation
                _PCR_Sim.add_annealings(hit_title, fwd_annealings, rev_annealings)
            #compute PCR products quantities
            _PCR_Sim.run()
            #remove empty query products dict
            if not _PCR_Sim: 
                del self._PCR_Simulations[record.query]
        if not self._PCR_Simulations:
            print '\nNo possible PCR products have been found using BLAST hits.'
            return False
        self._have_results = True
        return True
    #end def
    
    
    def blast_and_analyze(self, entrez_query=''):
        return self.blast_query(entrez_query) \
        and    self.simulate_PCR()
    #end def
    
    
    def write_PCR_report(self):
        if not self._have_results: return
        blast_report = self._open_report('BLAST PCR', self._PCR_report_filename)
        #header
        blast_report.write(time_hr())
        blast_report.write(wrap_text('For each hit all alignments are sorted '
                                     'into "forward" and "reverse" groups. '
                                     'Pairs of forward and reverse ' 
                                     'alignments comprise possible PCR products.\n'
                                     'If --no-exonuclease option was '
                                     'provided, alignments with mismatches on ' 
                                     "3'-end are ignored.\n"
                                     'Quantities of products are estimated '
                                     'using equilibrium equations and ' 
                                     'current PCR parameters.\n'))
        blast_report.write('\n')
        for record_name in self._PCR_Simulations:
            blast_report.write(hr(' query ID: %s ' % record_name, symbol='#'))
            blast_report.write(self._PCR_Simulations[record_name].format_report_header())
            blast_report.write(self._PCR_Simulations[record_name].format_quantity_explanation())
            if len(self._PCR_Simulations[record_name].hits()) == 1:
                #all products histogram
                hit = self._PCR_Simulations[record_name].hits()[0]
                blast_report.write(hr(' histogram of all possible PCR products ', symbol='='))
                blast_report.write(self._PCR_Simulations[record_name].per_hit_header(hit))
                blast_report.write(self._PCR_Simulations[record_name].per_hit_histogram(hit))
                blast_report.write('\n')
                blast_report.write(hr(' electrophorogram of all possible PCR products ', symbol='='))
                blast_report.write(self._PCR_Simulations[record_name].per_hit_electrophoresis(hit))
            else:
                #all products histogram
                blast_report.write(hr(' histogram of all possible PCR products ', symbol='='))
                blast_report.write(self._PCR_Simulations[record_name].all_products_histogram())
                blast_report.write('\n\n\n')
                #per hit histogram and phoresis
                blast_report.write(hr(' histograms and electrophorograms of PCR products of each hit ', symbol='='))
                blast_report.write(self._PCR_Simulations[record_name].all_graphs_grouped_by_hit())
        blast_report.close()
        print '\nPossible PCR products defined by hits from BLAST search were written to:\n   ' + \
            self._PCR_report_filename
        self._add_report('BLAST PCR', self._PCR_report_filename)
    #end def
    
        
    def write_hits_report(self):
        if not self._have_blast_results: return
        try:
            blast_report = open(self._hits_report_filename, 'w')
        except IOError, e:
            print '\nFailed to open blast hits report file for writing:\n   %s' % self._hits_report_filename
            print e.message
            return
        #header
        blast_report.write(time_hr())
        blast_report.write(wrap_text('All hits are filtered by dG of the annealing '
                                     'of alignments and, if --no-exonuclease '
                                     'option was provided, hits with mismatches on '
                                     "3'-end are also filtered.\n"))
        blast_report.write('\n')
        #filter parameters
        blast_report.write(hr(' filtration parameters '))
        if self._with_exonuclease:
            blast_report.write("DNA polymerase HAS 3'-5'-exonuclease activity\n")
        else: blast_report.write("DNA polymerase doesn't have 3'-5'-exonuclease activity\n")
        blast_report.write('Maximum dG of an alignment: %.2f kcal/mol\n' % SecStructures.max_dimer_dG)
        blast_report.write('\n')
        blast_report.write(hr(''))
        blast_report.write('\n\n')
        #print records
        for r in range(len(self._blast_results)):
            blast_record = self._blast_results[r]
            num_hits = len(blast_record.alignments)
            blast_report.write(hr(' query ID: %s '  % blast_record.query, symbol='#'))
            #filter hits by alignments and format report text for each hit
            hits = []
            for h in range(num_hits):
                hit  = blast_record.alignments[h]
                desc = blast_record.descriptions[h] 
                #check and format hsps
                hsps = []
                for hsp in hit.hsps:
                    hsp_duplex = self._duplex_from_hsp(hsp)
                    #check 3' mismatch
                    if not self._with_exonuclease \
                    and hsp_duplex.fwd_3_mismatch: continue
                    #check annealing enegry
                    if hsp_duplex.dG > SecStructures.max_dimer_dG: continue
                    #get primer concentration
                    if self._fwd_primer.has_subsequence(hsp_duplex.fwd_seq): 
                        hsp_primer_concentration = self._fwd_primer.concentration
                    else: hsp_primer_concentration = self._rev_primer.concentration
                    #format hsps representation
                    hsp_text  = ('score: %d; bits: %d; E-value: %.2e;\n\n'
                                 % (hsp.score, hsp.bits, hsp.expect))
                    hsp_text += str(hsp_duplex)
                    hsp_text += 'Conversion degree = %.2f%%\n\n' \
                    % (TD_Functions.primer_DNA_conversion_degree(hsp_primer_concentration, hsp_duplex.K)*100)
                    hsp_text += 'Template strand: '
                    if hsp.frame[1] == 1:    
                        hsp_text += 'antisense\n'
                        hsp_text += 'Position on template: %d ==> %d\n' \
                              % (hsp.sbjct_start, hsp.sbjct_end)
                    elif hsp.frame[1] == -1: 
                        hsp_text += 'sense\n'
                        hsp_text += 'Position on template: %d <== %d\n' \
                              % (hsp.sbjct_start, hsp.sbjct_end)
                    hsp_text += '\n\n'
                    hsps.append((hsp_duplex.dG, hsp_text))
                #no need to include weak hits to the report
                if not hsps: continue 
                #sort hsps by minimum dG
                hsps.sort(key=lambda(x): x[0])
                #format hit
                hit_str  = wrap_text(hit.title)+'\n'
                hit_str += 'Length:     %d\n' % hit.length
                hit_str += 'Max bits:   %d\n' % desc.bits
                hit_str += 'Alignments: %d\n' % len(hit.hsps)
                hit_str += '\n'
                hit_str += hr(' %d alignments after filtration ' % len(hsps))
                hit_str += ''.join(_hsp[1] for _hsp in hsps)[:-1]
                hits.append((hsps[0][0], desc.title, desc.score, desc.bits, desc.e, len(hsps), hit_str))
            if not hits:
                blast_report.write(wrap_text('All hits were filtered out.\n'))
            #write report to the file
            else:
                #print the short list of all hits
                num_hits     = len(hits)
                num_hits_len = len(str(num_hits))
                #sort hits by minimum dG
                hits.sort(key=lambda(x): x[0])
                #print header
                blast_report.write(hr(' %d hits ' % num_hits,  symbol='#'))
                for h, hit in enumerate(hits):
                    spacer    = ' '*(num_hits_len-len(str(h)))
                    blast_report.write(wrap_text('%d.%s %s\n' 
                                                 % (h+1, spacer, 
                                                   (hit[1].split('|')[-1]).strip())))
                    blast_report.write(('%s  min dG: %.2f; score: %d; bits: %d; '
                                        'E-value: %.2e; alignments: %d\n\n')
                                       % (' '*num_hits_len, hit[0], hit[2], 
                                          hit[3], hit[4], hit[5]))
                blast_report.write('\n')
                #print each formatted hit
                for h, hit in enumerate(hits):
                    blast_report.write(hr(' Hit #%d ' % (h+1), symbol='='))
                    blast_report.write(hit[6])
                    blast_report.write(hr('', symbol='='))
                    if h < num_hits-1: blast_report.write('\n\n')
            blast_report.write(hr('', symbol='#'))
            if r < len(self._blast_results)-1: blast_report.write('\n\n')
        blast_report.close()
        print '\nTop hits with top HSPs from BLAST results were written to:\n   ' + \
            self._hits_report_filename
        self._add_report('BLAST hits', self._hits_report_filename)
    #end def
#end class