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
# degen_primer is distributed in the hope that it will be useful, but
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

import os, re
try:
    from Bio.Blast import NCBIWWW, NCBIXML 
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
except ImportError:
    print'The BioPython must be installed in your system.'
    raise
from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from iPCR_Interface import iPCR_Interface
from ConfigParser import SafeConfigParser
from StringTools import hr, wrap_text, time_hr, print_exception
from SecStructures import Dimer, Duplex, max_dimer_dG
from WorkCounter import WorkCounter
import TD_Functions as tdf


class BlastPrimers(iPCR_Interface, MultiprocessingBase):
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
    
    _PCR_report_suffix = 'BLAST-PCR'
    
    def __init__(self, abort_event, *args, **kwargs):
        iPCR_Interface.__init__(self, abort_event, *args, **kwargs)
        MultiprocessingBase.__init__(self, abort_event) 
        self._blast_results    = None
        self._bounds           = None
        self._query            = None
        #PCR parameters
        self._PCR_Simulations  = dict()
        #query_id use hash of primers instead of all-config hash as in job_id
        self._primers_hash = (hash(tuple(self._primers)) & 0xFFFFFFF)
        self._query_id  = re.split('_[0-9]+\Z', self._job_id)[0]
        self._query_id += '_%s' % str(self._primers_hash)
        #results
        self._results_filename = self._query_id+'-blast.xml'
        self._query_filename   = self._query_id+'-blast.cfg'
        #reports
        self._hits_report_filename = '%s-%s-hits.txt' % (self._job_id, 
                                                         self._PCR_report_suffix.rstrip('-PCR'))
        #flags
        self._have_blast_results = False
        self._have_saved_results = False
    #end def
   
   
    def _format_query(self):
        '''Construct a concatenated query: primer1NNNNNNNprimer2NNNNNNprimer3...
        Also save primer boundaries positions in the concatenate.'''
        all_primers = []
        for primer in self._primers: all_primers.extend(primer.str_sequences)
        self._bounds = [[[1,1],]]
        query = ''
        for p, primer in enumerate(all_primers):
            #construct query
            query += primer
            #calculate boundaries
            self._bounds[-1][0][1] += len(primer)-1
            self._bounds[-1].append(primer) 
            #if not the last primer, add polyN spacer
            if p < len(all_primers)-1:
                query += self.spacer
                next_start = self._bounds[-1][0][1] + len(self.spacer)+1
                self._bounds.append([[next_start, next_start],])
        self._query = SeqRecord(Seq(query, IUPAC.ambiguous_dna), self._query_id)
    #end def
    
    
    def _save_query_config(self):
        #save query config
        query_config = SafeConfigParser()
        query_config.optionxform = str
        query_config.add_section('query')
        query_config.set('query', 'id', str(self._primers_hash))
        query_config.set('query', 'query', str(self._query.seq))
        query_config.set('query', 'boundaries', str(self._bounds))
        try:
            query_file = open(self._query_filename, 'wb')
            query_config.write(query_file)
            query_file.close()
        except IOError, e:
            print '\nUnable to write BLAST configuration file:\n   %s' % self._query_filename
            print e.message
        print '\nBLAST query configuration was written to:\n   %s' % self._query_filename
    #end def
    
    
    def _duplex_from_hsp(self, hsp):
        query_start = hsp.query_start
        query_end   = hsp.query_end
        for bounds,primer in self._bounds:
            if query_start < bounds[0] \
            or query_end   > bounds[1]:
                continue
            for i in xrange(len(hsp.match)):
                if hsp.match[i] == '|':
                    dimer = Dimer.from_sequences(primer, hsp.sbjct, query_start-bounds[0])
                    return Duplex(primer, hsp.sbjct, dimer, revcomp=True)
        return None
    #end def
    
    
    def blast_query(self, counter, entrez_query=''):
        counter.set_work(5)
        self._format_query(); counter.count()
        self._save_query_config(); counter.count()
        try:
            print '\nLaunching BLAST query #%d...' % self._primers_hash
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
            counter.count()
            #save results to a file
            results_file = open(self._results_filename, 'w')
            results_file.write(blast_results.read())
            results_file.close()
            blast_results.close()
            print '\nBLAST output was written to:\n   %s' % self._results_filename
            counter.count()
            #parse results
            results_file  = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            counter.count()
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


    def have_saved_results(self):
        if  not os.path.isfile(self._results_filename) \
        or  not os.path.isfile(self._query_filename): 
            print '\nNo results for BLAST query #%d were found.' % self._primers_hash
            return False
        self._have_saved_results = True
        return True
    #end def
    

    def load_results(self, counter):
        '''One should call have_saved_results() prior to this method.'''
        if not self._have_saved_results: return False
        print(('\nLoading previously saved BLAST results '
               'for query #%d...') % self._primers_hash)
        try:
            #check that results to be loaded are the results for the same query
            query_id = None
            with open(self._results_filename, 'r') as results:
                query_id_re = re.compile('\<BlastOutput_query-def\>(.*)\<\/BlastOutput_query-def\>')
                for line in results:
                    match = query_id_re.search(line)
                    if match is None: continue
                    query_id = match.group(1).split()[0]
                    break
            if query_id is None:
                print(('Error: <BlastOutput_query-def> tag was not found in %s'
                       'The file is corrupted.')
                      % self._results_filename)
                return False
            if query_id != self._query_id:
                print(('Error: query ID in saved results:\n   %s\n'
                       'does not match current query ID:\n   %s')
                      % (query_id, self._query_id))
                return False
            #load query configuration
            query_config = SafeConfigParser()
            query_config.read(self._query_filename)
            #parse query configuration
            self._bounds = eval(query_config.get('query', 'boundaries'))
            #load blast query results
            results_file = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            self._have_blast_results = len(self._blast_results)
        except Exception, e:
            print '\nFailed to load blast results:\n   %s\n   %s' \
                % (self._results_filename, self._query_filename)
            print_exception(e)
            return False
        print '\nPreviously saved BLAST results were loaded:\n   %s\n   %s' \
              % (self._results_filename, self._query_filename)
        counter.done()
        return True
    #end def
    
    
    def _find_products_in_alignment(self, alignment, products_finder):
        hit_title = (alignment.title.split('|')[-1]).strip()
        fwd_annealings = []
        rev_annealings = []
        #check and sort all hsps
        for hsp in alignment.hsps:
            if self.aborted(): return None
            #construct annealing duplex and check if it's stable
            hsp_duplex = self._duplex_from_hsp(hsp)
            if not hsp_duplex: continue
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
        mixture = products_finder.create_PCR_mixture(WorkCounter(),
                                                     hit_title, 
                                                     fwd_annealings, 
                                                     rev_annealings)
        if mixture is None: return hit_title, None
        return hit_title, mixture.save()
    #end def
    
    
    def simulate_PCR(self, counter):
        if not self._have_blast_results: return False
        print '\nSearching for possible PCR products in BLAST results...'
        P_Finder = self._new_PCR_ProductsFinder()
        counter.set_subwork(len(self._blast_results))
        for i, record in enumerate(self._blast_results):
            #setup counter
            counter[i].set_subwork(2, (10,1+5*self._include_side_annealings))
            #find_products
            mixtures = self.parallelize_work(0.1, self._find_products_in_alignment, 
                                             record.alignments, P_Finder, 
                                             counter=counter[i][0])
            if mixtures is None: return False
            #collect results and run the simulation
            PCR_Sim = self._new_PCR_Simulation()
            for hit_title, mixture_path in mixtures:
                if mixture_path is not None:
                    PCR_Sim.add_mixture(hit_title, mixture_path)
            PCR_Sim.run(counter[i][1])
            #add successful simulation to the dict 
            if PCR_Sim: self._PCR_Simulations[record.query] = PCR_Sim
        if not self._PCR_Simulations:
            print '\nNo possible PCR products have been found using BLAST hits.'
            return False
        self._have_results = True
        return True
    #end def
    
    
    def blast_and_analyze(self, counter, entrez_query=''):
        counter.set_subwork(2)
        return self.blast_query(counter[0], entrez_query) \
        and    self.simulate_PCR(counter[1])
    #end def
    
    def load_and_analyze(self, counter):
        counter.set_subwork(2)
        return self.load_results(counter[0]) \
        and    self.simulate_PCR(counter[1])
    #end def
    
    
    def _format_header(self): 
        header = wrap_text('For each hit all alignments are sorted '
                           'into "forward" and "reverse" groups. '
                           'Pairs of forward and reverse ' 
                           'alignments comprise possible PCR products.\n'
                           'If --no-exonuclease option was '
                           'provided, alignments with mismatches on ' 
                           "3'-end are ignored.\n"
                           'Quantities of products are estimated '
                           'using equilibrium equations and ' 
                           'current PCR parameters.\n')
        header += self._PCR_Simulations.values()[0].format_report_header()
        header += self._PCR_Simulations.values()[0].format_quantity_explanation()
        return header
    #end def
    
    
    def _format_report_body(self):
        body = ''
        for record_name in self._PCR_Simulations:
            body += hr(' query ID: %s ' % record_name, symbol='#')
            if len(self._PCR_Simulations[record_name].hits()) == 1:
                #all products histogram
                hit = self._PCR_Simulations[record_name].hits()[0]
                body += hr(' histogram of all possible PCR products ', symbol='=')
                body += self._PCR_Simulations[record_name].per_hit_header(hit)
                body += self._PCR_Simulations[record_name].per_hit_histogram(hit)
                body += '\n'
                body += hr(' electrophorogram of all possible PCR products ', symbol='=')
                body += self._PCR_Simulations[record_name].per_hit_electrophoresis(hit)
            else:
                #all products histogram
                body += hr(' histogram of all possible PCR products ', symbol='=')
                body += self._PCR_Simulations[record_name].all_products_histogram()
                body += '\n\n\n'
                #per hit histogram and phoresis
                body += hr(' histograms and electrophorograms of PCR products of each hit ', symbol='=')
                body += self._PCR_Simulations[record_name].all_graphs_grouped_by_hit()
        return body
    #end def
    
    
    def write_products_report(self):
        if not self._have_results: return
        #open report file
        ipcr_products = self._open_report('BLAST PCR products', self._PCR_products_filename)
        ipcr_products.write(time_hr())
        for record_name in self._PCR_Simulations:
            ipcr_products.write(hr(' query ID: %s ' % record_name, symbol='#'))
            ipcr_products.write(self._PCR_Simulations[record_name].format_products_report())
        ipcr_products.close()
        print '\nThe list of BLAST PCR products was written to:\n   ',self._PCR_products_filename
        self._add_report('BLAST PCR products', self._PCR_products_filename)
    #end def
    
        
    def write_hits_report(self):
        if not self._have_blast_results: return
        blast_report = self._open_report('blast hits report', self._hits_report_filename)
        if not blast_report: return
        #header
        blast_report.write(time_hr())
        blast_report.write(wrap_text('All hits are filtered by dG of the annealing '
                                     'of alignments and, if --no-exonuclease '
                                     'option was provided, hits with mismatches on '
                                     "3'-end are ignored.\n"))
        blast_report.write('\n')
        #filter parameters
        blast_report.write(hr(' filtration parameters '))
        if self._with_exonuclease:
            blast_report.write("DNA polymerase HAS 3'-5'-exonuclease activity\n")
        else: blast_report.write("DNA polymerase doesn't have 3'-5'-exonuclease activity\n")
        blast_report.write('Maximum dG of an alignment: %.2f kcal/mol\n' % max_dimer_dG)
        blast_report.write('\n')
        blast_report.write(hr(''))
        blast_report.write('\n\n')
        #print records
        for r, blast_record in enumerate(self._blast_results):
            num_hits = len(blast_record.alignments)
            blast_report.write(hr(' query ID: %s '  % blast_record.query, symbol='#'))
            #filter hits by alignments and format report text for each hit
            hits = []
            for h in xrange(num_hits):
                hit  = blast_record.alignments[h]
                desc = blast_record.descriptions[h] 
                #check and format hsps
                hsps = []
                for hsp in hit.hsps:
                    hsp_duplex = self._duplex_from_hsp(hsp)
                    #check if any stable dimers are formed by this duplex
                    if not hsp_duplex: continue
                    #check 3' mismatch
                    if not self._with_exonuclease \
                    and not hsp_duplex.have_3_matches: continue
                    #get primer concentration
                    for primer in self._primers:
                        if primer.has_subsequence(hsp_duplex.fwd_seq): 
                            hsp_primer_concentration = primer.concentration
                            break
                    #format hsps representation
                    hsp_text  = ('score: %d; bits: %d; E-value: %.2e;\n\n'
                                 % (hsp.score, hsp.bits, hsp.expect))
                    hsp_text += hsp_duplex.print_most_stable(include_fwd_3_mismatch=self._with_exonuclease)
                    hsp_text += 'Conversion degree = %.2f%%\n\n' \
                    % (tdf.primer_DNA_conversion_degree(hsp_primer_concentration, hsp_duplex.K)*100)
                    hsp_text += 'Template strand: '
                    if hsp.frame[1] == 1:    
                        hsp_text += 'antisense\n'
                        hsp_text += 'Position on template: %d ==> %d\n' \
                              % (hsp.sbjct_start, hsp.sbjct_end)
                    elif hsp.frame[1] == -1: 
                        hsp_text += 'sense\n'
                        hsp_text += 'Position on template: %d <== %d\n' \
                              % (hsp.sbjct_start, hsp.sbjct_end)
                    hsp_text += hr('', '.')
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
    
    def write_reports(self):
        self.write_hits_report()
        iPCR_Interface.write_reports(self)
    #end def
#end class