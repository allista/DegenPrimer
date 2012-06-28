'''
Created on Jun 26, 2012

@author: allis
'''

from Bio.Blast import NCBIWWW, NCBIXML


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


def write_blast_report(blast_results, blast_report_filename):
    blast_report = open(blast_report_filename, 'w')
    for blast_record in blast_results:
        blast_report.write('Query:        %s\n' % blast_record.query)
        blast_report.write('Query length: %s\n' % blast_record.query_length)
        blast_report.write('Hits:         %s\n' % len(blast_record.alignments))
        blast_report.write('-'*80+'\n')
        for a in range(len(blast_record.alignments)):
            alignment = blast_record.alignments[a]
            blast_report.write(('\n'+'*'*15+' Hit #%d '+'*'*15+'\n') % (a+1))
            blast_report.write(alignment.title+'\n')
            blast_report.write('Length:     %d\n' % alignment.length)
            blast_report.write('Max bits:   %d\n' % blast_record.descriptions[a].bits)
            blast_report.write('Alignments: %d\n' % len(alignment.hsps))
            for hsp in alignment.hsps:
                if hsp.align_length < blast_record.query_length/2 \
                and blast_record.query_length-hsp.query_end > 3:
                    blast_report.write('\nAlignment is shorter than half of the query and ends far from 3\'.\n')
                    blast_report.write('Length: %(len)d; End: %(al_end)d; Bits: %(bits)d; Location: [%(start)d, %(end)d]\n' % \
                                       {'len'   : hsp.align_length, 
                                        'start' : hsp.sbjct_start, 
                                        'end'   : hsp.sbjct_end, 
                                        'bits'  : hsp.bits,
                                        'al_end': hsp.query_end})
                    continue
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
        blast_report.write('-'*80+'\n\n')
    blast_report.close()
    print 'Blast report was written to: '+blast_report_filename