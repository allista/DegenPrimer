'''
Created on Jun 26, 2012

@author: allis
'''

from time import ctime
from Bio.Blast import NCBIWWW, NCBIXML
from SecStructures import hr

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


def write_blast_report(blast_results, blast_report_filename, top_hits=10, top_hsps=4):
    blast_report = open(blast_report_filename, 'w')
    blast_report.write(hr(' %s ' % ctime(), symbol='#'))
    for blast_record in blast_results:
        blast_report.write('\n')
        blast_report.write(hr(' %s ' % blast_record.query, symbol='#'))
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