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
Created on Jul 25, 2012

@author: Allis Tauri <allista@gmail.com>
'''

#imports
from DegenPrimer.OligoFunctions import generate_unambiguous, self_complement
from DegenPrimer import TD_Functions    
from DegenPrimer.TD_Functions import calculate_Tm, source_feature, add_PCR_conditions, format_PCR_conditions
from DegenPrimer.Blast import Blast
from DegenPrimer.SecStructures import SecStructures
from DegenPrimer.StringTools import hr, wrap_text, print_exception, time_hr
from DegenPrimer.iPCR import iPCR
try:
    from Bio import SeqIO
except Exception, e:
    print_exception(e)
    raise ImportError('The BioPython must be installed in your system.')
############################################################################### 


def degen_primer_pipeline(args):
    #set job ID and primers list
    primers = args.primers
    job_id  = args.job_id
    
    #set concentrations
    TD_Functions.C_Na   = args.Na
    TD_Functions.C_Mg   = args.Mg
    TD_Functions.C_dNTP = args.dNTP
    TD_Functions.C_DNA  = args.DNA
    TD_Functions.C_Prim = args.Primer
    TD_Functions.C_DMSO = args.DMSO
    TD_Functions.PCR_T  = args.PCR_T
    
    #check if at least one primer is provided
    if not args.sense_primer and not args.antisense_primer:
        print 'At least one primer (sense or antisense) should be provided.'
        return 1
    #test for self-complementarity
    for primer in primers:
        if primer and self_complement(primer[0]):
            print 'Error: %s primer "%s" [%s] is self-complementary.' % \
                    (primer[0].description, primer[0].id, primer[0].seq)
            print 'You should not use self-complementary oligonucleotides as primers.\n'
            return 1
    #save the configuration only after preliminary checks
    args.save_configuration()
    print ''
    #-----------------------------------------------------------------------------#
    
    
    #generate unambiguous primer sets
    for primer in primers:
        try:
            if primer: primer[1] = generate_unambiguous(primer[0])
        except Exception, e:
            print_exception(e)
            return 1
    #-----------------------------------------------------------------------------#
    
    
    #primer lists
    def pfam_primers(pfam):
        if not primers[pfam]: return [] 
        if not primers[pfam][1]:
            return [primers[pfam][0]]
        else: return primers[pfam][1]
    #end def
    
    def all_primers():
        return pfam_primers(0) + pfam_primers(1)
    #-----------------------------------------------------------------------------#
    
    
    #calculate melting temperatures for primers
    for primer in primers:
        if not primer: continue
        #if original primer is not a degenerate
        if not primer[1]: 
            primer_Tm = calculate_Tm(primer[0])
            primer[2]['Tm'] = primer_Tm
            continue
        #else...
        n_primers = len(primer[1])
        mean_Tm, min_Tm, max_Tm = 0, 10000, 0 
        for p in range(n_primers):
            primer_Tm = calculate_Tm(primer[1][p])
            if primer_Tm:
                primer_Tm = primer_Tm
                mean_Tm  += primer_Tm
                if min_Tm >= primer_Tm: min_Tm = primer_Tm
                if max_Tm <  primer_Tm: max_Tm = primer_Tm
            else: n_primers -= 1
        feature = source_feature(primer[0])
        add_PCR_conditions(feature)
        feature.qualifiers['Tm_min']  = str(min_Tm)
        feature.qualifiers['Tm_max']  = str(max_Tm)
        feature.qualifiers['Tm_mean'] = str(mean_Tm/n_primers)
        primer[2]['Tm_min']  = min_Tm
        primer[2]['Tm_max']  = max_Tm
        primer[2]['Tm_mean'] = mean_Tm/n_primers
    #-----------------------------------------------------------------------------#
    
    
    #fasta file with primers
    fmt = 'fasta'
    fmt_filename = job_id+'.'+fmt
    primers_list = list()
    for primer in primers:
        if not primer: continue
        primers_list.append(primer[0])
        primers_list += primer[1]
    SeqIO.write(primers_list, fmt_filename, fmt)
    print '\nUnambiguous primers were written to:\n   ',fmt_filename
    #-----------------------------------------------------------------------------#
    
    
    #check for hairpins, dimers and cross-dimers
    #write full and short reports
    structures_full_report_filename  = job_id+'-full-report.txt'
    structures_short_report_filename = job_id+'-short-report.txt'
    full_structures_file  = open(structures_full_report_filename, 'w')
    short_structures_file = open(structures_short_report_filename, 'w')
    #write header
    for f in (full_structures_file, short_structures_file):
        f.write(time_hr())
        f.write(wrap_text('For each degenerate primer provided, a set of unambiguous primers is generated. '
                          'For each such set the minimum, maximum and mean melting temperatures are calculated. '
                          'For each primer in each set stable self-dimers and hairpins are predicted. '
                          'For every possible combination of two unambiguous primers cross-dimers are also predicted. '
                          'If an unambiguous primer is provided, it is treated as a set with a single element.\n\n'))
        f.write(hr(' PCR conditions '))
        f.write(format_PCR_conditions()+'\n')
        f.write(hr(' primers and their melting temperatures '))
        #temperatures
        temps  = []
        
        if primers[0] and primers[1]:
            max_id = max(len(primers[0][0].id), len(primers[1][0].id))
        else: max_id = 0
        for primer in primers:
            if not primer: continue
            spacer = max_id-len(primer[0].id)
            if not primer[1]: #not degenerate
                f.write('%s:%s  %02db  %s\n' % (primer[0].id, ' '*spacer, len(primer[0].seq), str(primer[0].seq)) +\
                        '   Tm:      %.1f C\n' % primer[2]['Tm'])
                temps.append(primer[2]['Tm'])
            else:
                f.write('%s:%s  %02db  %s\n' % (primer[0].id, ' '*spacer, len(primer[0].seq), str(primer[0].seq)) +\
                        '   Tm max:  %.1f C\n' % primer[2]['Tm_max'] +\
                        '   Tm mean: %.1f C\n' % primer[2]['Tm_mean'] +\
                        '   Tm min:  %.1f C\n' % primer[2]['Tm_min'])
                temps.append(primer[2]['Tm_min'])
        #warning
        if len(temps) == 2:
            if abs(temps[0]-temps[1]) >= 5:
                f.write('\nWarning: lowest melting temperatures of sense and antisense primes \n'
                        '         differ more then by 5C\n')
        f.write('\n\n')
        f.write(hr(' stable secondary structures ', symbol='#'))
    #short report thresholds
    short_structures_file.write(hr(' only following structures are reported '))
    short_structures_file.write("   dimers with free energy below %.2f kcal/mol are reported\n" \
                                % (args.sec_dG   if args.sec_dG<= 0 else 0))
    short_structures_file.write("3'-dimers with free energy below %.2f kcal/mol are reported\n" \
                                % (args.sec_dG+1 if args.sec_dG<=-1 else 0))
    short_structures_file.write("   hairpins with free energy below %.2f kcal/mol are reported\n" \
                                % (args.sec_dG+2 if args.sec_dG<=-2 else 0))
    short_structures_file.write("3'-hairpins with free energy below %.2f kcal/mol are reported\n" \
                                % (args.sec_dG+3 if args.sec_dG<=-3 else 0))
    short_structures_file.write('\n')
    #check all self-dimers and hairpins
    all_primers = all_primers()
    for primer in all_primers:
        structs = SecStructures(primer, dG_threshold=args.sec_dG)
        full_structures_file.write(str(structs))
        if structs.belowThreshold():
            short_structures_file.write(structs.formatShort())
    #check ALL cross-dimers
    num_primers = len(all_primers)
    if num_primers > 1:
        for f in (full_structures_file, short_structures_file):
            f.write('\n'+hr(' stable cross dimers ', symbol='#'))
        for i in range(num_primers):
            for j in range(i+1,num_primers):
                structs = SecStructures(all_primers[i],all_primers[j], 
                                        dG_threshold=args.sec_dG)
                full_structures_file.write(str(structs))
                if structs.belowThreshold():
                    short_structures_file.write(structs.formatShort())
    full_structures_file.close()
    short_structures_file.close()
    print '\nFull report with all secondary structures was written to:\n   ',structures_full_report_filename
    print '\nShort report with a summary of secondary structures was written to:\n   ',structures_short_report_filename
    args.register_report('Tm and secondary structures', structures_short_report_filename)
    #-----------------------------------------------------------------------------#

    
    #ipcress program file and test for primers specificity by iPCR
    #this is only available for pairs of primers
    if args.sense_primer and args.antisense_primer:
        fwd_primers  = pfam_primers(0)
        rev_primers  = pfam_primers(1)
        ipcr = iPCR(job_id, fwd_primers, rev_primers, args.min_amplicon, args.max_amplicon, args.no_exonuclease)
        ipcr.writeProgram()
        #if target sequences are provided and the run_icress flag is set, run iPCR...
        if args.run_ipcress and args.fasta_files:
            ipcr.executeProgram(args.fasta_files, args.max_mismatches)
        else:
            #load previously saved results
            print '\nLoading raw iPCR results from the previous ipcress run.'
            ipcr.load_results()
        if ipcr.have_results:
            ipcr.write_PCR_report()
            ipcr.register_reports(args)
    #-----------------------------------------------------------------------------#
    
    
    #test for primers specificity by BLAST
    blast = Blast(job_id)
    #if --do-blast command was provided, make an actual query
    if args.do_blast:
        #construct Entrez query
        entrez_query = ''
        if args.organisms:
            for organism in args.organisms:
                if entrez_query: entrez_query += ' OR '
                entrez_query += organism+'[organism]'
        #do the blast
        blast.blast_primers(all_primers(), entrez_query)
    #otherwise, reparse previously saved results with current parameters
    else:
        #load blast results
        blast.load_results()
    if blast.have_results:
        #report results in human readable form
        blast.write_hits_report(args.sec_dG,
                                args.no_exonuclease)
        blast.write_PCR_report(args.min_amplicon, 
                               args.max_amplicon, 
                               args.no_exonuclease)
        blast.register_reports(args)
    #-----------------------------------------------------------------------------#
    
    print '\nDone\n'
    return 0
#end def