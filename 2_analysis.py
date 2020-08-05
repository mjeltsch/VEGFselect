#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script extracts the n most informative sequences from a multiple fasta file
# of protein sequences (e.g. VEGFA_cleaned_intersection.fasta) and generates an
# alignemnt.
#
# The alignment is trimmed manually (from the first to the last cysteine of the PDGF/VEGF cysteine signature +15 aa N-terminally
# and 5 aa C-terminally. The the alignments are merged (to identify the common VEGF homology domain) and in part 2,
# the analysis for purifying and adaptive evolution is prepared.
# The analysis is done by another software (HyPhy, Mac only) using a maximum-likelyhood (ML) approach to infer nonsynonymous
# and synonymous substitution rates on a per-site basis for the coding alignment.
# Then the resulting data are visualized with the seaborn package by the script 4_grafic.py

import os, re, sys, glob, shutil
from Bio import SeqIO
from phylolib import execute_subprocess
from shutil import copyfile

def run():
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    n = 50
    MODE = 'mcoffee'
    input_genes = ['PGF', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD']
    # Make directory where all final alignments are gathered
    if not os.path.isdir('VEGF_all'):
            os.mkdir('VEGF_all')
    final_string_of_alignments = ''
    for input_gene in input_genes:
        input_file = APPLICATION_PATH+'/'+input_gene+'_refseq_protein_intersection.fasta'
        # Make a directory to hold all data for a specific protein
        outputdirectory = APPLICATION_PATH+'/'+input_gene
        if not os.path.isdir(outputdirectory):
            os.mkdir(outputdirectory)
        # The sequences in th ekeep file are ALWAYS included in the alignment (default: mouse, human, rat)
        keepfile = input_gene+'_keep.fasta'
        # Generate filenames
        outputfile = outputdirectory+'/'+input_gene+'_informative_subset.fasta'
        encoded_file = outputdirectory+'/'+input_gene+'_encoded.fasta'
        decoded_file = outputdirectory+'/'+input_gene+'_decoded.fasta'
        alignment_file = outputfile+'_aln'
        mRNA_file = APPLICATION_PATH+'/'+input_gene+'_refseq_transcript_intersection.fasta'
        mRNA_file_informative = outputdirectory+'/'+input_gene+'_refseq_transcript_informative.fasta'
        tree_file = outputdirectory+'/'+input_gene+'_final_mRNA_alignment.phylip_phyml_tree.txt'
        codename_file = outputdirectory+'/code_names.list'
        trimmed_file = outputdirectory+'/'+input_gene+'_trimmed.fasta'
        print('Input file is: {0}'.format(input_file))

        # TRIMMING PART 1: MOST INFORMATIVE PROTEIN SEQUENCES (This might not be ok! Test with and without this option!)
        if os.path.isfile(keepfile):
            bash_command = 'nice -n 19 t_coffee -other_pg seq_reformat -in {0} -in2 {1} -action +trim _seq_n{2} -output fasta_seq > {3}'.format(input_file, keepfile, n, outputfile)
        else:
            bash_command = 'nice -n 19 t_coffee -other_pg seq_reformat -in {0} -action +trim _seq_n{2} -output fasta_seq > {3}'.format(input_file, keepfile, n, outputfile)
        comment = 'Extracting the {0} most informative sequences from {1}'.format(n, input_file)
        result, error = execute_subprocess(comment, bash_command)
        # TRIMMING PART 2: REMOVING CORRESPONDING mRNA SEQUENCES
        protein_records = list(SeqIO.parse(outputfile, "fasta"))
        transcript_records = list(SeqIO.parse(mRNA_file, "fasta"))
        new_transcript_records = []
        for protein_record in protein_records:
            organism = re.search('\[(.*)\]', str(protein_record)).group(0)[1:-1]
            for transcript_record in transcript_records:
                if organism in transcript_record.description:
                    new_transcript_records.append(transcript_record)
                    break
        SeqIO.write(new_transcript_records, mRNA_file_informative, "fasta")

        # ALIGNMENT
        bash_command = 'nice -n 19 t_coffee -seq {0} -output=html,fasta -mode {1}'.format(outputfile, MODE)
        comment = 'Making MSA for file {0}'.format(outputfile)
        result, error = execute_subprocess(comment, bash_command, working_directory = outputdirectory)
        
        # # MAKE CODE LIST
        # bash_command = 't_coffee -other_pg seq_reformat -in {0} -output code_name > {1}'.format(alignment_file, codename_file)
        # comment = 'Converting fasta descriptions part 1 (creating code list) with t_coffee.'
        # result, error = execute_subprocess(comment, bash_command, working_directory = outputdirectory)
        # # ENCODE
        # bash_command = 't_coffee -other_pg seq_reformat -code {0} -in {1} > {2}'.format(codename_file, alignment_file, encoded_file)
        # comment = 'Converting fasta descriptions part 2 (replacing fasta descriptions with codes) with t_coffee.'
        # result, error = execute_subprocess(comment, bash_command, working_directory = outputdirectory)
        # # REMOVE EVERYTHING FROM FASTA DESCRIPTION EXCEPT MACHINE-READABLE CODE
        # records = list(SeqIO.parse(encoded_file, "fasta"))
        # new_records = []
        # for seq_record in records:
        #     print('New sequence id (encoded): {0}'.format(seq_record.id))
        #     print('Description: {0}'.format(seq_record.description))
        #     seq_record.description = ''
        #     new_records.append(seq_record)
        # SeqIO.write(new_records, encoded_file, "fasta")

        # TRIMMING
        trim_dictionary = {'PGF': ['NP_002623.2', 36, 135], 'VEGFA': ['NP_003367.4', 216, 315], 'VEGFB': ['NP_001230662.1', 31, 129], 'VEGFC': ['NP_005420.1', 114, 216], 'VEGFD': ['NP_004460.1', 94, 196]}
        reference_seq = trim_dictionary[input_gene][0]
        start = trim_dictionary[input_gene][1]
        end = trim_dictionary[input_gene][2]
        bash_command = 't_coffee -other_pg seq_reformat -in {0} -action +extract_block \'{1}\' {2} {3} > {4}'.format(alignment_file, reference_seq, start, end, trimmed_file)
        comment = 'Trimming sequences according to trim_dictionary'
        result, error = execute_subprocess(comment, bash_command, working_directory = outputdirectory)
        final_string_of_alignments += trimmed_file+' '
        copy_target = 'VEGF_all/{0}.fasta'.format(input_gene)
        copyfile(trimmed_file, copy_target)

    # MERGING ALIGNMENTS WITH MAFFT (https://mafft.cbrc.jp/alignment/software/merge.html)
    bash_command = 'cat {0}> VEGF_all/VHD.fasta'.format(final_string_of_alignments)
    comment = 'Concatenating MSAs...'
    result, error = execute_subprocess(comment, bash_command)
    bash_command = 'ruby makemergetable.rb {0}> subMSAtable.txt'.format(final_string_of_alignments)
    comment = 'Preparing for alignment merging...'
    result, error = execute_subprocess(comment, bash_command)
    bash_command = 'mafft --merge subMSAtable.txt VEGF_all/VHD.fasta > VEGF_all/VHD_aligned.fasta'.format(final_string_of_alignments)
    comment = 'Merging individual alignments with MAFFT...'
    result, error = execute_subprocess(comment, bash_command)

    #
    # PART 2 
    #    

    n = 50
    MODE = 'mcoffee'
    input_genes = ['PGF', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD']
    # Before running this script, you need to manually specify (trim) the part of the alignment that you want 
    # to use. The full (original) alignment is in the file with the naming format (example)
    # "VEGFB_refseq_protein_intersection_encoded.fasta"
    # The trimmed alignment shoul dbe in the format
    # "VEGFB_refseq_protein_intersection_encoded_trimmed.fasta"
    # Trimming can be done e.g. with t_coffee:
    # t_coffee -other_pg seq_reformat -in VEGFB_refseq_protein_intersection_encoded.fasta  -action +extract_block cons 395 499 > VEGFB_refseq_protein_intersection_encoded_trimmed.fasta
    # In this script, the trimming is done in part 1 based on the trim_dictionary. 
    #
    for input_gene in input_genes:
        codename_file = APPLICATION_PATH+'/'+input_gene+'/'+'code_names.list'
        gblocks_output = APPLICATION_PATH+'/'+input_gene+'/'+input_gene+'_trimmed.fasta'
        decoded_file = APPLICATION_PATH+'/'+input_gene+'/'+input_gene+'_intersection_decoded.fasta'
        outputdirectory = APPLICATION_PATH+'/'+input_gene+'/'
        mRNA_file_informative = outputdirectory+input_gene+'_refseq_transcript_informative.fasta'
        final_fasta = outputdirectory+input_gene+'_final_mRNA_alignment.fasta'
        final_fasta_phylip = outputdirectory+input_gene+'_final_mRNA_alignment.phylip'
        tree_file = outputdirectory+input_gene+'_final_mRNA_alignment.phylip_phyml_tree.txt'
        paml_output_file_baseml = outputdirectory+input_gene+'_final_baseml.paml'
        paml_output_file_codeml = outputdirectory+input_gene+'_final_codeml.paml'
        paml_control_file_baseml = outputdirectory+input_gene+'_paml_baseml'
        paml_control_file_codeml = outputdirectory+input_gene+'_paml_codeml'
        codename_file = outputdirectory+'code_names.list'
        yn00_control_file = outputdirectory+'yn00.ctl'
        yn00_infile = outputdirectory+'yn'
        yn00_outfile = outputdirectory+input_gene+'.yn'
        print('Input gene is: {0}'.format(input_gene))

        ## converting gblocks output into human-readable ids 
        #bash_command = 't_coffee -other_pg seq_reformat -decode {0} -in {1} > {2}'.format(codename_file, gblocks_output, decoded_file)
        #comment = 'Decoding tree file file into human-readable format using the following command:\n'
        #result, error = execute_subprocess(comment, bash_command)
        # pal2nal
        bash_command = 'pal2nal.pl -output fasta {0} {1} > {2}'.format(gblocks_output, mRNA_file_informative, final_fasta)
        comment = 'Generating mRNA alignment corresponding to protein alignment using pal2nal:\n'
        result, error = execute_subprocess(comment, bash_command)

        # TREEBUILDING
        comment = 'Convert into phylip using the following command:\n'
        bash_command = 't_coffee -other_pg seq_reformat -in {0} -output phylip_aln > {1}'.format(final_fasta, final_fasta_phylip)
        result, error = execute_subprocess(comment, bash_command)

        # Detect whether parallel bootstrapping should be performed
        mpirun_path = shutil.which('mpirun')
        phymlmpi_path = shutil.which('phyml-mpi')
        if mpirun_path != '' and phymlmpi_path != '':
            # b = 100: with 100 bootstrap replicates
            # The core count should be done dynamically (e.g. by interrogating with lscpu)
            bash_command = 'mpirun -n 2 phyml-mpi -i {0} -d nt -b 100'.format(final_fasta_phylip)
        else:
            # b = -1: approximate likelihood ratio test returning aLRT statistics
            bash_command = 'phyml -i {0} -d nt -b -1'.format(final_fasta_phylip)
            #
            # other options:
            # b = -5: approximate Bayes branch supports

        comment = 'Make tree with the following command:\n'
        result, error = execute_subprocess(comment, bash_command)

        # PAML (use baseml and codeml with seqtype =1, which specifies nucleotide)
        #
        # Create control file for PAML baseml
        control_text ='''       seqfile = {0}
      treefile = {1}
      outfile = {2}   *in result file

        noisy = 9   * 0,1,2,3: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI

        model = 4   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu

        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.3   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 8   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK

        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:TipDate
        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*        ndata = 5
*        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
*    readfpatt = 0  * read site pattern frequences instead of sequences
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
        method = 0  * 0: simultaneous; 1: one branch at a time'''.format(final_fasta, tree_file, paml_output_file_baseml)
        with open(paml_control_file_baseml, 'w+') as file:
            file.write(control_text)
        comment = 'Running PAML baseml:\n'
        bash_command = 'baseml {0}'.format(paml_control_file_baseml)
        result, error = execute_subprocess(comment, bash_command)
        #
        # Create control file for PAML codeml
        control_text ='''       seqfile = {0}
      outfile = {1}
      treefile = {2}

        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1  * 1: detailed output, 0: concise output
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 0
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0    * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0  * 0:rates, 1:separate;

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0 * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

        clock = 3   * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
*    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*        ndata = 10
        method = 1   * 0: simultaneous; 1: one branch at a time'''.format(final_fasta, paml_output_file_codeml, tree_file)
        with open(paml_control_file_codeml, 'w+') as file:
            file.write(control_text)
        comment = 'Running PAML codeml:\n'
        bash_command = 'codeml {0}'.format(paml_control_file_codeml)
        #result, error = execute_subprocess(comment, bash_command)

        # Create control file for yn00
        control_text = '''seqfile = {0}   * sequence data file name
    outfile = {1}   * main result file
    verbose = 1   * 1: detailed output (list sequences), 0: concise output
    icode = 0   * 0: universal code; 1: mammalian mt; 2-10: see below
    weighting = 0   * weighting pathways between codons (0/1)?
    commonf3x4 = 0   * use one set of codon freqs for all pairs (0/1)?
    * ndata = 1

    * Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
    * 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
    * 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
    * 10: blepharisma nu.
    * These codes correspond to transl_table 1 to 11 of GENEBANK.'''.format(final_fasta, yn00_outfile)
        with open(yn00_control_file, 'w+') as file:
            file.write(control_text)
        # PAML (yn00)
        bash_command = 'yn00 {0}'.format(yn00_control_file)
        comment = 'Using Yang & Nielsen (2000) method to estimate synonymous and nonsynonymous substitution rates:\n'
        result, error = execute_subprocess(comment, bash_command, working_directory = outputdirectory)
        #os.rename(yn00_infile, yn00_outfile)

if __name__ == '__main__':
    run()
