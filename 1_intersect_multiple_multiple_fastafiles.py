#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script takes multiple multiple fasta files downloaded from
# NCBI (a search for one specific gene for all othologs)
# and identifies the genes that are common to all files.
# It generates a common file with all fasta files and
# individual multiple fasta files for each input file.
# For the VEGFs and the purpose of this analysis, we
# have downloaded from NCBI https://www.ncbi.nlm.nih.gov/gene/XXXX/ortholog,
# XXXX = 7422, 7423, 7424, 2277, 5228 for VEGF-A, -B, -C, -D and PlGF, respectively

import os, re, sys
from Bio import SeqIO

def extract_organism(seq_record, VERBOSE = True):
    # get organism from square brackets (if there are such)
    re_search_result = re.search('\[(.*)\]', str(seq_record))
    if re_search_result:
        #if VERBOSE: print('Protein sequence detected.')
        organism = re_search_result.group(0)[1:-1]
        # get rid of three-word organisms
        organism = ' '.join(organism.split(' ')[:2])
        if VERBOSE: print('Organism (from protein record {0}): {1}.'.format(seq_record.id, organism))
    # Check whether it is a protein sequence or mRNA sequence
    elif seq_record.id[:3] in ['sp|', 'tr|']:
        #if VERBOSE: print('Protein sequence detected.')
        organism = re.search('OS=(.*) OX=', str(seq_record)).group(1)
        if VERBOSE: print('Organism (from protein record {0}): {1}.'.format(seq_record.id, organism))
    # The rest should be mRNA sequences
    else:
        #if VERBOSE: print('mRNA sequence detected.')
        if 'PREDICTED' in str(seq_record):
            new_str = str(seq_record).split('PREDICTED: ')[1]
        else:
            new_str = seq_record.description.split(' ', 1)[1]
            #print('new_str: {0}'.format(new_str))
        organism = ' '.join(new_str.split(' ', 2)[:2])
        if VERBOSE: print('Organism (from mRNA record {0}): {1}.'.format(seq_record.id, organism))
    return organism

def run():
    # This list contains species where the reference mRNA sequence and the
    # protein sequence do not agree or where some abberant splicing has been predicted, etc.
    #
    # VEGF-C issues with:
    #   Monodelphis domestica (?)
    #   Propithecus coquereli (?)
    # VEGF-A issues
    #   Jaculus jaculus (only short form)
    #   Bos indicus (low quality protein)
    #   Condylura cristata (low quality, only short form)
    #   Trichechus manatus latirostris (only short form)
    #   Neomonachus schauinslandi (only short form)
    #   Mesocricetus auratus (only short form)
    #   Elephantulus edwardii (only short form)
    # Only test the first two words of the scientific name!!!
    blacklist = ['Monodelphis domestica',
        'Propithecus coquereli',
        'Jaculus jaculus',
        'Bos indicus',
        'Condylura cristata',
        'Echinops telfairi',
        'Ceratotherium simum',
        'Fukomys damarensis',
        'Myotis davidii',
        'Myotis brandtii',
        'Physeter catodon',
        'Phascolarctos cinereus', # VEGF-A blacklist starts after this one
        'Trichechus manatus',
        'Piliocolobus tephrosceles',
        'Puma concolor',
        'Mesocricetus auratus',
        'Pteropus alecto',
        'Ochotona princeps',
        'Elephantulus edwardii',
        'Equus caballus',
        'Neophocaena asiaeorientalis',
        'Mandrillus leucophaeus',
        'Sus scrofa',
        'Equus asinus',
        'Orycteropus afer',
        'Phocoena sinus',
        'Meriones unguiculatus',
        'Cavia porcellus',
        'Ursus maritimus',
        'Marmota flaviventris',
        'Panthera tigris',
        'Pteropus vampyrus',
        'Rousettus aegyptiacus',
        'Theropithecus gelada',
        'Heterocephalus glaber',
        'Octodon degus',
        'Hipposideros armiger',
        'Chinchilla lanigera',
        'Neomonachus schauinslandi',
        'Colobus angolensis',
        'Sorex araneus',
        'Carlito syrichta',
        'Macaca mulatta',
        'Ailuropoda melanoleuca', # VEGF-D blacklist starts after this line
        'Marmota marmota',
        'Manis javanica',
        'Hylobates moloch',
        'Chelonia mydas',
        'Rattus rattus', # VEGF-B blacklist starts after this line
        'Rhinopithecus bieti']
    # Determine directories of script (in order to load & save the data files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    FASTA_COMMON_OUTFILE = ''
    #print('Input files are:')
    # remove result files from list of files
    inputfiles = ['PGF_refseq_protein.fasta',
                'PGF_refseq_transcript.fasta',
                'VEGFA_refseq_protein.fasta',
                'VEGFA_refseq_transcript.fasta',
                'VEGFB_refseq_protein.fasta',
                'VEGFB_refseq_transcript.fasta',
                'VEGFC_refseq_protein.fasta',
                'VEGFC_refseq_transcript.fasta',
                'VEGFD_refseq_protein.fasta',
                'VEGFD_refseq_transcript.fasta']
    # Generate new filename for output
    for file in inputfiles:
        FASTA_COMMON_OUTFILE += os.path.splitext(file)[0]+'_'
    FASTA_COMMON_OUTFILE = FASTA_COMMON_OUTFILE[:-1]+'_intersection.fasta'
    # Make a list of lists that has as many elements as we have input files
    dictionary_of_lists = {}
    for file in inputfiles:
        dictionary_of_lists[file] = []
        FASTA_INFILE = file
        records = list(SeqIO.parse(FASTA_INFILE, "fasta"))
        print('{0}: {1} proteins'.format(FASTA_INFILE, len(records)))
        for seq_record in records:
            # The following line contains the actual search for the organism
            # FOR PROTEIN FILES FROM UNIPROT.ORG:
            #organism = re.search('OS=(.*) OX=', str(seq_record)).group(1)
            # FOR mRNA FILES FROM NCBI:
            organism = extract_organism(seq_record, VERBOSE = False)
            if organism not in blacklist:
                if organism not in dictionary_of_lists[file]:
                    dictionary_of_lists[file].append(organism)
    # Convert into sets for intersection operation
    new_set = set(dictionary_of_lists[inputfiles[0]])
    for i in range(len(inputfiles)-1):
        new_set = new_set.intersection(set(dictionary_of_lists[inputfiles[i+1]]))
    print('Common organisms among all files ({0}):'.format(len(new_set)))
    for item in new_set:
        print(item)
    list_of_organisms = list(new_set)
    # Write all entries to a common fasta file
    new_common_records = []
    for file in inputfiles:
        FASTA_INFILE = file
        FASTA_OUTFILE = os.path.splitext(FASTA_INFILE)[0]+'_intersection.fasta'
        records = list(SeqIO.parse(FASTA_INFILE, "fasta"))
        new_records = []
        # This is more complicated because we want to have all organisms
        # in the same order in both protein and mRNA file
        for organism in list_of_organisms:
            for seq_record in records:
                extracted_organism = extract_organism(seq_record, VERBOSE = False)
                if extracted_organism == organism:
                    new_common_records.append(seq_record)
                    new_records.append(seq_record)
        SeqIO.write(new_records, FASTA_OUTFILE, "fasta")
        print('Output written to {0}'.format(FASTA_OUTFILE))
    # In case we want to have all sequences in one file
    #SeqIO.write(new_common_records, FASTA_COMMON_OUTFILE, "fasta")
    #print('Output written to {0}'.format(FASTA_COMMON_OUTFILE))
    print('Common organisms among all files: {0}.'.format(len(new_set)))

if __name__ == '__main__':
    run()
