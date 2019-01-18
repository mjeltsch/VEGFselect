#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#
#
# UBUNTU PACKAGES
#
# This script has been developed to run under Linux and was specifically tested ONLY under Ubuntu 16.04 and 18.04.
# It requires the installation of the following ubuntu packages:
#
# python3-setuptools
# clustalw
# mafft
# dialign-tx
# poa
# probcos
# muscle
# kalign
# amap-align
# proda
# prank
# t-coffee
# phyml
#
# sudo apt -y --show-progress install inkscape python3-setuptools python3-pyqt4 python3-pyqt4.qtopengl python3-pip autoconf t-coffee clustalw mafft dialign-tx poa probcons muscle kalign amap-align proda prank phyml t-coffee imagemagick build-essential libblas-dev liblapack-dev zlib1g-dev libcairo2-dev libcurl4-openssl-dev python3-numpy python3-lxml python3-six
#
#
# For some reason the Ubuntu package names some of the alignment executables differently and some links need to be
# created in order for t-coffee to find them:
#
# sudo ln -s /usr/bin/dialign-tx /bin/dialign-t
# sudo ln -s /usr/bin/clustalw /bin/clustalw2
#
#
# MANUAL INSTALLATIONS
#
# In addition, the pcma executable need to be manually downloaded from http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# and compiled because it is not available from the Ubuntu repository:
#
# wget http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# tar -xvzf pcma.tar.gz
# cd pcma
# make
# sudo cp pcma /bin
#
#
# OPTIONAL INSTALLS
#
# The multicore-enabled version of phyml (phyml-mpi) is not available as a precompiled Ubuntu package and needs
# to be installed manually, but they single-core version works as well (is just slower). The command to execute
# the multicore version is: "mpirun -n 4 phyml-mpi -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b -1"
# If this script is run on a (headless) server, the xvfb package is required since the ete3 package requires the presence of x.org.
#
#
# PYTHON MODULES
#
# The following Python modules need to be installed:
#
# biopython
# ete3
#
# sudo pip3 install biopython ete3

import argparse, subprocess, Bio, os, sys, shutil, re
#from Bio.Blast import NCBIWWW
#from Bio.Blast import NCBIXML
#from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio import Entrez
from Bio import SeqIO
from Bio import Phylo
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, SequenceFace, ImgFace, SVGFace, faces, add_face_to_node

def print_subprocess_result(name, out, err):
    if out.decode('utf-8') != '':
        print("Output from " + name + ": " + str(out))
#    if err != None:
    if err.decode('UTF-8') != '':
        print("Error from " + name + ": " + str(err))

def execute_subprocess(comment, bash_command):
    print("\n" + comment, bash_command)
    process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    process_status = process.wait()
    if output.decode('utf-8') != '':
        print("Output: " + str(output))
    if error.decode('UTF-8') != '':
        print("Error: " + str(error))

# Gather all sequences into one fasta file for t_coffee
list_of_sequences = []
filelist = os.listdir('cDNAs/')
for file in filelist:
    if file[-2:] == "gb":
        SeqIO.convert('cDNAs/'+file, "genbank", 'cDNAs/'+file[:-2]+'fa', "fasta")
    if file[-7:] == "_ref.fa":
        sequence = SeqIO.read('cDNAs/'+file, "fasta")
        sequence_length = len(sequence.seq)
        REFERENCE_SEQUENCE = sequence.id

FASTA = "sequences.fasta"
FASTA_ALIGNED = "sequences_aligned.fasta"
FASTA_ALIGNED_TRIMMED = "sequences_aligned_trimmed.fasta"
bash_command = 'cat cDNAs/*.fa > {0}'.format(FASTA)
execute_subprocess("Concatenating all fasta sequeces into one file.", bash_command)

def align():
    execute_subprocess(
        "Generating multiple sequence alignment with the following command:",
        't_coffee {0} -outfile {1} -output=fasta_aln -mode mcoffee'.format(FASTA, FASTA_ALIGNED))

    execute_subprocess(
        "Trimming multiple sequence alignment with the following command:",
        't_coffee -other_pg seq_reformat -in {0} -action +extract_block {1} 1 {2} > {3}'.format(FASTA_ALIGNED, REFERENCE_SEQUENCE, sequence_length, FASTA_ALIGNED_TRIMMED))

def treebuilding():
    execute_subprocess(
        "Converting fasta descriptions part 1 (creating code list) with t_cofeee using the following command:",
        "t_coffee -other_pg seq_reformat -in " + FASTA_ALIGNED_TRIMMED + " -output code_name > code_names.list")

    execute_subprocess(
        "Converting fasta descriptions part 2 (replacing fasta descriptions with codes) with t_cofeee using the following command:",
        "t_coffee -other_pg seq_reformat -code code_names.list -in " + FASTA_ALIGNED_TRIMMED + " > " + FASTA_ALIGNED_TRIMMED_CODED)

    execute_subprocess(
        "Convert into phylip using the following command:",
        "t_coffee -other_pg seq_reformat -in " + FASTA_ALIGNED_TRIMMED_CODED + " -output phylip_aln > " + PHYLIP_ALIGNED_TRIMMED_CODED)

    # Detect whether parallel bootstrapping should be performed
    mpirun_path = shutil.which('mpirun')
    phymlmpi_path = shutil.which('phyml-mpi')
    if mpirun_path != '' and phymlmpi_path != '':
        phylo_command = "mpirun -n 4 phyml-mpi -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b 100"
    else:
        phylo_command = "phyml -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b -1"

    # The gene tree building is actually never used since the species tree is used for the tree drawing.
    # We anyway calculate it to be able to compare gene and species trees.
    execute_subprocess(
        "Make tree with the following command:",
        phylo_command)

    # phyml adds or doesn't add the .txt extension to the output file (depending on the version) and we need to check for this!
    phyml_output_file = PHYLIP_ALIGNED_TRIMMED_CODED + "_phyml_tree"
    if os.path.isfile(phyml_output_file):
        os.rename(phyml_output_file, phyml_output_file + ".txt")
    execute_subprocess(
        "Decoding tree file file into human-readable format using the following command:",
        "t_coffee -other_pg seq_reformat -decode code_names.list -in " + phyml_output_file + ".txt > " + PHYLIP_ALIGNED_TRIMMED_DECODED)

def run():
    sequence_dictionary =  {'XP_022098898.1':["Starfish.svg", "Starfish", "Acanthaster planci", "Asteroidea"],
             'XP_007894115.1':["Callorhinchus_milii.svg", "Ghost shark", "Callorhinchus milii", "Chondrichthyes"],
    #         'NP_002599.1':["Homo_sapiens.svg", "Human PDGF-B", "Homo sapiens", "Mammalia"],
             'XP_020376152.1':["Rhincodon_typus.svg", "Whale shark", "Rhincodon typus", "Chondrichthyes"],
             'ENSEBUT00000000354.1':["Eptatretus_burgeri.svg", "Hagfish", "Eptatretus burgeri", "Myxini", '''>ENSEBUT00000000354.1 PREDICTED: VEGF-C Inshore hagfish [Eptatretus burgeri]
    LAIDVLHLHIHPDYLQDNEDIQTDHDPWEIIDTDTFPKGALGPKRIERLTRRLLAASSVD
    DLLTLLYPWPEEATAQRCRRGHRTEPQFQAAVVNINWEAIELEWSNTLCAPRQACVPTGP
    DSHSVERSLHYRPPCVSLHRCTGCCNDPRRSCTSTAVQHVSKTVIEISLFPELVIRPVTI
    SYKNHTECHCLTIPFHNVRPPRSVSKTWRDGGQCGPTSGSCAKGTSWNVEACRCVAQQGV
    GEVCGPGMKWNEEMCNCVCWRVCPRGQRLHTQSCGCECALNTRDCFLRARRFDRRKCRCV
    TAPCPGAPEVCPVGLGFSEELCRCVPQDWIQGLQRNGG\n'''],
             'LS-transcriptB2-ctg17881':['Leucoraja_erinacea.svg', 'Skate', "Leucoraja erinacea", "Chondrichthyes", '''>LS-transcriptB2-ctg17881 PREDICTED: VEGF-C Little skate [Leucoraja erinacea]
    RDQAHSQGQATSQLEQQLRSAASIIELMDIFYPEYRRIQECLQRRSTMAKHARREVEEEQ
    EEEEEEEWTEAAAFTVLWREEDLRNIELEWERTQCKPREVCLDLGRELGTATNNFYKPPC
    VSVHRCGGCCNNEGFQCINVSTAFVSKTLMEITIPQVGLSRPVVISFINHTACGCHPRHI
    FSHSHSIIRRSFHVSPTSCVMGNETCPRGHHWDPHHCGCVSVHEVAAPPASTAEPDVTEG
    EFDDFCGPYMVFDEDSCSCVCTNRPSSCHPSKEFDENTCRCVCFNRQHRGLCREEEQEEW
    DDDACQCVCRKSCPRHLPLNTNTCTCECSESPASCFRRGKKFDPYTCRCYRLPC\n'''],
             'XP_006632034.2':["Spotted_gar.svg", "Spotted gar", "Lepisosteus oculatus", "Actinopterygii"],
             'NP_001167218.1':["Salmo_salar.svg", "Atlantic salmon", "Salmo salar", "Actinopterygii"],
             'NP_991297.1':["Danio_rerio.svg", "Zebrafish", "Danio rerio", "Actinopterygii"],
             'XP_011610643.1':["Takifugu_rubripes.svg", "Fugu fish", "Takifugu rubripes", "Actinopterygii"],
             'XP_020464669.1':["Monopterus_albus.svg", "Swamp eel", "Monopterus albus", "Actinopterygii"],
             'XP_015809835.1':["Nothobranchius_furzeri.svg", "Killifish", "Nothobranchius furzeri", "Actinopterygii"],
             'XP_023189044.1':["Platyfish.svg", "Platy fish", "Xiphophorus maculatus", "Actinopterygii"],
             'XP_007564695.1':["Poecilia_formosa.svg", "Amazon molly", "Poecilia formosa", "Actinopterygii"],
             'XP_006006690.1':["Coelacant.svg", "Coelacant", "Latimeria chalumnae", "Sarcopterygii"],
             'XP_002933363.1':["Xenopus_tropicalis.svg", "Xenopus", "Xenopus tropicalis", "Amphibia"],
             'XP_018419054.1':["Nanorana_parkeri.svg", "Tibet frog", "Nanorana parkeri", "Amphibia"],
             'XP_015283812.1':["Gekko_japonicus.svg", "Gekko", "Gekko japonicus", "Reptilia"],
             'ETE60014.1':["Ophiophagus_hannah.svg", "King cobra", "Ophiophagus hannah", "Reptilia"],
             'XP_003221689.1':["Lizard.svg", "Anole lizard", "Anolis carolinensis", "Reptilia"],
             'XP_005304228.1':["Turtle.svg", "Painted turtle", "Chrysemys picta bellii", "Reptilia"],
             'XP_006276984.1':["Alligator_mississippiensis.svg", "American Alligator", "Alligator mississippiensis", "Reptilia"],
             'XP_009329004.1':["Pygoscelis_adeliae.svg", "Penguin", "Pygoscelis adeliae", "Aves"],
             'XP_013045797.1':["Anser_cygnoides_domesticus.svg", "Goose", "Anser cygnoides domesticus", "Aves"],
             'XP_420532.3':["Gallus_gallus.svg", "Chicken", "Gallus gallus", "Aves"],
             'XP_009486688.1':["Pelecanus_crispus.svg", "Pelican", "Pelecanus crispus", "Aves"],
             'XP_008490178.1':["Calypte_anna.svg", "Hummingbird", "Calypte anna", "Aves"],
             'XP_009564005.1':["Cuculus_canorus.svg", "Cuckoo", "Cuculus canorus", "Aves"],
             'XP_002189592.1':["Taeniopygia_guttata.svg", "Zebra finch", "Taeniopygia guttata", "Aves"],
             'XP_003415871.1':["Loxodonta_africana.svg", "African elephant", "Loxodonta africana", "Mammalia"],
             'XP_540047.2':["Canis_familiaris.svg", "Dog", "Canis lupus familiaris", "Mammalia"],
             'XP_526740.1':["Pan_troglodytes.svg", "Chimpanzee", "Pan troglodytes troglodytes", "Mammalia"],
             'NP_005420.1':["Homo_sapiens.svg", "Human", "Homo_sapiens", "Mammalia"],
             'NP_776913.1':["Bos_taurus.svg", "Cattle", "Bos taurus", "Mammalia"],
             'XP_019777186.1':["Tursiops_truncatus.svg", "Dolphin", "Tursiops truncatus", "Mammalia"],
             'XP_004280970.1':["Orcinus_orca.svg", "Orca", "Orcinus orca", "Mammalia"],
             'NP_033532.1':["Mus_musculus.svg", "Mouse", "Mus musculus", "Mammalia"],
             'NP_446105.1':["Rattus_norvegicus.svg", "Rat", "Rattus norvegicus", "Mammalia"],
             'XP_007496150.2':["Monodelphis_domestica.svg", "Opossum", "Monodelphis domestica", "Mammalia"],
             'XP_004465018.1':["Dasypus_novemcinctus.svg", "Armadillo", "Dasypus novemcinctus", "Mammalia"],
             'XP_002709527.1':["Oryctolagus_cuniculus.svg", "Rabbit", "Oryctolagus cuniculus", "Mammalia"],
             'XP_017897598.1':["Capra_hircus.svg", "Goat", "Capra hircus", "Mammalia"]}
    align()

if __name__ == '__main__':
    run()
