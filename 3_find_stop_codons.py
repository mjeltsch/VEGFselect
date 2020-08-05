#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script tries to find stop codons in aligned mRNA sequences

import sys
from Bio import SeqIO

def run():
    stop_codons = ["TGA", "TAG", "TAA"]
    is_stop_codon = lambda x : any(x[i:i+3] in stop_codons for i in range(0,len(x),3))
    tmp_files = sys.argv[1:]
    for file in tmp_files:
        records = list(SeqIO.parse(file, "fasta"))
        for seq_record in records:
            #print(str(seq_record.seq))
            if is_stop_codon(str(seq_record.seq)):
                print('Sequence {0} has a stop codon.'.format(seq_record.description))

if __name__ == '__main__':
    run()
