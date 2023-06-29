#!/usr/bin/env python3

import sys
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

align = AlignIO.read(sys.argv[1], "phylip")
cds = AlignIO.read(sys.argv[2], "fasta")
outname = sys.argv[1].split("/")[-1].split(".")[0]

#codons = int(len(align[1].seq)/3)
for l in align:
    if l.id == "hg":
        aseq = str(l.seq)
#for c in range(1,codons):
#    iset = []
#    oset = []
true_codon = 0
for l in cds:
    if l.id == "hg":
        hseq = str(l.seq)
        nogap = hseq.replace("-","")
        alignments = pairwise2.align.globalms(nogap, aseq,2, -1, -.5, -.1)
        nogap_aln = str(alignments[0][0])
        aseq_aln = str(alignments[0][1])
        #print(format_alignment(*alignments[0]))

        for c in range(1,int(len(aseq_aln)/3)+1):
            bases = aseq_aln[c*3-3:c*3]
            if bases != "---":
                true_codon += 1
            print(*[outname,c,true_codon,bases],sep="\t")
            





