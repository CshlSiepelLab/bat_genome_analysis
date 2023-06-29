#!/usr/bin/env python3

import sys
from Bio import AlignIO
#from Bio.Seq import Seq

def main(align):
    clean = align[:,0:0]
    clean = clean + align[:,0::3]
    return clean

align = AlignIO.read(sys.argv[1], "fasta")
outalign = main(align)
outmsa = "codons1.mfa"
with open(outmsa, "w") as handle:
    AlignIO.write(outalign, handle, "fasta")
