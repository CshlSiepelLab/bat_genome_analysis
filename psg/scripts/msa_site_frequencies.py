#!/usr/bin/env python3

import sys
from Bio import AlignIO
ingroup = ["jam","meso","HLpipKuh2","HLmyoMyo6","mydav","mybra","drot","HLphyDis3","HLrouAeg4","HLrhiFer5","efus","HLmolMol2","myluc","mnat","palec"]
outgroup = ["mm","hg","ecab","sscrofa","canfam"]
align = AlignIO.read(sys.argv[1], "phylip")
outname = sys.argv[1].split("/")[-1].split(".")[0]

def translate(seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein


codons = int(len(align[1].seq)/3)
for c in range(1,codons):
    iset = []
    oset = []
    for l in align:
        name = l.id
        cds = str(l.seq[c*3-3:c*3])
        if "N" in cds or "n" in cds:
            aa = "X"
        else:
            aa = translate(cds)
        if name in outgroup:
            oset.append(aa)
        elif name in ingroup:
            iset.append(aa)
    # major allele freq
    major_allele = max(set(iset), key=iset.count)
    major_allele_o = max(set(oset), key=oset.count)
    maf = iset.count(major_allele) / len(iset)
    mafo = oset.count(major_allele_o) / len(oset)
    # intersecting residues
    intersect_frac = len([value for value in iset if value in oset]) / len(iset)
    outline = [outname,c,maf,mafo,intersect_frac]
    print(*outline,sep="\t")

