import sys
import os.path
import itertools
from Bio import AlignIO
from collections import defaultdict

csdict = defaultdict(list)
with open(sys.argv[2]) as inbed:
    for line in inbed:
        line = line.strip()
        line = line.split("\t")
        chrom = line[0]
        start = int(line[1])
        csdict[chrom].append(start)

inmaf = AlignIO.parse(sys.argv[1], "maf")

for multiple_alignment in inmaf:
    # require human alignment
    for record in multiple_alignment:
        species = str(record.id).split(".")[0]
        chrom = str(record.id).split(".")[1]
        start = int(record.annotations["start"])
        if species == "hg" and chrom in csdict:
            hgchrom = str(record.id).split(".")[1]
            hgstart = str(record.annotations["start"])
            if start in csdict[chrom]:
                # print record as fasta
                outfile = hgchrom + "__" + hgstart + ".fa"
                if not os.path.isfile(outfile):
                    for record in multiple_alignment:
                        species = str(record.id).split(".")[0]
                        if "Anc" not in species:
                            seq = record.seq
                            chrom = str(record.id).split(".")[1]
                            start = str(record.annotations["start"])
                            strand = str(record.annotations["strand"])
                            header = ">"+species + "__" + chrom +"__" + start + "__" + strand
                            with open(outfile,'a') as out:
                                out.write(header+"\n")
                                out.write(str(seq)+"\n")
    
