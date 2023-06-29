#!/usr/bin/env python3

import sys
from os import path
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
my_parser = argparse.ArgumentParser(description='Add dummy sequence of missing taxa to a FASTA alignment')
my_parser.add_argument('infasta',
                       metavar='infasta',
                       type=str,
                       help='input fasta file')
my_parser.add_argument('taxa',
                       metavar='taxa',
                       type=str,
                       help='File with newline separated list of all taxa')
my_parser.add_argument('out',
                       metavar='out',
                       type=str,
                       help='prefix for output file')
my_parser.add_argument('-p', '--phylip', action='store_true', 
    help="export alignment in phylip format")

args = my_parser.parse_args()
fasta = args.infasta
taxa = args.taxa
if not path.exists(fasta):
    print('The input fasta file does not exist')
    sys.exit()
elif not path.exists(taxa):
    print('Input taxa file does not exist')
    sys.exit()

align = AlignIO.read(fasta, "fasta")
print("Fasta import complete")

codon_stop_array = ["TAG", "TGA", "TAA", "UGA", "UAA", "UAG"]

alltaxa = []
with open(taxa,'r') as intaxa:
    for line in intaxa:
        line = line.strip()
        alltaxa.append(line)    
fataxa = []
align_list = []
for record in align:
    fataxa.append(record.id)
    tempRecordSeq = list(record.seq)
    for index in range(0, len(record.seq), 3):
        codon = record.seq[index:index+3]
        if codon in codon_stop_array:
            tempRecordSeq[index:index+3] = 'N','N','N'
    record.seq = Seq("".join(tempRecordSeq))
    a = SeqRecord(record.seq, id=record.id)
    align_list.append(a)
align_out = MultipleSeqAlignment(align_list)

positions = len(align_out[0].seq)
dummyseq = "N" * positions
for taxon in alltaxa:
    if taxon not in fataxa:
        align_out.add_sequence(taxon, dummyseq)

if args.phylip:
    outmiss = args.out + ".alltaxa.phy"
    with open(outmiss, "w") as handle:
        count = AlignIO.write(align_out, handle, "phylip-relaxed")
else:
    outmiss = args.out + ".alltaxa.mfa"
    with open(outmiss, "w") as handle:
        count = AlignIO.write(align_out, handle, "fasta")

