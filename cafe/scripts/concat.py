#!/usr/bin/python

from Bio import SeqIO
from collections import defaultdict

sequence_map = defaultdict(str)

for sequence in SeqIO.parse('all.fasta', "fasta"):
  sequence_map[sequence.name] += str(sequence.seq)

for key,value in sorted(sequence_map.items()):
  print('>' + key)
  print(sequence_map[key])
