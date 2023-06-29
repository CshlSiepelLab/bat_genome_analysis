#!/usr/bin/env python3

import sys
from os import path
import argparse

my_parser = argparse.ArgumentParser(description='Mask indel positions in a multifasta alignment')
my_parser.add_argument('intable',
                       metavar='intable',
                       type=str,
                       help='input Orthogroups.tsv file from Orthofinder')
my_parser.add_argument('maxmissing',
                       metavar='maxmissing',
                       type=int,
                       help='maximum missing taxa from orthogroup')
my_parser.add_argument('blacklist',
                       metavar='blacklist',
                       type=str,
                       help='File with newline separator of taxa that are never allowed to be missing')
#my_parser.add_argument('dupelist',
#                       metavar='dupelist',
#                       type=str,
#                       help='File with newline separator of taxa that are allowed to be duplicated')

args = my_parser.parse_args()
intable = args.intable
maxm = args.maxmissing
blist = args.blacklist
#dlist = args.dupelist


ingroup = set(["HLmolMol2","HLmyoMyo6","HLphyDis3","HLpipKuh2","HLrhiFer5","HLrouAeg4","drot","efus","jam","meso","mnat","mybra","mydav","myluc","palec"])
outgroup = set(["canfam","ecab","hg","mm","sscrofa"])

if not path.exists(intable):
    print('The input ortholog table does not exist')
    sys.exit()
if not path.exists(blist):
    print('The input blacklist does not exist')
    sys.exit()
bl = set()
#dl = []
with open(blist,'r') as bl_taxa:
    for line in bl_taxa:
        line = line.strip()
        bl.add(line)    
#with open(dlist,'r') as dl_taxa:
#    for line in dl_taxa:
#        line = line.strip()
#        dl.append(line)    
  
with open(intable) as orthogroups:
    taxanames = next(orthogroups)
    taxanames = taxanames.strip()
    taxanames = taxanames.split("\t")[1:]
    for ortholog in orthogroups:
        ortholog = ortholog.strip()
        ortholog = ortholog.split("\t")
        name = ortholog[0]
        tmpdict = {}
        cleantaxa = set()
        misstaxa = set()
        multaxa = set()
        count = 0
        multicopy = False
        uniqdup = False
        for genes in ortholog[1:]:
            taxon = taxanames[count]
            # strip all whitespace
            genes =  "".join(genes.split())
            genes = genes.split(",")
            if genes != [''] and len(genes) == 1:
                tmpdict[taxon] = genes
                cleantaxa.add(taxon)
            elif genes == ['']:
                misstaxa.add(taxon)
            elif  len(genes) > 1: #and taxon not in dl:
                multaxa.add(taxon)
            #elif len(genes) > 1 and taxon in dl:
            #    cleantaxa.add(taxon)
            #if len(genes) > 1 and taxon in dl:
            #    uniqdup = True
            count += 1
               
        missing_taxa = len(taxanames) - len(cleantaxa)
        ingroup_present = len(ingroup.intersection(cleantaxa))
        outgroup_present = len(outgroup.intersection(cleantaxa))
        #if len(coretaxa_present) == len(bl) and missing_taxa <= maxm and not multicopy and uniqdup:

        if ingroup_present >= 12 and outgroup_present >= 3 and len(multaxa) <= 1:
            genelist = []
            for sp,gene in tmpdict.items():
                if sp not in list(multaxa):
                    genelist.append(gene[0])
                
            genelist = ",".join(genelist)
            print(name,genelist)
