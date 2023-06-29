#!/usr/bin/env python3

import sys
import os
from os import path
import itertools
import argparse
import csv
from Bio import AlignIO
from Bio.Seq import Seq


def count_codons(algin):
    for counter,record in enumerate(align):
        tempRecordSeq = list(record.seq)
        print(record.id, "Length before : ",len(tempRecordSeq))

def consecutive_ranges(integer_list):
    '''Take in list of integers and return list of tuple 
    pairs with the range of consecutive integers'''
    integer_list = list(map(int, integer_list))
    integer_list = sorted(set(integer_list))
    gaps = [[s, e] for s, e in zip(integer_list, integer_list[1:]) if s+1 < e]
    edges = iter(integer_list[:1] + sum(gaps, []) + integer_list[-1:])
    return list(zip(edges, edges))

def guidance2list(infile):
    '''Convert Seqs.Orig_DNA.fas.FIXED.MSA.MAFFT.Removed_Col file
    to a list of 1-based positions with low scores'''
    rem_list = []
    with open(infile) as removed:
        for line in removed:
            line = line.strip().split()
            #line = line.split()
            position = line[2]
            rem_list.append(position)
    return rem_list

def mask_columns(alignment,targetrange):
    '''Replaces all bases at a list of positions in an
    alignment with N character'''
    # get consecutive masked ranges
    # loop over masked ranges
    for target in targetrange:
        for counter,record in enumerate(alignment):
            tempRecordSeq = list(record.seq)
            record_idx = counter
            targetlen = (target[1] - target[0]) + 1 
            numn = ['N'] * targetlen
            tempRecordSeq = tempRecordSeq[:target[0]-1] + numn + tempRecordSeq[target[1]:]
            alignment[record_idx].seq = Seq(''.join(tempRecordSeq))
    return alignment


def mask_indel(alignment, taxon, indelrange, numneighbor):
    '''Masks codon for a record in a multiple sequence alignment
    object by replacing the codon with NNN. The codon is passed as 
    an integer start position. A specified n neighboring codons
    are also masked upstream and downstream. Returns masked 
    alignment object'''
    for counter,record in enumerate(alignment):
        if taxon in record.id:
            tempRecordSeq = list(record.seq)
            #print("Lengh before : ",len(tempRecordSeq))
            record_idx = counter
    indellen = indelrange[1] - indelrange[0] + 1
    numneighborl = numneighbor * 3
    numneighborr = numneighbor * 3
    if indelrange[0] - numneighborl < 0:
        numneighborl = indelrange[0]
    if indelrange[1] + numneighborr > len(tempRecordSeq):
        print(len(tempRecordSeq),indelrange[1])
        numneighborr = len(tempRecordSeq) - (indelrange[1] + 1)
    #print(indellen,numneighborl,numneighborr)
    numn = ['N'] * (indellen + numneighborl + numneighborr)
    #print(len(numn))
    #print(tempRecordSeq[indelrange[0]-numneighborl:indelrange[1]+1+numneighborr])
    tempRecordSeq = tempRecordSeq[:indelrange[0]-numneighborl] + numn + tempRecordSeq[indelrange[1]+1+numneighborr:]
    #print(tempRecordSeq[indelrange[0]-numneighborl:indelrange[1]+1+numneighborr])
    alignment[record_idx].seq = Seq(''.join(tempRecordSeq))
    #print("Length after : ",len(tempRecordSeq))
    return alignment
       
def find_proximal(querynum,numlist):
    '''Return closest number to querynum in unsorted list of numbers '''
    closest = min(numlist, key=lambda x:abs(x-querynum))
    return closest
       
def classify_indels(indeldict,juncdict,alignment):
    '''Takes in dict of gapped codons for alignment
    indices corresponding to taxa. Outputs a dict with a nested
    list with each deletion's features [del_start, del_stop, attr1, attr2] 
    for each taxon. Requires alignment to get full taxon list'''
    # list of all taxa
    alltaxa = []
    for taxon in alignment:
        alltaxa.append(taxon.id)
    if 'hg' in alltaxa:
        reference = 'hg'
    else:
        reference = alltaxa[0]
    juncproxdict = {}
    for target,indels in indeldict.items():
        inlist =  [list(elem) for elem in indels[1]]
        #concatenate deletions and insertions
        targetindels = indels[0] + inlist
        indelswithjunc = []
        for indel in targetindels:
            try:
                indelstart = indel[0]
                indelstop = indel[1]
                #prox1 = abs(find_proximal(indelstart,juncdict[target]) - indelstart)
                prox1 = find_proximal(indelstart,juncdict[target])
                # sort to always get range in ascending order for slicing sequence
                prox1list = sorted([indelstart,prox1])
                # get real distance after removing gaps
                prox1dist = 0
                for taxon in align:
                    if taxon.id == target:
                        prox1seq = list(taxon.seq)[prox1list[0]:prox1list[1]]
                        for base1 in prox1seq:
                            if base1 not in ["N","-"]:
                                prox1dist += 1
                prox2 = find_proximal(indelstop,juncdict[target])
                prox2list = sorted([indelstop,prox2])
                # get real distance after removing gaps
                prox2dist = 0
                for taxon in align:
                    if taxon.id == target:
                        prox2seq = list(taxon.seq)[prox2list[0]:prox2list[1]]
                        for base2 in prox2seq:
                            if base2 not in ["N","-"]:
                                prox2dist += 1
                               
                #prox2 = abs(find_proximal(indelstop,juncdict[target]) - indelstart)
                closest = min(prox1dist,prox2dist)
            except:
                # single exon genes have no closest junction
                closest = 100
            indelinfo = [indelstart,indelstop,closest]
            indelswithjunc.append(indelinfo)
        juncproxdict[target] = indelswithjunc
    return juncproxdict 
    ## binary deletion attributes
    # 1 overlaps_exon_boundary
    # 3 species_specific_deletion
    # 4 species_specific_insertion
    # exon boundary overlap is only applicable to species_specific indels 
    # return two dicts
    # one with all information for a VCF table
    # the other just with single-species indels that need flank masking    
def find_indel(alignment):
    '''Take in alignment object and return dictionary with
    taxa as keys and a nested list of 0-based gap start and stop
    positions as values.'''
    indeldict = {}
    for taxon in alignment:
        gaplist = []
        tempRecordSeq = list(taxon.seq)
        for idx,base in enumerate(tempRecordSeq):
            # Account for first and last positions being gaps
            if idx == 0 or idx == (len(tempRecordSeq) - 1):
                if "-" in base:
                    gaplist.append(idx)
            # start of gap    
            elif "-" in base and tempRecordSeq[idx-1] != base:
                gaplist.append(idx)
            # end of gap
            elif "-" in base and tempRecordSeq[idx+1] != base:
                gaplist.append(idx)
        i=0
        pairlist=[]
        while i<len(gaplist):
            pairlist.append(gaplist[i:i+2])
            i+=2
        indeldict[taxon.id] = [pairlist]
    # Next we need to add unique insertions
    totalcodons = len(align[0].seq)
    totalrows = len(align)
    n = 0
    indict = {}
    while n < totalcodons:
        # change to also count N characters
        c1 = align[:,n].count("-")
        c2 = align[:,n+1].count("-") 
        c3 = align[:,n+2].count("-") 
        gapc = c1 + c2 + c3
        # codon with unique insertion
        if gapc == (totalrows * 3) - 3:
            coords = [n,n+1,n+2]
            for index,base in enumerate(align[:,n]):
                if base != "-":
                    taxon = align[index].id
                    if taxon in indict:
                        for coord in coords:
                            indict[taxon].append(coord)
                    else:
                        indict[taxon] = coords
        n += 3 
    for sp in indict:
        coordlist = indict[sp]
        # merge consecutive unique inserts
        insrange = consecutive_ranges(coordlist)
        if sp in indeldict:
            indeldict[sp].append(insrange)
        else:
            indeldict[sp] = [[],insrange]
    # If there are no indels, add empty list   
    for key,value in indeldict.items():
        if len(value) == 1:
            value.append([])
    return indeldict

def unique_indel(indeldict,alignment,targets):
    ''' Take in dictionary with nested list values and 
    return dictionary with consecutive numbers as keys. 
    The value is a nested list of previous keys (taxa) and a single
    indel start stop position shared by the taxa.'''
    uniquedict = {}
    for target in targets:
        uniqdels = []
        uniqins = []
        alldels =  [] 
        allins = []
        for taxon,value in indeldict.items():
            if taxon != target:
                # If taxon has deletions
                for dels in value[0]:
                    alldels.append(dels)
                # If taxon has insertions
                if len(value)>1:
                    for ins in value[1]:
                        allins.append(ins)
        try:
            targetdels = indeldict[target][0]
            targetins = indeldict[target][1]
            for candel in targetdels:
                if candel not in alldels:
                    uniqdels.append(candel)
            for canins in targetins:
                if canins not in allins:
                    uniqins.append(canins)
        except:
            pass 
        uniquedict[target] = [uniqdels,uniqins]
    return uniquedict

def Sort(sub_li): 
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using first element of  
    # sublist lambda has been used 
    sub_li.sort(key = lambda x: x[0]) 
    return sub_li 

def gff_to_cds_coords(ingff):
    '''Take path to gff file with CDS sequence of a single 
    taxon and output a list of 0-based coordinates of the
    codon start position before each exon junction. The position
    0 is at the start codon of the gene.'''
    # make nested list of col 4 and 5
    # sort nested list with Sort()
    # start and stop pos list
    sslist = []
    # list of lengths
    llist = []
    # junction list
    jlist = []
    with open(ingff,'r') as cdsgff:
        for line in cdsgff:
            line = line.strip()
            if not line.startswith("#"):
                line = line.split("\t")
                if line[2] == "CDS":
                    cdspos = [int(line[3]),int(line[4])]
                    # can be + or - strand but assumes order is appropriate
                    sslist.append(cdspos)
    for coordpair in sslist:
        # do not add 1 because 0-based
        # getting the true number of bases requires + 1
        cdslen = coordpair[1] - coordpair[0]
        llist.append(cdslen)
    startpos = 0
    for length in llist:
        juncpos = startpos + length
        jlist.append(juncpos)
        startpos = juncpos
    return jlist

def lift_cds_coords_to_aln(alignment,targets):
    '''Take an alignment object and returns a dict of junction
    coordinate lists per species lifted over to the gapped
    alignment. This requires the workdir to contain gff
    files named as <taxon>.gff with CDS positions for the gene.'''
    juncdict = {}
    for taxon in alignment:
        if taxon.id in targets:
            # lifted over junctions
            ljunc = []
            taxonid = taxon.id
        # my custom headers have species name last after __ separator
            #cleanid = taxonid.split("__")[-1]
            gffname = taxonid + ".gff"
            juncout = gff_to_cds_coords(gffname)
            tempRecordSeq = list(taxon.seq)
            # ignore gaps to lift over coords
            for junction in juncout:
                ungap_pos = 0
                for idx, pos in enumerate(tempRecordSeq,start=1):
                    if "-" not in pos:
                        ungap_pos += 1
                    # hit junction position in alignment sequence
                    if ungap_pos == int(junction):  
                        ljunc.append(idx)
                        break
            juncdict[taxonid] = ljunc
    print(juncdict)
    return juncdict

def filter_short_block(align,targets):
    recnum = 0
    for record in align:
        if record.id in targets:
            baselist = []
            tempRecordSeq = list(record.seq)
            for idx,base in enumerate(tempRecordSeq):
                if base not in ["N","-"]:
                    baselist.append(idx)
            alnregions = consecutive_ranges(baselist)
            for region in alnregions:
                alnlen = (region[1] - region[0]) + 1 
                if alnlen <= 30:
                    numn = ['N'] * alnlen
                    tempRecordSeq = tempRecordSeq[:region[0]] + numn + tempRecordSeq[region[1]+1:]
                    align[recnum].seq = Seq(''.join(tempRecordSeq))
        recnum += 1
    return align

def rowfilter(align,targets):
    totalcodons = len(align[0].seq)
    recnum = 0
    for record in align:
        if record.id in targets:
            tempRecordSeq = list(record.seq)
            missingcount = tempRecordSeq.count("N") + tempRecordSeq.count("-")
            fracmissing = float(missingcount) / float(totalcodons) 
            if fracmissing > 0.5:
                print("Fully masking ", record.id, " due to high missingness")
                numn = ['N'] * totalcodons
                tempRecordSeq = numn
                align[recnum].seq = Seq(''.join(tempRecordSeq))
        recnum += 1
    return align

def main(align):
    # get exon junctions
    totalcodons = len(align[0].seq)
    totalrows = len(align)
    # init empty alignment
    clean = align[:,0:0]
    n = 0
    mask_codon_dict = {}
    print("Parsing ", totalcodons / 3, " codons for ", totalrows, " taxa...") 
    while n < totalcodons:
        # change to also count N characters
        c1 = align[:,n].count("-") + align[:,n].count("N") 
        c2 = align[:,n+1].count("-") + align[:,n+1].count("N") 
        c3 = align[:,n+2].count("-") + align[:,n+2].count("N") 
        gapc = c1 + c2 + c3
        # gaps in codon aln
        # if statements will fail for python 2 due to integer division
        # over half of rows gapped
        #if (c1 / totalrows) > 0.2 or (c3 / totalrows) > 0.2 or (c3 / totalrows) > 0.2:
        if int(c1) > 3:
            recnum = 0
            for record in align:
                tempRecordSeq = list(record.seq)
                codon = record.seq[n:n+3]
                # mask codon
                tempRecordSeq = tempRecordSeq[:n] +['N','N','N'] + tempRecordSeq[n+3:]
                align[recnum].seq = Seq(''.join(tempRecordSeq))
                recnum += 1
        # less than half of rows gapped    
        #elif (c1 / totalrows) < 0.5 or (c3 / totalrows) < 0.5 or (c3 / totalrows) < 0.5:
        else:
            recnum = 0
            # loop over codon in each row
            for record in align:
                tempRecordSeq = list(record.seq)
                codon = record.seq[n:n+3]
                # codons contains gap
                if any(x in codon for x in ('-')):
                    #print(codon)
                    #print(recnum, n, num_neighbor)
                    #align = mask_record_codon(align, recnum, n, num_neighbor)
                    # mask codon
                    tempRecordSeq = tempRecordSeq[:n] + ['N', 'N', 'N'] + tempRecordSeq[n+3:]
                    align[recnum].seq = Seq(''.join(tempRecordSeq))
                    if recnum in mask_codon_dict:
                        mask_codon_dict[recnum].append(n)
                    else:
                        mask_codon_dict[recnum] = [n]
                recnum += 1
            # incrementally add good codons to clean alignment
            clean = clean + align[:,n:n+3]
        # move forward by one 3bp codon
        n += 3
    # Function to remove codon pos from dict that has neighbors masked
    # retain only start and stop codons of gaps
    # mask

    print("Analysis complete. Exporting masked alignments...")
    # export masked mfa with only good codons
    outclean = args.out + ".clean.mfa"
    with open(outclean, "w") as handle:
        AlignIO.write(clean, handle, "fasta")
    # export masked mfa with all codons
    #outmask = args.out + ".mask.mfa"
    #with open(outmask, "w") as handle:
    #    count = AlignIO.write(align, handle, "fasta")


my_parser = argparse.ArgumentParser(description='Mask indel positions in a multifasta alignment')
my_parser.add_argument('infasta',
                       metavar='infasta',
                       type=str,
                       help='input fasta file')
my_parser.add_argument('out',
                       metavar='out',
                       type=str,
                       help='prefix for output files')
my_parser.add_argument('-m',
                       '--mask_size',
                       #type=int,
                       action='store_true',
                       default=1,
                       help='number of neighboring codons to mask after indel masking (default:1)')
my_parser.add_argument('-b',
                      '--minimum_block', 
                      action='store_true',
                      help='Mask aligned blocks of bases of less than 30bp in target species (default:false)')
my_parser.add_argument('-r',
                      '--row_filter', 
                      action='store_true',
                      help='fully mask sequence of species with over 50% missing characters (default:false)')


args = my_parser.parse_args()
print("Mask size set to ", args.mask_size)
fasta = args.infasta
if not path.exists(fasta):
    print('The input fasta file does not exist')
    sys.exit()
targets = ['jam', 'meso']
align = AlignIO.read(fasta, "fasta")
print("Fasta import complete")
indeldict = find_indel(align)
uindel = unique_indel(indeldict,align,targets)
juncdict = lift_cds_coords_to_aln(align,targets)

results = classify_indels(uindel,juncdict,align)
#uniqindels = results[0]
#allindels = results[1]
newalign = align
count_codons(newalign)
for taxon,uniqindels in results.items():
    for uniqindel in uniqindels:
        startpos = uniqindel[0]
        stoppos = uniqindel[1]
        juncprox = uniqindel[2]
        mask_size = args.mask_size
        if int(juncprox) <5:
            mask_size = 3
        indelrange = [startpos,stoppos]
        newalign = mask_indel(newalign,taxon,indelrange,mask_size)
# Mask guidance columns
guidefile = "Seqs.Orig_DNA.fas.FIXED.MSA.MAFFT.Removed_Col"
pguidefile = "Seqs.Orig_DNA.fas.FIXED.MSA.PRANK.Removed_Col"
#count_codons(newalign)
try:
    try:
       if os.stat(guidefile).st_size != 0:
           targetlist = guidance2list(guidefile)
           targetrange = consecutive_ranges(targetlist)
           guide_aln = mask_columns(newalign,targetrange)
       else:
           guide_aln = newalign
    except:
       if os.stat(pguidefile).st_size != 0:
           targetlist = guidance2list(pguidefile)
           targetrange = consecutive_ranges(targetlist)
           guide_aln = mask_columns(newalign,targetrange)
       else:
           guide_aln = newalign
except:
    guide_aln = newalign
print(results)
count_codons(guide_aln)
if args.minimum_block:
    print("Masking short blocks...")
    guide_aln = filter_short_block(guide_aln,targets)
if args.row_filter:
    print("Masking rows with >0.5 missing characters...")
    guide_aln = rowfilter(guide_aln,targets)

main(guide_aln)

print("All done!")

# Begin task
#main(align,numneighbor,uniqindels)
