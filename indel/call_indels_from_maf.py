import sys
import itertools
from Bio import AlignIO

##############################
# Configurable parameters
inmaf = sys.argv[1]
# ann species must always be in alignment
# coordinates are based on this species
annspecies = sys.argv[2]
refspecies = sys.argv[3]
ingroup = set(sys.argv[4].split(','))
outgroups = set(sys.argv[5].split(','))
# indel_type can be "del" or "ins"
indel_type = sys.argv[6]
min_alignment_depth = 4
min_alignment_len = 10
# min overlap fraction betweem two deletions to count as the same
min_overlap = 0.9
# Hardcoded deletion calling parameters
#annspecies = "hg"
#refspecies = "hg"
#ingroup = {"phydis","meso","jam","rhifer","drotund","mymy"}
#outgroups = {"hg","mm","sscrofa"}

# Uncomment below lines to output insertions
#refspecies = "mymy"
#outgroups = {"phydis","meso","jam","rhifer","drotund","mymy"} 
#ingroup = {"hg","mm","sscrofa"}
##############################

def callvar(ref,alt):
    # convert to uppercase to avoid false polymorphisms
    ref = ref.upper()
    alt = alt.upper()
    if ref == '-':
        if alt != '-':
            return 'I'  # insertion
    else:
        if alt == '-':
            return 'D'  # deletion
        elif alt != ref:
            return 'P'  # polymorphism
    return 'U'          # unchanged

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def indelpos(rangelist,minlen,maxlen,refseq,refstart):
    # input rangelist example: [(0, 4), (7, 9), (11, 11)]
    outlist = []
    for myrange in rangelist:
        # bases before indel
        if myrange[1] - myrange[0] > minlen and myrange[1] - myrange[0] < maxlen:
            bases = len(refseq[0:myrange[0]])
            gaps = refseq[0:myrange[0]].count("-")
            blen = bases - gaps
            indelbases = len(refseq[myrange[0]:myrange[1]]) - refseq[myrange[0]:myrange[1]].count("-")
            truestart = refstart + blen 
            trueend = truestart + indelbases
            outlist.append((truestart,trueend))
    return outlist
def merge_intervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged





for multiple_alignment in AlignIO.parse(inmaf, "maf"):
    #print("printing a new multiple alignment")
    #print("Alignment length %i" % multiple_alignment.get_alignment_length())
    #print("Alignment depth",multiple_alignment.__len__())
    # if record longer 20bp and has over 4 species aligned
    if multiple_alignment.__len__() >= min_alignment_depth and multiple_alignment.get_alignment_length() >= min_alignment_len:
        spset = set()
        tmpdict = {}
        for record in multiple_alignment:
            # start pos
            start = record.annotations["start"]
            start_id = start
            # ungapped length
            size = record.annotations["size"]
            # strand where 1 is positive and -1 is negative
            strand = record.annotations["strand"]
            src = record.annotations["srcSize"]
            # To Do: check for off-by-one errors
            if strand == -1:
                start = src - (start + size)
                end = src - start
            elif strand == 1:
                end = start + size
            # the end position can be calculated based on the above three values
            #print(record.annotations)
            spname = record.id
            seqid = ".".join(spname.split(".")[1:])
            spname = spname.split(".")[0]
            spset.add(spname)
            seq = str(record.seq)
            rstart = start
            rend = end
            alnlist = [seqid,seq,rstart,rend,strand,start_id,size]
            tmpdict[spname] = alnlist


        aln_out = list(spset.intersection(outgroups))
        aln_in = list(spset.intersection(ingroup))
        al_sp = aln_out + aln_in
        try:
            if annspecies in al_sp:
                deldict = {}
                for a in al_sp:
                    counts = dict.fromkeys('DIPU', 0)
                    indelseq = ""
                    oseq = tmpdict[a][1]
                    rseq = tmpdict[refspecies][1]
                    # Compare outgroup sp sequence to ref
                    for ref, alt in zip(rseq, oseq):
                        indel = callvar(ref, alt)
                        indelseq = indelseq + indel
                        # count variant types
                        counts[indel] += 1
                    # Get positions of all deletions in outgroup
                    delidx = list(ranges(list_duplicates_of(indelseq, 'D')))
                    deldict[a] = delidx
        
                refdel = deldict[aln_out[0]]
        
                # cluster dict
                cdict = {}
                aranges = []
                for key,value in deldict.items():
                    # find largest clusters based on all deletion ranges
                    for item in value:
                        aranges.append(item)
                clusters = merge_intervals(aranges) 
        
                for clust in clusters:
                    cname = str(clust[0]) + "_" + str(clust[1])
                    cdict[cname] = {}
                    for key,value in deldict.items():
                         if key not in cdict[cname]:
                             cdict[cname][key] = []   
                         for deletion in value:
                            isectr = range(max(deletion[0], clust[0]), min(deletion[-1], clust[-1])+1)
                            isectlen = len(isectr)
                            if isectlen > 0:
                                cdict[cname][key].append(deletion)
                                # Add deletion to cluster
                            # remove duplicates
                            #c_overlap = [t for t in (set(tuple(i) for i in c_overlap))]
                        # Check if overlap is from same species due to double alignment
                # In each cluster get unique ranges
                # For each unique range, report ingroup, outgroup perfect/partial matches
                # a perfect match has the same start and end position
                for c in cdict:
                    aranges = []
                    udict = {}
                    for key,value in cdict[c].items():
                        # if species has deletion
                        if value:
                           for v in value:
                               aranges.append(v)
                     # unique ranges
                    udels = [t for t in (set(tuple(i) for i in aranges))]
                    for u in udels:
                        match = []
                        lap = []
                        lapmatch = []
                        for key,value in cdict[c].items():
                            # if species has deletion
                            if value:
                                for v in value:
                                    isectr = range(max(u[0], v[0]), min(u[-1], v[-1])+1)
                                    isectlen = len(isectr)
                                    vlen = len(range(v[0],v[1])) + 1
                                    ulen = len(range(u[0],u[1])) + 1
                                    if isectlen > 0:
                                        lap.append(key)
                                    if isectlen == vlen and vlen == ulen:
                                        match.append(key)
                                    if isectlen / ulen > min_overlap and  isectlen / vlen > min_overlap:
                                        lapmatch.append(key)
                        # print the position of the deletion in hg,mm
                        cstart = tmpdict[annspecies][2]
                        dchrom = tmpdict[annspecies][0]
                        delstart = cstart + u[0]
                        delend = cstart + u[1]
                        dstrand = tmpdict[annspecies][4]

                        start_tag = tmpdict[annspecies][5] 
                        size_tag = tmpdict[annspecies][6]
                        #aln_out
                        #aln_in
                        tot_in = len(list(set(aln_in)))
                        tot_out = len(list(set(aln_out)))
                        in_match = len(list(set(aln_in).intersection(set(match))))
                        in_lapmatch = len(list(set(aln_in).intersection(set(lapmatch))))
                        in_lap = len(list(set(aln_in).intersection(set(lap))))
                        out_match = len(list(set(aln_out).intersection(set(match))))
                        out_lapmatch = len(list(set(aln_out).intersection(set(lapmatch))))
                        out_lap = len(list(set(aln_out).intersection(set(lap))))
                        delrange_out = ",".join(map(str, u))
                        matchsp = ",".join(map(str,match))
                        lapmsp = ",".join(map(str,lapmatch))
                        lapsp = ",".join(map(str,lap))

                        # start_tag : start position of alignment
                        # size_tag : length of alignment
                        # tot_in : number of ingroup species in alignment
                        # in_match : number of ingroup species with matching indel genotype
                        # in_lapmatch : number of ingroup species with overlapping indel and matching genotype 
                        # in_lap: ingroup species with an overlapping indel
                        # tot_out: outgroup species in alignment
                        # out_match: outgroup species with matching indel genotype
                        # out_lapmatch: outgroups species with matching indel genotype and overlapping position
                        # What is most important is having mostly overlaping indels in one group, and no overlapping indels in the other group.
                        # For added stringency, we can require that the ingroup indels have the exact same sequence
                        # Simple output would thus be counts of in and out species, and counts of overlapping indels in each group and counts of matching indels in each group
                        #if int(in_lapmatch) >1 and tot_out> 1 and out_lap == 0:
                        in_cnt = 0
                        out_cnt = 0
                        for i in ingroup:
                            if i in tmpdict:
                                in_cnt += 1
                        for i in outgroups:
                            if i in tmpdict:
                                out_cnt += 1

                        if indel_type == "ins":
                            splist = []
                            for s in tmpdict.keys():
                                if s not in lap:
                                    splist.append(s)
                            splist = ','.join(splist)
                            outline = [annspecies,dchrom,delstart,delstart,delend-delstart+1,indel_type,in_cnt,out_cnt,in_lap,out_lap,in_lapmatch,out_lapmatch,splist]
                        elif indel_type == "del":
                            outline = [annspecies,dchrom,delstart,delend,delend-delstart+1,indel_type,in_cnt,out_cnt,in_lap,out_lap,in_lapmatch,out_lapmatch,lapmsp]
                        print(*outline,sep='\t')
                        # print the alignment for debugging
                        for sp,alseq in tmpdict.items():
                            if "Anc" not in sp:
                                print(sp,alseq[1])
        except:
            pass
