import sys

with open(sys.argv[1]) as orthogroups:
    taxanames = next(orthogroups)
    taxanames = taxanames.strip()
    taxanames = taxanames.split("\t")[1:]
    headerline = ["Desc","Family ID"] + taxanames
    print(*headerline, sep = '\t')
    for ortholog in orthogroups:
        ortholog = ortholog.strip()
        ortholog = ortholog.split("\t")
        name = ortholog[0]
        count = 0
        missc = 0
        countlist = []
        for genes in ortholog[1:]:
            taxon = taxanames[count]
            # strip all whitespace
            if not genes:
                countlist.append(0)
                missc += 1
            else:
                genes =  "".join(genes.split())
                genes = genes.split(",")
                countlist.append(len(genes))
            count += 1
        #outlist = ["(null)",name] 
        outlist = ["(null)",name] + countlist
        print(*outlist,sep="\t")
