name="$1"
# PARSE ABSREL
ABSREL="${name}.absrel.out"
# Median omega
omega=`grep "Branch-level non-synonymous/synonymous rate ratio" $ABSREL| sed 's/.*distribution has median  //' | sed 's/,.*//'`
#Codon number
codons=`grep "Loaded a multiple sequence alignment with" $ABSREL| sed 's/.*sequences, \*\*//'| sed 's/\*\* codons.*//'`

about=`grep -H -E 'Node10' $ABSREL| grep "|" | tail -1|sed 's/ //g'| tr '|' '\t'|sed 's@^./@@'| sed 's@/.*:@@'| sed 's/.absrel.out//'| cut -f1,2,6 | sed 's/Node10/bat/'| sed 's/://'`
cat branchsite_delta_likelihood.txt|grep -v big | while read l; do
    og=`echo "$l"| cut -f1`
    lineage=`echo "$l"| cut -f2`
    val=`echo "$l"| cut -f3| sed 's/-//'`
    pval=`chi2 1 $val| sed 's/.*prob = //'| sed 's/ .*//'| sed '/^$/d'`
    absrelpval=`echo "$about"| grep $lineage | cut -f3`

    if [[ -z "${absrelpval}" ]]
    then
        absrelpval="NA"
    fi
    if [[ -z "$codons" ]]
    then
        codons="NA"
    fi
    if [[ -z "$omega" ]]
    then
        omega="NA"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$og" "$lineage" "$pval"  "$absrelpval" "$codons" "$omega"

done


