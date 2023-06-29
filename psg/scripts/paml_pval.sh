cat $1 | while read l; do
    og=`echo "$l"| cut -f1`
    lineage=`echo "$l"| cut -f2`
    val=`echo "$l"| cut -f3| sed 's/-//'`
    pval=`chi2 1 $val| sed 's/.*prob = //'| sed 's/ .*//'| sed '/^$/d'`
    if [ -z "$val" ]
    then
        echo "MISSING"
    else
        echo "$og $lineage $pval"
    fi
done
