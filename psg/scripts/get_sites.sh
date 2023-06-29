OG="$1"
# Num meme sites
msites=`grep " Yes, p =" ${OG}.meme.out| tr -s ' '| cut -d' ' -f2,17| wc -l`
# Mean pval
if [[ $msites -gt 0 ]]
    then
    mallsites=`paste <(grep " Yes, p =" ${OG}.meme.out| tr -s ' '| cut -d' ' -f2| sed 's/^/meme\t/') <(grep " Yes, p =" ${OG}.meme.out|  sed 's/.*Yes, p =  //'|sed 's/ .*//')`
    mpval=`grep " Yes, p =" ${OG}.meme.out|  sed 's/.*Yes, p =  //'|sed 's/ .*//'| awk '{sum+=$1}END{print sum /NR}'`
    else
    mpval="NA"
    mallsites=`printf "meme\tNA\tNA\NA\n"`
fi
# Num busted sites
bsites=`python scripts/parse_busted.py ${OG}_busted.json| wc -l`
if [[ $bsites -gt 0 ]]
    then
    ballsites=`python scripts/parse_busted.py ${OG}_busted.json | sed 's/^/busted\t/'`
    ber=`python scripts/parse_busted.py ${OG}_busted.json| cut -d' ' -f2| awk '{sum+=$1}END{print sum /NR}'`
    else
    ber="NA"
    ballsites=`printf "busted\tNA\tNA\NA\n"`
fi
# Background
bg=`grep -A 3 "Improving branch lengths" ${OG}.busted.out| tail -2| sed 's/.* //'| head -1`
# Foreground
fg=`grep -A 3 "Improving branch lengths" ${OG}.busted.out| tail -2| sed 's/.* //'| tail -1`

echo "$mallsites"| sed 's/ /\t/g'| sed "s/^/${1}\t/"
echo "$ballsites"| sed 's/ /\t/g' | sed "s/^/${1}\t/"
#printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$OG" "$msites" "$mpval" "$bsites" "$ber" "$bg" "$fg"
