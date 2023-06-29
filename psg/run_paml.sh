OG="$1"
trimal -in ${OG}.alltaxa.phy -phylip_paml > ${OG}.alltaxa.paml.phy
sp=`grep '>' ${OG}.best.fas | sed 's/>//'| tr '\n' ',' | sed 's/,$//'`
# for each paml tree to test
lineage="bat"
tree_doctor -N --prune-all-but $sp -n data/paml.tre | sed 's/ # 1/#1/' > ${OG}.${lineage}.tre
mkdir altbs_$lineage
cd altbs_$lineage
cp ../data/codeml.ctl .
sed -i "s/PLACEHOLDER/${OG}/g" codeml.ctl
sed -i "s/LINEAGE/$lineage/" codeml.ctl
sed -i "s/TEST/alt/" codeml.ctl
codeml
cd ..
mkdir nullbs_$lineage
cd nullbs_$lineage
cp ../data/codeml.ctl .
sed -i "s/PLACEHOLDER/${OG}/g" codeml.ctl
sed -i "s/LINEAGE/$lineage/g" codeml.ctl
sed -i "s/fix_omega = 0/fix_omega = 1/" codeml.ctl
sed -i "s/TEST/null/" codeml.ctl
codeml
cd ..
altlike=`grep "lnL" ${OG}.${lineage}.bs.alt.paml.out | sed 's/.*-/-/'| sed 's/ .*//'`
nulllike=`grep "lnL" ${OG}.${lineage}.bs.null.paml.out| sed 's/.*-/-/'| sed 's/ .*//'`
delta=`echo "2*($altlike - $nulllike)" | bc -l`
printf "%s\t%s\t%s\n" "${OG}" "$lineage" "$delta" >> branchsite_delta_likelihood.txt
