OG="$1"
sp=`grep '>' ${OG}.best.fas | sed 's/>//'| sed 's/meso/meso{Foreground}/'| sed 's/jam/jam{Foreground}/'| tr '\n' ',' | sed 's/,$//'`
tree_doctor -N --prune-all-but $sp -n data/absrel.tre > ${OG}.tre
touch meso.gff
touch jam.gff
python scripts/clean_msa.py -r -b ${OG}.best.fas ${OG} > clean.out
touch dummy.txt
python scripts/add_missing_taxa.py -p ${OG}.clean.mfa dummy.txt ${OG}
hyphy aBSREL --tree ${OG}.tre --alignment ${OG}.alltaxa.phy --branches Foreground --output ${OG}_ABSREL.json > ${OG}.absrel.out
FILE=errors.log
if test -f "$FILE"; then
    echo "aBSREL failed!"
fi
