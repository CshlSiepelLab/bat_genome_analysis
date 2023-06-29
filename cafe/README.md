
Calibration points from timetree.org: data/mcmctree_calibration.nwk

Extract 1st and 2nd codon positions from alignment.

split_msa_by_codon.py
concat.py



first calc gradient and hessian of the branch lengths with "usedata = 0" in the control file, then rename out.BV to in.BV and set "usedata = 2" and run the final analysis.

mcmctree

The final tree with branch lengths is output in nexus format (see data/timetree.tre). I converted this to a plainer newick tree to use with CAFE (data/bat_tree.nwk) 

python of2cafe.py data/Orthogroups.tsv > Orthogroups.cafe.tsv


cafe5 -i data/Orthogroups.min3.max90dif.tsv -t data/bat_tree.nwk &> 1lambda.out
cafe5 -i data/Orthogroups.min3.max90dif.tsv -t data/bat_tree.nwk -y data/bat_tree_lambda2.txt &> 2lambda.out
cafe5 -i data/Orthogroups.min3.max90dif.tsv -t data/bat_tree.nwk -y data/bat_tree_lambda.txt &> 3lambda.out

# Use full set of orthogroups for final analysis
cafe5 -i data/Orthogroups.cafe.txt -t data/bat_tree.nwk -y data/bat_tree_lambda.txt -m 0.00054549662537203,0.00081172003147854,0.0017662770480184 -o results -e data/Base_error_model.txt

Resource
http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf
trimal -gapthreshold 0.9 -in all_concat.fasta -out all_codons12_concat.py -phylip
paste <(cut -f1 -d' ' Orthogroups_sco_min12bat_min3out_max1multi.tsv ) <(cut -d' ' -f2 Orthogroups_sco_min12bat_min3out_max1multi.tsv | while read l; do echo "$l" |tr ',' '\n'| wc -l;done)| awk '$2=="20"' > og_allsp.txt
cut -f1 og_allsp.txt | while read l; do ls mfa/$l/*fas;done > og_allsp_complete.list
while read l; do cat $l >>all.fasta;done<og_allsp_complete.list
python concat.py > all_concat.fa
python split_msa_by_codon.py all_concat.fa
python split_msa_by_codon.py all_concat.fa
cat codons*mfa > codons_all.mfa
python concat.py > all_concat.fasta 
trimal -gapthreshold 0.9 -in all_concat.fasta -out all_codons12_concat.py -phylip
