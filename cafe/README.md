# Gene family size expansion and contraction analysis with CAFE5

The [CAFE5](https://github.com/hahnlab/CAFE5) analysis was carried out in four steps.

1. Generate an alignment of 1st and 2nd codon positions from conserved genes
2. Infer an ultrametric time tree with PAML
3. Generate a matrix of gene counts by species from Orthofinder output
4. Infer expanded and contracted gene families with CAFE5


Firstly, alignments of conserved single-copy orthologs (occuring in >11 bat species and >2 outgroup species with <2 species with a multicopy gene) were written to a single file and concatenated into one sequence per species with `scripts/concat.py`. 1st and 2nd codon positions were extracted from the alignments using `scripts/split_msa_by_codon.py`. The resulting file was converted to phylip and sites with <90% of species present were removed using `trimal`.

```
trimal -gapthreshold 0.9 -in codons12.fasta -out codons12_concat.py -phylip
```  

Secondly, calibration points from [timetree](http://timetree.org/) were used to produce the calibrated constraint tree `data/mcmctree_calibration.nwk`. To infer the timetree with PAML, we ran `mcmctree` two times. This program relies on the configuration file `data/mcmctree.ctl`. The initial run calculated the gradient and hessian of the branch lengths with `usedata = 0` in the configuration file, then rename out.BV to `data/in.BV` and set `usedata = 2` and run the final analysis.

```
mcmctree
```

The final tree with branch lengths is output in nexus format (see `data/timetree.tre`). This file can be reformatted to a plainer newick tree to use with CAFE5 (`data/bat_tree.nwk`).

Thirdly, we converted the default Orthofinder output file `data/Orthogroups.tsv`, containing comma-separated lists of genes for each each species in each orthogroup, to a count matrix. 

```
python of2cafe.py data/Orthogroups.tsv > Orthogroups.cafe.tsv
```

For inferring parameters with CAFE5, we generated a filtered version `data/Orthogroups.min3.max90dif.tsv` excluding Orthogroups with fewer than three species or a larger maximimum delta gene count between two species than 90. We then calculated the likelihood of our tree having 1 to 3 independent evolutionary rates (lambda).

```
cafe5 -i data/Orthogroups.min3.max90dif.tsv -t data/bat_tree.nwk &> 1lambda.out
cafe5 -i data/Orthogroups.min3.max90dif.tsv -t data/bat_tree.nwk -y data/bat_tree_lambda2.txt &> 2lambda.out
cafe5 -i data/Orthogroups.min3.max90dif.tsv -t data/bat_tree.nwk -y data/bat_tree_lambda3.txt &> 3lambda.out
```

We found that three independent lambdas provided the best fit to the data. The final CAFE5 analysis was thus executed with three pre-estimated lambdas.

```
cafe5 -i data/Orthogroups.cafe.txt -t data/bat_tree.nwk -y data/bat_tree_lambda3.txt -m 0.00054549662537203,0.00081172003147854,0.0017662770480184 -o results -e data/Base_error_model.txt
```

## Resources

* http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf
* https://github.com/hahnlab/CAFE5

## Dependencies

* PAML
* CAFE5
* biopython
