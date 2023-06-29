library(phytools)

tree <- phytools::read.newick("data/RAxML_bestTree.prdm9v3_10bs_tidy_20062023")
domains <- read.delim("data/prdm9_pfam_simple_20062023.txt",sep="\t")
domains$KRAB <- as.numeric(domains$KRAB)
domains$SSXRD <- as.numeric(domains$SSXRD)
domains$SET <- as.numeric(domains$SET)
domains2 <- domains[,-1]
rownames(domains2) <- domains[,1]
dmat <- matrix(as.numeric(unlist(domains2)),nrow=nrow(domains2))
rownames(dmat) <- domains$Gene__Species

#par(mar=c(16,0,0,0)) ## reset margins to default
pdf(file="prdm9_tree_with_domains.pdf",width=10,height=10)
phylo.heatmap(tree,dmat,
              split=c(0.9,0.1),fsize=c(0.1,0.1,0.1),
              standardize=FALSE)
dev.off()

