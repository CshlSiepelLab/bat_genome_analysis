library(topGO)
library(dplyr)
library(GO.db)

args = commandArgs(trailingOnly=TRUE)
# Significant genes
coreGenes <- read.table(args[2])
# Gene to GO mappings
geneID2GO <- readMappings(args[1])
geneNames <- names(geneID2GO)

myInterestingGenes <- as.character(coreGenes$V1)

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGO.weight <- runTest(GOdata, algorithm = "weight", statistic = "Fisher" )
resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )

allRes <- GenTable(GOdata, elimKS = resultTopGO.elim,
                   orderBy = "elimKS",
                   topNodes = 356)
myterms <- allRes$GO.ID
mygenes <- genesInTerm(GOdata, myterms)
# Remove low p-val 
allRes <- allRes[which(allRes$elimKS<=0.01),]
outname = paste(args[2],"_topGO_BP_elim_fisher001.txt",sep = "")
write.table(allRes, file = outname, sep = "\t", quote = F, col.names = T, row.names = F)

var=c()
for (i in 1:length(myterms))
{
myterm=myterms[i]
mygenesforterm= mygenes[myterm][[1]]
mygenesforterm=paste(mygenesforterm, collapse=',')
var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
outname2 = paste(args[2],"_elim_genetoGOmapping.txt",sep = "")
write.table(var,outname2,sep="\t",quote=F)

# Repeat above code for weight method
allRes <- GenTable(GOdata, elimKS = resultTopGO.weight,
                   orderBy = "elimKS",
                   topNodes = 356)
myterms <- allRes$GO.ID
mygenes <- genesInTerm(GOdata, myterms)
# Remove low p-val 
allRes <- allRes[which(allRes$elimKS<=0.01),]
outname = paste(args[2],"_topGO_BP_weight_fisher001.txt",sep = "")
write.table(allRes, file = outname, sep = "\t", quote = F, col.names = T, row.names = F)

var=c()
for (i in 1:length(myterms))
{
myterm=myterms[i]
mygenesforterm= mygenes[myterm][[1]]
mygenesforterm=paste(mygenesforterm, collapse=',')
var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
outname2 = paste(args[2],"_weight_genetoGOmapping.txt",sep = "")
write.table(var,outname2,sep="\t",quote=F)

