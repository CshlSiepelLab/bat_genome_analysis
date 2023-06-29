library("genoPlotR")

# A simple BLAST search of the target region generates the key input data
#blastn -outfmt 6 -num_threads 6 -query seq.fa  -db db.fa >  out.blastn
# results are filtered to retain only matches with >90% identity

long_short <- read_comparison_from_blast("data/Ajamaicensis_IFN_scaffold_1964_GCF_014825515.1_IFN.blastn.switch_coords.pid90.txt")
main <- "Comparison of IFN in long and short read Ajam"
long <- read_dna_seg_from_tab("data/Ajamaicensis.seg.txt")
short <- read_dna_seg_from_tab("data/GCF_014825515.1_IFN.concat.seg.txt")
mydna_segs <- list(long,short)
my_blast <- list(long_short
)
xlims <- list(c(0,258000),
              c(0,535636)
)
mydna_segs <- list(short,long)
my_blast <- list(long_short)
xlims <- list(c(0,535636),
              c(0,258000)
)
pdf(file="IFN_short_vs_long_synteny_plot.pdf",width=10,height=10)
plot_gene_map(mydna_segs, my_blast,
              dna_seg_scale=TRUE,
              xlims = xlims,
              main=main,
              scale=FALSE)
dev.off()
