# Plotting scripts for bat genome analyses

The following set of R scripts can be used to reproduce figures from the manuscript. The scripts refer to input files provided in the `data` directory. 

* `genoplot_IFN_long_short.R`: Plot synteny between the Type I Interferon region in a short-read and our long-read genome assembly of A. jamaicensis
* `plot_ifn_ifitm.R`: Plot a species tree showing the syntenic Type I Interferon and IFITM regions with their genes for each species
* `plot_repeat_landscape.R`: Plot a barchart to illustrate RepeatMasker outputs on repeat landscape and expansion in bats
* `prdm9_tree.R`: Plot a PRDM9 gene tree with a presence-absence heatmap for the three key PRDM9 domains

## Dependencies
The plotting scripts depend on four R packages.

* ggplot2
* dplyr
* genoPlotR
* phytools
