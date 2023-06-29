
## Calling insertions and deletions from a multiple alignment format (MAF) file 

There are limited tools to detect insertions and deletions directly from a multiple alignment, rather than from aligned reads. After testing tools such as [smartsv2](https://github.com/EichlerLab/smrtsv2) and [syri](https://schneebergerlab.github.io/syri/) for possible adaption to this task, I decided to take a simpler approach that I hoped would be better at detecting small indels <~1000bp from highly accurate Cactus multiple alignments. The basic idea is to use a single genome as a reference coordinate system and then to find indels that occur in a specific group of species (bats) relative to this reference (human). I assumed that any bat-specific indel occuring in multiple bats that is supported by the absence of this indel in human, mouse and pig (abbreviated as hg, mm and sscrofa) is likely to be a true lineage-specific indel. I did model potential impacts of alignment error, and I necessarily split any indels that overlap multiple alignment blocks as indels are called indepently for each block.

**Warning**: This custom script for calling lineage-specific small indels is not thoroughly benchmarked or tested. It is provided primarily for reproducibility of a custom analysis. I manually validated all indels discussed in the manuscript. Any users are cautioned to validate their results.  

To call deletions from a MAF file, we can use the following command on the provided example file that uses the human genome as the reference. 
```
# General usage relies on six fixed-order user arguments
python call_indels_from_maf.py <alignment.maf> <coordinate_species> <alignment_reference> <ingroup species list>  <outgroup species list> <ins|del> > output
# For example:
python call_indels_from_maf.py data/human_example.maf hg hg phydis,meso,jam,rhifer,drotund,mymy  hg,mm,sscrofa del > deletions.txt
```

Similarly, to call insertions, switch the in and outgroups and use a MAF file with a bat (*Myotis*, abbreviated as mymy). Note that indels are always called by identifying alignment gaps in the ingroup species, so to find bat-specific insertions the ingroups is `hg,mm,sscrofa`.

```
python call_indels_from_maf.py bat_example.maf hg mymy hg,mm,sscrofa phydis,meso,jam,rhifer,drotund,mymy ins > insertions.txt
```

The output files have thirteen tab-separated columns, which are explained below.

```
1) Reference coordinate species on which the indel coordinates are based (not necessarily the reference species in the MAF)
2) Chromosome
3) Indel start
4) Indel end
5) Indel length (bp)
6) Indel type (insertion or deletion)
7) Number of aligned ingroup species
8) Number of aligned outgroup species
9) Number of ingroup species with an overlapping indel
10) Number of outgroup species with an overlapping indel
11) Number of ingroup species with an overlapping indel with matching start and end positions
12) Number of outgroup species with an overlapping indel with matching start and end positions
13) Comma-separated list of species with the indel
```

One of the challenges with this method of indel detection is that a shared lineage-specific indel at a genomic position may not be perfectly identical. This can occur for both biological as well as technical reasons including convergent evolution, misalignment, or misassembly. We therefore tried to include sufficient information to either include or exclude overlapping indels that do not have an exact match in the start and end positions.

Using the raw output of the script, we can now filter cases where any of the outgroups (hg,mm,sscrofa) have any overlapping deletion and we can require that >5 ingroups and >2 outgroups are present. This will strictly identify only the most confident lineage-specific indels with exact matches.

```
awk '$7>4 && $8>2 && $11 >5 && $12 == 0 && $10==0' deletions.txt > deletions_filt.txt
```

Similarly, for the indels we can filter them like so.

```
awk '$7>2 && $8>5 && $11>2 && $12==0 && $10==0' insertions.txt > insertions_filt.txt
```

A script can then be used to extract alignment blocks containing the indels for validation or further analysis. The indel can be converted to bed by extracting columsn 2,3, and 4. Note that this script will write one file per row in the bed file that can be found in the MAF input file.
```
cut -f2-4 deletions_filt.txt > deletions.bed
python get_indel_fasta.py alignment.maf deletions.bed
```

To extract the indel-containing alignment block as a MAF file, the binary [mafsInRegion](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) can be used.
