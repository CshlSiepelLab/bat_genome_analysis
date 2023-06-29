# Setting up a custom genome browser resource for de novo genome assemblies

Using UCSC Track Hubs, it is straightforward to set up a custom genome browser with a diverse set of annotation tracks. Here we provide the configuration files for the *Artibeus jamaicensis* custom genome browser set up on our labs's bat genome browser [landing page](http://compgen.cshl.edu/bat/). We used the following data to set up the browser, converting it to the appropriate formats as shown in the section below.

* Genome assembly in FASTA format
* Gene annotation in GFF3 format
* Repeat annotation in GFF3 format
* Interspecies multiple alignment in MAF format
* Mapped RNAseq in BAM format

The configuration files needed are shown below.

* `Ajamaicensis_hub.txt` : The hub file is the main control file, which refers to a genomes file, and provides some basic metadata.
* `Ajamaicensis_genomes.txt` : The genomes file is another control file with the genome name and the file paths to the sequence and trackDb.
* `data/Ajamaicensis_trackDb.txt`: The trackDb contains the configuration of all of the tracks in the browser, including the file paths to the underlying data such as the gene annotation.

When the configuration files are set up, the browser can be accessed simply by adding the location of the hub file to a UCSC URL, for example: [*A. jamaicensis* genome browser](http://genome.ucsc.edu/cgi-bin/hgHubConnect?hgHub_do_redirect=on&hgHubConnect.remakeTrackHub=on&hgHub_do_firstDb=1&hubUrl=http://compgen.cshl.edu/bat/Ajamaicensis_hub.txt).

## Converting standard genomic formats to UCSC track-compatible formats

The genome assembly is converted to 2bit format.

```
faToTwoBit ajamaicensis.fa ajamaicensis.2bit
```

A standard MAF alignment file was converted to binary big bed format.
```
mafToBigMaf ajamaicensis ajamaicensis.maf ajamaicensis.bed
bedToBigBed -tab -type=bed3+1 -as=bigMaf.as ajamaicensis.bed ajamaicensis_chromsizes.txt ajamaicensis.bb
```
The ajamaicensis_chromsizes.txt file is a two-column tab-separarted file with scaffolds and scaffold sizes in bp, for instance:
```
head -4 ajamaicensis_chromsizes.txt
scaffold_1964   73234996
scaffold_916    68176794
scaffold_474    50465402
scaffold_537    46405026
```

All gene annotation GFF3 files were converted to big bed format in three steps.
```
gff3ToGenePred ajamaicensis.gff3 ajamaicensis.gp
genePredToBed  ajamaicensis.gp ajamaicensis.bed
bedToBigBed ajamaicensis.gp ajamaicensis.bb
```

## Software

All binaries used above can be downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).

## Resources 

https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Setup
