library("genoPlotR")

main <- "Type I IFN locus"
hgseg <- read_dna_seg_from_tab("data/ifn/hg_NC_000009_12.seg.txt")
mmseg <- read_dna_seg_from_tab("data/ifn/mm_NC_000070_7.seg.txt")
ecabseg <- read_dna_seg_from_tab("data/ifn/ecab.seg.txt")
ssseg <- read_dna_seg_from_tab("data/ifn/ss_NC_010443_5.seg.txt")
cfseg  <- read_dna_seg_from_tab("data/ifn/dog.seg.txt")
jamseg <- read_dna_seg_from_tab("data/ifn/jam.seg.cor.txt")
physeg <- read_dna_seg_from_tab("data/ifn/GCF_004126475.2_mPhyDis1.pri.v3_genomic.seg.txt")
myseg <-read_dna_seg_from_tab("data/ifn/GCF_014108235.1_mMyoMyo1.p_genomic.seg.txt")
molseg <-read_dna_seg_from_tab("data/ifn/GCF_014108415.1_mMolMol1.p_genomic.seg.txt")
rouseg <-read_dna_seg_from_tab("data/ifn/GCF_014176215.1_mRouAeg1.p_genomic.seg.txt")
palseg <-read_dna_seg_from_tab("data/ifn/palecto.seg.txt")
pipseg <-read_dna_seg_from_tab("data/ifn/GCF_014108245.1_mPipKuh1.p_genomic.seg.txt")
drotseg <- read_dna_seg_from_tab("data/ifn/drot_NW_020091304.seg.txt")
mesoseg <- read_dna_seg_from_tab("data/ifn/meso_congtig_73.seg.txt")
rhiseg <- read_dna_seg_from_tab("data/ifn/rhifer.seg.txt")

mydna_segs <- list(hgseg,mmseg,ssseg,ecabseg,cfseg,palseg,rouseg,rhiseg,molseg,pipseg,myseg,mesoseg,drotseg,physeg,jamseg)
names(mydna_segs) <- c("hg", "mm","sscrofa","ecab","canfam","palec","rou","rhi","mol","pip","mymy","meso","drot","phy","jam")
#my_blast <- list(hg_ss,
#                 ss_ecab,
#                 ecab_jam)
xlims <- list(c(21077869,21481068), # hg
              c(88440463,88798416), # mm
              c(201219848,201851495), # sscrofa
              c(39688443,40066564), # ecab
              c(41610249,41823034), # canfam
              c(1,247889), # palec
              c(61766381,62156918), # rousette
              c(18761480,18996012), # rhifer
              c(71936042,73024397), # mol
              c(47653629,47803515), # pip
              c(60629328,60914371),# my
              c(35119553,35283727),# meso
              c(14517538,14715153),# drot
              c(147238348,147543661), #phydis
              c(57573939,57830248) # jam
) 

my.phylog <- newick2phylog("((mm:0.89074,hg:0.89074):0.0421,(((canfam:0.710929,ecab:0.710929):0.025277,sscrofa:0.736205):0.015982,((rhi:0.575772,(rou:0.178771,palec:0.178771):0.397001):0.054143,((mol:0.468789,(pip:0.224725,mymy:0.224725):0.244065):0.067426,(meso:0.373131,(drot:0.252477,(phy:0.19779,jam:0.19779):0.054686):0.120655):0.163084):0.093699):0.122272):0.180653);")

pdf(file="ifn_synteny_plot.pdf",width=10,height=10)
plot_gene_map(mydna_segs,
              dna_seg_scale=TRUE,
              xlims = xlims,
              main=main,
              scale=FALSE,
              dna_seg_labels=NULL,
              tree=my.phylog)
dev.off()
# Now several genomes have a fragmented IFN locus, so we need to replot with the other scaffold
drot1 <- read_dna_seg_from_tab("data/ifn/drot_NW_020076331.seg.txt")
drot2 <- read_dna_seg_from_tab("data/ifn/drot_NW_020090274.seg.txt")
meso1 <- read_dna_seg_from_tab("data/ifn/meso_congtig_1174.seg.txt")
meso2 <- read_dna_seg_from_tab("data/ifn/meso_congtig_774.seg.txt")
mydna_segs <- list(meso1,meso2,drot1,drot2,molseg)
names(mydna_segs) <- c("meso_c1174","meso_c774","drot_NW020076331","drot_NW020090274","scaler")
my.phylog <- newick2phylog("((meso_c1174:0.89074,meso_c774:0.89074):0.0421,(drot_NW020076331:0.710929,drot_NW020090274:0.710929):0.025277,scaler:0.01):0.01;")

pdf(file="ifn_synteny_plot_fragments.pdf",width=10,height=10)
plot_gene_map(mydna_segs,
              dna_seg_scale=TRUE,
              main=main,
              scale=FALSE,
              dna_seg_labels=NULL,
              tree=my.phylog)
dev.off()
## IFITM

main <- "IFITM locus"
hgseg <- read_dna_seg_from_tab("data/ifitm/hg_NC_000011_10.seg.txt")
mmseg <- read_dna_seg_from_tab("data/ifitm/mm_NC_000073_7.seg.txt")
ecabseg <- read_dna_seg_from_tab("data/ifitm/ecab_NC_009155_3.seg.txt")
ssseg <- read_dna_seg_from_tab("data/ifitm/ss_NC_010444_4.seg.txt")
cfseg  <- read_dna_seg_from_tab("data/ifitm/dog.seg.txt")
jamseg <- read_dna_seg_from_tab("data/ifitm/jam.seg.txt")
physeg <- read_dna_seg_from_tab("data/ifitm/phy_NC_040908_2.seg.txt")
myseg <-read_dna_seg_from_tab("data/ifitm/myo_NW_023416332_seg.txt")
molseg <-read_dna_seg_from_tab("data/ifitm/mol_NW_023425368_1.seg.txt")
rouseg <-read_dna_seg_from_tab("data/ifitm/rou_NW_023416308_1.seg.txt")
pipseg <-read_dna_seg_from_tab("data/ifitm/pipkuh_NW_023425595.seg.txt")
drotseg <- read_dna_seg_from_tab("data/ifitm/drot.seg.txt")
mesoseg <- read_dna_seg_from_tab("data/ifitm/meso_contig_349.seg.txt")
rhiseg <- read_dna_seg_from_tab("data/ifitm/rhifer_JACAGC010000011.seg.txt")

mydna_segs <- list(hgseg,mmseg,ssseg,ecabseg,cfseg,rouseg,rhiseg,molseg,pipseg,myseg,mesoseg,drotseg,physeg,jamseg)
names(mydna_segs) <- c("hg", "mm","sscrofa","ecab","canfam","rou","rhi","mol","pip","mymy","meso","drot","phy","jam")

xlims <- list(c(298200,1750595), # hg
              c(140528871,141927526), # mm
              c(107035,1183309), # sscrofa
              c(34106806,36878811), # ecab
              c(25978345,26079310), # canfam # IFITM10 end: 46689944
              c(120834877,121733267), # rousette
              c(89008084,90080203), # rhifer
              c(2245859,2309830), # mol
              c(263967,283269), # pip
              c(2915712,2989554),# my
              c(304576,353595),# meso
              c(9361270,8354539),# drot
              c(170608691,171598164), #phydis
              c(854397,1728051) # jam
)
# Exclude IFITM10
xlims2 <- list(c(298200,320860), # hg
              c(140528871,140596805), # mm
              c(107035,157650), # sscrofa
              c(36814927,36878811), # ecab
              c(25978345,26079310), # canfam # IFITM10 end: 46689944
              c(121690472,121733267), # rousette
              c(90041735,90080203), # rhifer
              c(2245859,2309830), # mol
              c(263967,283269), # pip
              c(2915712,2989554),# my
              c(304576,353595),# meso
              c(9290545,9362049),# drot
              c(171432975,171598164), #phydis
              c(1641509,1728051) # jam
) 
my.phylog <- newick2phylog("((mm:0.89074,hg:0.89074):0.0421,(((canfam:0.710929,ecab:0.710929):0.025277,sscrofa:0.736205):0.015982,((rhi:0.575772,rou:0.178771):0.054143,((mol:0.468789,(pip:0.224725,mymy:0.224725):0.244065):0.067426,(meso:0.373131,(drot:0.252477,(phy:0.19779,jam:0.19779):0.054686):0.120655):0.163084):0.093699):0.122272):0.180653);")

pdf(file="ifitm_synteny_plot_fragments.pdf",width=10,height=10)
plot_gene_map(mydna_segs,
              dna_seg_scale=TRUE,
              xlims = xlims2,
              main=main,
              scale=FALSE,
              dna_seg_labels=NULL,
              tree=my.phylog)
dev.off()
# Now plot the extra contigs

meso1 <- read_dna_seg_from_tab("data/ifitm/meso_contig_214.seg.txt")
mm1 <- read_dna_seg_from_tab("data/ifitm/mm_NC_000082_7.seg.txt")
mol1 <- read_dna_seg_from_tab("data/ifitm/mol_NW_023425348.seg.txt")
pip1 <- read_dna_seg_from_tab("data/ifitm/pipkuh_NW_023425441.seg.txt")
pip2 <- read_dna_seg_from_tab("data/ifitm/pipkuh_NW_023425597.seg.txt")
my1 <- read_dna_seg_from_tab("data/ifitm/myo_NW_023416313_seg.txt")
my2 <- read_dna_seg_from_tab("data/ifitm/myo_NW_023416324_seg.txt")
my3 <- read_dna_seg_from_tab("data/ifitm/myo_NW_023416370_seg.txt")

mydna_segs <- list(hgseg,meso1,mm1,mol1,pip1,pip2,my1,my2,my3)
names(mydna_segs) <- c("hg", "meso_contig_214","mm_NC_000082_7","mol_NW_023425348","pipkuh_NW_023425441","pipkuh_NW_023425597","myo_NW_023416313","myo_NW_023416324","myo_NW_023416370")

my.phylog <- newick2phylog("((hg:0.89074):0.0421,(((meso_contig_214:0.710929,mm_NC_000082_7:0.710929):0.025277,mol_NW_023425348:0.736205):0.015982,((pipkuh_NW_023425441:0.575772,pipkuh_NW_023425597:0.178771):0.054143,((myo_NW_023416313:0.468789,(myo_NW_023416324:0.224725,myo_NW_023416370:0.224725):0.244065):0.067426):0.163084):0.122272):0.180653);")

pdf(file="ifitm_synteny_plot.pdf",width=10,height=10)
plot_gene_map(mydna_segs,
              dna_seg_scale=TRUE,
              main=main,
              scale=FALSE,
              dna_seg_labels=NULL,
              tree=my.phylog)

dev.off()

# Differentiating IFN alpha and omega
# Omega is ~5aa longer and has different STOP position
# Omega has a few distinct aa substitutions shared by multiple species near end
# Sequence similarity to related bat alpha/omega protein seqs

