library(data.table)


chrsz <- data.table(
  chr = c("1", "2", "3", "4", "5", "6", "7", "8", 
          "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
          "20", "21", "22", "X", "Y"), 
  starts = c(0, 253902921, 500669562, 
             703689723, 898025604, 1084280925, 1259870526, 1424144247, 1574191848, 
             1717288929, 1855896810, 1995944331, 2134165212, 2253485013, 2363996694, 
             2470839615, 2565902256, 2654091777, 2739309138, 2802869979, 2871941700, 
             2923598421, 2979326382, 3140008743), 
  length = c(248902921, 241766641, 
             198020161, 189335881, 181255321, 170589601, 159273721, 145047601, 
             138097081, 133607881, 135047521, 133220881, 114319801, 105511681, 
             101842921, 90062641, 83189521, 80217361, 58560841, 64071721, 
             46656721, 50727961, 155682361, 56827081)
)
# from https://genome.ucsc.edu/cgi-bin/hgTable
allgenes <- fread('~/Data/BTBdata//resources/Ensembl/Canonical_Gencode.txt')
allgenes$chr <- substr(allgenes$`#hg38.knownCanonical.chrom`,4,6)
allgenes <- allgenes[chr %in% c(1:22,'X','Y'),-1]
colnames(allgenes)[1:3]=c('start','end','name')
allgenes$cumstart <- NA
allgenes$cumend <- NA
for (i in 1:nrow(chrsz)) {
  ix <- allgenes$chr==chrsz$chr[i]
  allgenes$cumstart[ix] <- allgenes$start[ix]+chrsz$starts[i]
  allgenes$cumend[ix] <- allgenes$end[ix]+chrsz$starts[i]
}
# Cosmic genes have Tiers, Role in Cancer, Fusion, Germline/Somatic and much more....!!!!
tumorgenes <- fread('~/Data//BTBdata//resources/COSMIC/cancer_gene_census.csv')
local_tumorgenes <- fread('~/Data//BTBdata//resources/Teresita gene list/Brain_tumor_list_of_genes_2017.txt',header = F)
hotspots_inframe <- fread('~/Data//BTBdata//resources/MSK hotspots V2 (25000 samples)/hotspots_v2_inframe.csv')
setkey(hotspots_inframe,Hugo_Symbol,Amino_Acid_Position)
hotspots_snv <- fread('~/Data//BTBdata//resources/MSK hotspots V2 (25000 samples)/hotspots_v2_snv.csv')
# to add: hotspots not in MSK-impact?? Are there any?
oncokb_all <- fread('~/Data/BTBdata/resources/OncoKB/OncoKB-allAnnotatedVariants.txt')
oncokb_act <- fread('~/Data/BTBdata/resources/OncoKB/OncoKB-allActionableVariants.txt')
# also to add: our panel of tumors and normals for recurrent FPs
## Save them as one file:
save(chrsz,allgenes,tumorgenes,local_tumorgenes,hotspots_inframe,hotspots_snv,oncokb_act,oncokb_all,file='~/reports/mutations_genes_ref.Rdata')

## Small/structural variant filtering
## 2+ observations in SWEgen/other ref

## Copy number variant filtering
## Somatic amplification: remove too large?
## Somatic homozygous deletion: remove too large
## Germline deletion: remove if known CNV







