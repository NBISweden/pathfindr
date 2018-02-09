library(data.table)

allgenes <- fread('~/Data/BTBdata//resources/Ensembl/genes.txt')
allgenes <- allgenes[`Chromosome/scaffold name` %in% c(1:22,'X','Y')]
allgenes <- allgenes[order(`Chromosome/scaffold name`,`Transcript start (bp)`)]
# Cosmic genes have Ties, Role in Cancer, Fusion, Germline/Somatic and much more....!!!!
tumorgenes <- fread('~/Data//BTBdata//resources/COSMIC/cancer_gene_census.csv')
local_tumorgenes <- fread('~/Data//BTBdata//resources/Teresita gene list/Brain_tumor_list_of_genes_2017.txt',header = F)
hotspots_inframe <- fread('~/Data//BTBdata//resources/MSK hotspots V2 (25000 samples)/hotspots_v2_inframe.csv')
hotspots_snv <- fread('~/Data//BTBdata//resources/MSK hotspots V2 (25000 samples)/hotspots_v2_snv.csv')
# to add: hotspots not in MSK-impact?? Are there any?
oncokb_all <- fread('~/Data/BTBdata/resources/OncoKB/OncoKB-allAnnotatedVariants.txt')
oncokb_act <- fread('~/Data/BTBdata/resources/OncoKB/OncoKB-allActionableVariants.txt')
## Prepare hotspot mutation definitions


## Small/structural variant filtering
## 2+ observations in SWEgen/other ref

## Copy number variant filtering
## Somatic amplification: remove too large?
## Somatic homozygous deletion: remove too large
## Germline deletion: remove if known CNV







