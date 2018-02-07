library(readr)

allgenes <- read_tsv('~/Dropbox/BTB/resources/Ensembl/genes.txt')
tumorgenes <- read_csv('~/Dropbox/BTB/resources/COSMIC/cancer_gene_census.csv')
hotspots_inframe <- read_delim('~/Dropbox/BTB/resources/MSK hotspots V2 (25000 samples)/hotspots_v2_inframe.csv',delim = ';')
hotspots_snv <- read_delim('~/Dropbox/BTB/resources/MSK hotspots V2 (25000 samples)/hotspots_v2_snv.csv',delim=';')



## Prepare hotspot mutation definitions



## Small variant filtering
## Apply to both T and N small variants
## 2+ or 5+ obervations in SweGen

## Allocate small variants to tier 1, 2 and 3
## Rules: Tier 1: known hotspots or truncating in cancer census gene.
##        Tier 2: other mutation with effect in cancer census gene, not "benign" in ClinVar
##        Tier 3: other mutation with effect
##        All: Annotate with somatic LOH

## Structural variant filtering
## Apply to both T and N SVs
## 2+ or 5+ observations in SweGen

## Allocate structural variants to tiers
## Rules: Tier 1: Small "Deletion" or "Duplication" includes a cancer census gene. or either end of a SV is within or adjacent to cancer census gene
##        Tier 2: "Imprecise" variant with the above
##        Tier 3: Either end is near a coding gene, or small "Deletion" or "Amplification" with genes inside
##        Small/Local: Annotate with somatic LOH

## Copy number variant filtering
## Somatic amplification: remove too large
## Somatic homozygous deletion: remove too large
## Germline deletion: remove if known CNV

## Allocate copy number variants to tiers
## Rules: Tier 1: Amplification or Homozygous deletion implied by both ASCAT and Freec and contains a cancer census gene. Germline deletion of cancer census gene
##        Tier 2: Amplification or Homozygous deletion implied by either ASCAT or Freec and contains a cancer census gene. Germline duplication of cancer census gene
##        Tier 3: Amplification or Homozygous deletion implied by either ASCAT or Freec and contains any coding gene. Germline (non cnv) deletion of any coding gene






