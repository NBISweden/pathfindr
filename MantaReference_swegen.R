

## This script was run using all the diploid vcfs available within each project @ Bianca

library(VariantAnnotation)
library(data.table)



# For SWEGEN 
folders=dir(path = 'Manta_SAMPLES/',include.dirs = T)
files=paste0('Manta_SAMPLES/',folders,'/results/variants/diploidSV.vcf.gz') 

allkeys=NULL
passkeys=NULL
for (i in 376:length(files)) if (!i %in% c(372,373,375)) {
  file=files[i]
  cat(i,':->',file,'\n')
  
  nstruc_vcf=readVcf(file = file,genome = 'GRCh38')
  g=geno(nstruc_vcf)
  inf=info(nstruc_vcf)
  rr=rowRanges(nstruc_vcf)
  
  # manipulate into a data frame with relevant data
  nstruc=as.data.table(rr)
  nstruc$chr <- substr(x = nstruc$seqnames,start = 4,stop = 6)
  #nstruc=nstruc[chr %in% c(1:22,'X','Y')]
  
  ## Key has only chr,start,end
  key=nstruc[,c('chr','start','end')]
  key$imprecise='(pr)'
  ## If imprecise, round the pos to 10
  ix=inf$IMPRECISE==T
  key$imprecise[ix]='(impr)'
  key$start[ix]=round(key$start[ix]/10)*10
  key$end[ix]=round(key$end[ix]/10)*10
  key=paste(key$chr,key$start,key$end,key$imprecise)
  fname=basename(file)
  # add unique keys from this pat to the long vectors (to be tabled below)
  passkeys=c(passkeys,unique(key[nstruc$FILTER=='PASS']))
  allkeys=c(allkeys,unique(key))
}

save(passkeys,allkeys,file='~/mantakeys.Rdata')

## Processing of keysets:
# make tables
manta_reference_all=as.data.table(sort(table(allkeys),decreasing = T))
manta_reference_pass=as.data.table(sort(table(passkeys),decreasing = T))
colnames(manta_reference_all) <- colnames(manta_reference_pass) <- c('name','value')
swegen_hg38_manta_all <- manta_reference_all
setkey(swegen_hg38_manta_all,name)
swegen_hg38_manta_pass <- manta_reference_pass
setkey(swegen_hg38_manta_pass,name)
save(swegen_hg38_manta_all,swegen_hg38_manta_pass,file='~/Swegen_hg38_Manta_counts.Rdata')

fwrite(swegen_hg38_manta_all,'~/reports/swegen_sv_counts.csv')


