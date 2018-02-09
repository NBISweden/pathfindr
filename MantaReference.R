

## This script was run using all the diploid vcfs available within each project @ Bianca

library(VariantAnnotation)

allkeys=NULL
passkeys=NULL

files="~/Data/BTBdata/VEP/Manta_P2233_103T_vs_P2233_122N.diploidSV.ann.vcf" # dummy
files=dir(pattern='diploidSV.vcf.vep.ann.vcf') # if annotated
files=dir(pattern='diploidSV.vcf.gz') # if not annotated

for (file in files) {
  
  cat('now',file,'\n')
  
  nstruc_vcf=readVcf(file = file,genome = 'GRCh38')
  g=geno(nstruc_vcf)
  inf=info(nstruc_vcf)
  rr=rowRanges(nstruc_vcf)
  
  # manipulate into a data frame with relevant data
  nstruc=data.frame(rr,stringsAsFactors = F)
  nstruc$chr <- substr(x = nstruc$seqnames,start = 4,stop = 6)
  nstruc$altchr=nstruc$chr
  nstruc$altpos=nstruc$end
  # loop through all and extract endpoint chr and pos
  for (i in 1:nrow(nstruc)) {
    t=strsplit(x = nstruc$ALT[[i]],split = ':')[[1]]
    if (length(t)>1 & t[1]!="<DUP") {
      nstruc$altchr[i]=strsplit(x = t[1],split = 'chr')[[1]][2]
      tt=strsplit(t[2],'\\[')[[1]][1]; tt=strsplit(tt[1],'\\]')[[1]][1]
      nstruc$altpos[i]=as.numeric(tt)
    }
  }
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
  # add unlque keys from this pat to the lists
  passkeys=c(passkeys,unique(key[nstruc$FILTER=='PASS']))
  allkeys=c(allkeys,unique(key))
}

#done looping over VCFs
save(passkeys,allkeys,file='manta_keys_2016xxx.Rdata')
#save(passkeys,allkeys,file='manta_keys_2017xxx.Rdata')

## Local processing of keysets:
setwd('~/Dropbox/BTB/resources/Manta ref files/')
load('manta_keys_2016xxx.Rdata')
all_2016=allkeys
pass_2016=passkeys
load('manta_keys_2017xxx.Rdata')
all_2017=allkeys
pass_2017=passkeys
# make tables
manta_reference_all=data.frame(sort(table(c(all_2016,all_2017)),decreasing = T),stringsAsFactors = F)
manta_reference_pass=data.frame(sort(table(c(pass_2016,pass_2017)),decreasing = T),stringsAsFactors = F)
colnames(manta_reference_all) <- colnames(manta_reference_pass) <- c('name','value')
btb_hg38_manta_all <- data.table(manta_reference_all)
setkey(swegen_hg38_manta_all,name)
btb_hg38_manta_pass <- data.table(manta_reference_pass)
setkey(btb_hg38_manta_pass,name)
save(btb_hg38_manta_all,btb_hg38_manta_pass,file='BTB_hg38_manta_counts.Rdata')


