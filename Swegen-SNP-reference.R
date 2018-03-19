




# counts=read.table('counts.frq.count',header = F,skip=1,sep='\t',fill=T,col.names = paste0('V',1:11),stringsAsFactors = F)
# 
# head(counts)
# dim(counts)
# colnames(counts)=c('chr','pos','alleles','total_fq',paste0('A',seq_len(7)))
# unique(counts$chr)
# counts$chr=substr(counts$chr,4,6)


# file=readLines('counts.frq.count')
# file=file[-1]
# file=strsplit(file,'\t')
# save(file,file='counts.frq.count.Rdata')
# load('counts.frq.count.Rdata')


## this was used @ Bianca
load('counts.frq.count.Rdata')
names=rep('',150e6)
values=rep(0,150e6)
i=1
p0=proc.time()[3]
n=0
for ( entry in file ) {
  n=n+1
  if (n %% 10000 ==0) {
    p1=proc.time()[3]
    dp=p1-p0
    cat(date(),'average',round(10000/dp), 'per second. Progress is', n,'\n')
    p0=p1
  }
  for (allele in entry[-(1:4)]) {
    temp=strsplit(allele,':')[[1]]
    names[i]=paste(substr(entry[1],start = 4,stop = 6),entry[2],temp[1])
    values[i]=as.numeric(temp[2])
    i=i+1
  }
}
ix <- names!=''
save(names,values,file='names_values.Rdata',compress = F)
values=values[ix]
names(values)=names[ix]
save(values,file='named_values.Rdata',compress = F)
## fanns inte på den R-versionen: hoppat till förra
library(hash)
snp_hash=hash(keys=names(values),values=values)
save(snp_hash,file='snp_hash.Rdata',compress = F)

# parallell vector version
load('counts.frq.count.Rdata')
n=length(file)
starts=seq(1,n,by=round(n/16))
ends=c(starts[-1]-1,n)
library(doParallel)
registerDoParallel(cores=16)
results <- foreach(i = 1:16,.combine = c) %dopar% {
  names=rep('',50e6)
  values=rep(0,50e6)
  i=1
  for (entry in file[starts[i]:ends[i]]) { #
    for (allele in entry[-(1:4)]) {
      temp=strsplit(allele,':')[[1]]
      names[i]=paste(substr(entry[1],start = 4,stop = 6),entry[2],temp[1])
      values[i]=as.numeric(temp[2])
      i=i+1
    }
  }
  names(values)=names
  values[values>0]
}
save(results,file='snp_results.Rdata',compress = F)
library(hash)
snp_hash=hash(keys=names(results),values=results)
save(snp_hash,file='snp_hash.Rdata',compress = F)




install.packages('hash')
library(hash)
snp_hash=hash(keys=names,values=values)
save(snp_hash,file='snp_hash.Rdata',compress = F)


# parallel hash version...
snps=hash()
registerDoParallel(cores=2)
foreach(entry = file[1:10],.export = 'snps') %dopar% {
  for (allele in entry[-(1:4)]) {
    temp=strsplit(allele,':')[[1]]
    snps[[paste(substr(entry[1],start = 4,stop = 6),entry[2],temp[1])]] <- as.numeric(temp[2])
  }
}


## put in filehash
load('~/Data/BTBdata/resources/SWEGEN snps/named_values.Rdata')
names=names(values)
names(values)=NULL
library(filehash)
dbCreate('myDB')
db <- dbInit('myDB')
p0=proc.time()[3]
for (i in 1:1000) {
  db[[names[i]]] <- values[i]
  if (i %% 100==0) { 
    p1=proc.time()[3]
    dp=p1-p0
    cat(date(),'average',round(100/dp), 'per second. Progress is', i,'\n')
    p0=p1
  }
} # 1kps...

indb=dbList(db)
dbMultiFetch(indb)


## in data table
library(data.table)
swegen_hg38_allele_counts=data.table(name=names,value=values)
setkey(swegen_hg38_allele_counts,name)
tables()
save(swegen_hg38_allele_counts,file='../resources/SWEGEN snps/swegen_hg38_allele_counts.Rdata',compress = F)





### New code for reading individual VCFs
library(VariantAnnotation)
library(data.table)
setwd('~')
files=dir(path = 'swegenVCF',pattern = 'vcf.gz$',full.names = T)

vcf=readGeno(files[1],x='GT')
snptable=data.table(name=unique(rownames(vcf)),value=1)
setkey(snptable,name)

for (i in 2:length(files)) {
  cat(i,files[i])
  vcf=readGeno(files[i],x='GT')
  cat('...read VCF.\n')
  r=unique(rownames(vcf))
  
  m = match(r, snptable[,name])
  
  snptable[m,'value']=snptable[m,'value']+1
  
  new=data.table(name=r[is.na(m)],value=1)
  snptable=rbind(snptable,new)
}


save(snptable,file='snptable_swegen_fromVCF.Rdata',compress = F)

setkey(snptable,name)
fwrite(snptable,file = '../resources/SWEGEN snps/snptable_swegen_fromVCF.csv')





