




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

