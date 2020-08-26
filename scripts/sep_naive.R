library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
dir2 = args[2]
inner_orig = args[3]
inner_naive = args[4]
inner_mem = args[5]
inner_naive0 = args[6]
inner_mem0 = args[7]

#dir = "fastq_final_vdj/"
#dir2 = "fastq_final_results/"
 # inner_orig = "orig/"
 # inner_orig_full = "orig_full/"
 # inner_naive = "naive3/"
 # inner_mem = "mem3/"
 # inner_naive0 = "naive0/"
 # inner_mem0 = "mem0/"
 

file.pattern = 'vdj\\.corr\\.func\\.csd\\.'
files = list.files(dir,pattern = "vdj.corr.func.csd.[Ll]")
dir.create(paste(dir,inner_orig,sep=""))
dir.create(paste(dir,inner_naive,sep=""))
dir.create(paste(dir,inner_mem,sep=""))
dir.create(paste(dir,inner_naive0,sep=""))
dir.create(paste(dir,inner_mem0,sep=""))

f1 = function(file.name){
  #f.orig = fread(paste(dir,"vdj.corr.func.csd.lom-01_S1_L001_IGK.txt",sep=""))
  #nrow(f.orig)
  f = fread(paste(dir,file.name,sep=""))
  prename = gsub(file.pattern,'',file.name)
  name = gsub('_.*','',prename)
  chain = gsub('.*_','',prename)
  chain = gsub('\\.txt','',chain)
  print(file.name)
  mut.file = list.files(dir2,pattern = paste(name,".*",chain,".seed1.cdr13.txt",sep=""))
  m = fread(paste(dir2,mut.file[1],sep=""))
  if (nrow(m[!is.na(m$mut.n),])!=nrow(m)) print("NAs in mutations file!")
  mf = merge(f,m,by.x="cdr3nt",by.y="nSeqCDR3")
  mf = mf %>% select(count,freq,cdr3nt,cdr3aa,v,d,j,VEnd,DStart,DEnd,JStart,mut.n)
  if (nrow(mf)!=nrow(f)) print("Hmm, lost ccdr13 clonotypes")
  #print(nrow(mf))
  if (nrow(mf)<10) {
    print("SMALL file!")
  }
  #naive3
  mf.naive3= mf %>% filter(mut.n<=3)
  mf.mem3 = mf %>% filter(mut.n>3)
  naive.sum3 = sum(mf.naive3$count)
  mem.sum3 = sum(mf.mem3$count)
  mf.naive3 = mf.naive3 %>% mutate(freq = count/naive.sum3) %>% select(-mut.n)
  mf.mem3 = mf.mem3 %>% mutate(freq = count/mem.sum3) %>% select(-mut.n)
  #naive0
  mf.naive0= mf %>% filter(mut.n<=0)
  mf.mem0 = mf %>% filter(mut.n>0)
  naive.sum0 = sum(mf.naive0$count)
  mem.sum0 = sum(mf.mem0$count)
  mf.naive0 = mf.naive0 %>% mutate(freq = count/naive.sum0) %>% select(-mut.n)
  mf.mem0 = mf.mem0 %>% mutate(freq = count/mem.sum0) %>% select(-mut.n)
  
  mf = mf %>% select(-mut.n)
  write.table(f,paste(dir,inner_orig,file.name,sep=""),
              sep="\t",row.names = F,quote = F)
  write.table(mf.naive3,paste(dir,inner_naive,file.name,sep=""),
              sep="\t",row.names = F,quote = F)
  write.table(mf.mem3,paste(dir,inner_mem,file.name,sep=""),
              sep="\t",row.names = F,quote = F)
  write.table(mf.naive0,paste(dir,inner_naive0,file.name,sep=""),
              sep="\t",row.names = F,quote = F)
  write.table(mf.mem0,paste(dir,inner_mem0,file.name,sep=""),
              sep="\t",row.names = F,quote = F)
  return(data.frame(name=name,chain=chain,reads.total = sum(mf$count),reads.naive3=naive.sum3,reads.mem3=mem.sum3,
                    reads.naive0=naive.sum0,reads.mem0=mem.sum0))
}
  
res = lapply(files,f1) %>% rbindlist
print("Separation of naive and memory finished!")
res.sum = res %>% group_by(chain) %>% summarize(inner_orig=min(reads.total),inner_naive0=min(reads.naive0),
                                                inner_naive=min(reads.naive3),
                                      inner_mem0=min(reads.mem0),inner_mem=min(reads.mem3))
print(res.sum)
m = melt(res.sum,id.vars = c("chain"), variable.factor=FALSE)
m$variable = as.character(m$variable)
for (i in c(1:nrow(m))){
  ch=m[i,"chain"]
  var = m[i,"variable"]
  #print(var)
  val = m[i,"value"]
  file = paste(dir,eval(as.name(var)),"min.count.",ch,".txt",sep="")
  #print(file)
  write(val,file,sep="\t")
}
# print(res %>% filter(chain=="IGH"))
# print(res %>% filter(chain=="IGK"))
# print(res %>% filter(chain=="IGL"))



