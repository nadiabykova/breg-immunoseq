library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(stringdist)

multi_x = function(x,t,all_seqs,df){
  res.x = one_x(x,t,all_seqs,df)
  new.x = setdiff(res.x,x)
  sum.x = union(x,res.x)
  iter = 0
  while(length(new.x)>0){
    print(paste(iter," new.x:",length(new.x),sep=""))
    res.x = lapply(new.x,one_x,t=t,all_seqs=all_seqs,df=df) %>% unlist %>% unique
    new.x = setdiff(res.x,sum.x)
    sum.x = union(sum.x,res.x)
  }
  return(sum.x)
}
one_x = function(x,t,all_seqs,df){
  res.seqs = c()
  res.x = stringdist(x,all_seqs,method="hamming")
  seq_ids = which(res.x<=t)
  in_data = unique((df %>% filter(cdr13==x))$nMutations.total)
  for (xi in seq_ids){
    xs = all_seqs[xi]
    in_data.res = unique((df %>% filter(cdr13==xs))$nMutations.total)
    for (a in in_data){
      for (b in in_data.res){
        if (tnv(a,b)) res.seqs = c(res.seqs,xs)
      }
    }
  }
  return(unique(res.seqs))
}
tnv = function(a,b){
  a1=strsplit(a,";")[[1]]
  b1=strsplit(b,";")[[1]]
  res.l=0
  suma=0
  for (i in seq_along(a1)){
    a2=strsplit(a1[i],",")[[1]]
    b2=strsplit(b1[i],",")[[1]]
    l=intersect(a2,b2)
    res.l=res.l+length(l)
    suma=suma+length(a2)
  }
  if (res.l==0&(suma)>0) return(FALSE)
  else return(TRUE)
}

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
antibodies.file = args[2]

#antibodies.file = "antibodies/antibodies.aln.txt"

a = fread(antibodies.file)
a = a %>% filter(cloneId>=0) %>% mutate(chain = substr(bestVGene,1,3),
                                        cdr13 = paste(nSeqCDR1,nSeqFR2,nSeqCDR2,nSeqFR3,nSeqCDR3,sep=""),
                                        nMutations.total = paste(nMutationsInCDR1RelativeVRegion,nMutationsInFR2RelativeVRegion,
                                                                 nMutationsInCDR2RelativeVRegion,nMutationsInFR3RelativeVRegion,sep=";")) %>%
  select(chain,bestVGene,bestJGene,aaSeqCDR3,cdr13,nMutations.total)

a$mut.n = sapply(a$nMutations.total,function(x){
  grex = gregexpr("[A-Z0-9]+",x)[[1]]
  return(ifelse(grex[1]==-1,0,length(grex)))
})

## load data
files = list.files(dir,pattern = ".seed1.cdr13.txt")
res.df = lapply(files,function(file.name,dir){
  name = gsub('\\.seed1\\.cdr13\\.txt',"",file.name)
  print(name)
  chain = gsub('.*\\.',"",name)
  f = fread(paste(dir,file.name,sep=""))
  f$name=name
  f$chain=chain
  return(f)
},dir=dir) %>% rbindlist

## find bad clones
t = 5
bad_clones = c()
for (z in c(1:nrow(a))){
  x=a[z,"cdr13"]
  ch=a[z,"chain"]
  vgene=a[z,"bestVGene"]
  print(a[z,"aaSeqCDR3"])
  all_seqs = unique((res.df %>% filter(chain==ch,bestVGene==vgene))$cdr13)
  #res = one_x(x,t,all_seqs,res.df)
  res = multi_x(x,t,all_seqs,res.df)
  print(length(res))
  bad = res.df %>% filter(cdr13 %in% res)
  print(nrow(bad))
  bc = unique(bad$nSeqCDR3)
  print(length(bc))
  bad_clones = c(bad_clones,bc)
}
print(length(bad_clones))
print(res.df %>% filter(nSeqCDR3 %in% bad_clones))

# filter
files = list.files(dir,pattern = "seed1.clones.withCdr13.txt")
lapply(files,function(file.name,dir){
  print(file.name)
  name = gsub('\\.txt',"",file.name)
  print(name)
  f = fread(paste(dir,file.name,sep=""))
  #print(nrow(f))
  #print(sum(f$cloneCount))
  f.filtered = f %>% filter(!(nSeqCDR3 %in% bad_clones))
  print(nrow(f.filtered)/nrow(f))
  print(sum(f.filtered$cloneCount)/sum(f$cloneCount))
  out.clones.file = paste(dir,name,".decont.txt",sep="")
  write.table(f.filtered,out.clones.file,quote = F,sep="\t",row.names=F,na = "")
},dir=dir)
