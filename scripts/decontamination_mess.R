library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(stringdist)

a = fread("antibodies/antibodies.aln.txt")
a = a %>% filter(cloneId>=0) %>% mutate(chain = substr(bestVGene,1,3),
                                        cdr13 = paste(nSeqCDR1,nSeqFR2,nSeqCDR2,nSeqFR3,nSeqCDR3,sep=""),
                                        nMutations.total = paste(nMutationsCDR1,nMutationsFR2,
                                                                 nMutationsCDR2,nMutationsFR3,sep=";")) %>%
  select(chain,bestVGene,aaSeqCDR3,cdr13,nMutations.total)

a$mut.n = sapply(a$nMutations.total,function(x){
  grex = gregexpr("[A-Z0-9]+",x)[[1]]
  return(ifelse(grex[1]==-1,0,length(grex)))
})

dir = "fastq_final_results/"
files = list.files(dir,pattern = ".seed1.aln.txt")
res.df = lapply(files,function(file.name,dir){
  name = gsub('\\.seed1\\.aln\\.txt',"",file.name)
  print(name)
  chain = gsub('.*\\.',"",name)
  f = fread(paste(dir,file.name,sep=""))
  f = f %>%
    filter(!(nSeqCDR1==""|nSeqFR2==""|nSeqCDR2==""|nSeqFR3==""|nSeqCDR3=="")) %>% 
    filter(cloneId>=0) %>%
    mutate(cdr13 = paste(nSeqCDR1,nSeqFR2,nSeqCDR2,nSeqFR3,nSeqCDR3,sep=""),
           nMutations.total = paste(nMutationsCDR1,nMutationsFR2,
                                    nMutationsCDR2,nMutationsFR3,sep=";"))
  nrow(f)
  f$mut.n = sapply(f$nMutations.total,function(x){
      grex = gregexpr("[A-Z0-9]+",x)[[1]]
      return(ifelse(grex[1]==-1,0,length(grex)))
  })
  fa = f %>% group_by(cloneId,bestVGene,aaSeqCDR3,cdr13,mut.n,nMutations.total) %>% summarize(n = n()) %>% arrange(cloneId,desc(n)) %>%
    group_by(cloneId) %>% slice(1)
  fa$name=name
  fa$chain=chain
  return(fa)
},dir=dir) %>% rbindlist

multi_x = function(x,t,all_seqs,df){
  sum.x = x
  res.x = one_x(x,t,all_seqs,df)
  new.x = setdiff(res.x,sum.x)
  sum.x = union(sum.x,res.x)
  while(length(new.x)>0&length(new.x)<100){
    print(length(new.x))
    print(length(sum.x))
    res.x = one_x(new.x,t,all_seqs,df)
    new.x = setdiff(res.x,sum.x)
    sum.x = union(sum.x,res.x)
    if (length(new.x)>=100) print("exit by 100")
  }
  return(sum.x)
}
one_x = function(x,t,all_seqs,df){
  res.seqs = c()
  for (j in seq_along(x)){
    xj = x[j]
    res.x = stringdist(xj,all_seqs,method="hamming")
    seq_ids = which(res.x<=t)
    in_data = unique((df %>% filter(cdr13==xj))$nMutations.total)
    if (length(in_data)>1) {
      print(paste("in_data ",in_data))
      print(xs)
    }
    for (xi in seq_ids){
      xs = all_seqs[xi]
      xd = res.x[xi]
      in_data.res = unique((df %>% filter(cdr13==xs))$nMutations.total)
      if (length(in_data.res)>1) {
        print(paste("in_data res",in_data.res))
        print(xs)
      }
      for (a in in_data){
        for (b in in_data.res){
          if (tnv(a,b)) {
            res.seqs = c(res.seqs,xs)
         #   print("add")
          }
        }
      }
    }
  }
  return(unique(res.seqs))
}
final_x = function(x,ch,vgene,df,maxt,t){
  print(ch)
  print(x)
  all_seqs = unique((df %>% filter(chain==ch,bestVGene==vgene))$cdr13)
  xs = stringdist(x,all_seqs,method="hamming")
  seqs_maxt = all_seqs[which(xs<=maxt)]
  res = multi_x(x,t,seqs_maxt,df)
  print(length(res))
  return(res)
}
tnv = function(a,b){
  a1=strsplit(a,";")[[1]]
  b1=strsplit(b,";")[[1]]
 # print(a1)
#  print(b1)
  res.l=0
  suma=0
  sumb=0
  for (i in seq_along(a1)){
    a2=strsplit(a1[i],",")[[1]]
    b2=strsplit(b1[i],",")[[1]]
    l=intersect(a2,b2)
    res.l=res.l+length(l)
    suma=suma+length(a2)
    sumb=sumb+length(b2)
  }
  if (res.l==0&(suma)>0) return(FALSE)
  else return(TRUE)
}
total.df = rbind(a %>% select(chain,bestVGene,aaSeqCDR3,cdr13,nMutations.total),
                 res.df %>% select(chain,bestVGene,aaSeqCDR3,cdr13,nMutations.total))
for (z in c(1:nrow(a))){
  x=a[z,"cdr13"]
  ch=a[z,"chain"]
  vgene=a[z,"bestVGene"]
  print(a[z,"aaSeqCDR3"])
  res = final_x(x,ch,vgene,total.df,20,5)
  print(length(res))
}
maxt = 20
for (z in c(1:nrow(a))){
  x=a[z,"cdr13"]
  ch=a[z,"chain"]
  vgene=a[z,"bestVGene"]
  print(a[z,"aaSeqCDR3"])
  all_seqs = unique((total.df %>% filter(chain==ch,bestVGene==vgene))$cdr13)
  xs = stringdist(x,all_seqs,method="hamming")
  seqs_maxt = all_seqs[which(xs<=maxt)]
  res2 = one_x(x,5,seqs_maxt,total.df)
  print(length(res))
}


a %>% filter(aaSeqCDR3=="CQSYDSSNHGVF")

res.df %>% filter(cdr13 %in% setdiff(res,res2))

st = res.df %>% group_by(bestVGene,cdr13,nMutations.total) %>% summarise(n1=n()) %>% 
  group_by(bestVGene,cdr13) %>% summarise(n2=n()) %>% filter(n2>=2)
for (s in st$cdr13){
  print(res.df %>% filter(cdr13==s))  
}


res.df %>% filter(chain=="IGL") %>% group_by(cdr13) %>% nrow
