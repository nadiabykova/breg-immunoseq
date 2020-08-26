library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(edgeR,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
aln.file = args[1]
#aln.file = "fastq_2.2_seed_results/lom-01_S1_L001.IGH.cdr3.seed1.cdr3.aln.txt"
f = read.csv(aln.file,sep="\t",header=T,stringsAsFactors = F)
f = f[f$cloneId!=-1,]
f$cdr13 = paste(f$nSeqCDR1,f$nSeqFR2,f$nSeqCDR2,f$nSeqFR3,f$nSeqCDR3,sep="")
f$nMutations.total = paste(f$nMutationsCDR1,f$nMutationsFR2,f$nMutationsCDR2,f$nMutationsFR3,sep=";")
f$mut.n = sapply(f$nMutations.total,function(x){
  grex = gregexpr("[A-Z0-9]+",x)[[1]]
  return(ifelse(grex[1]==-1,0,length(grex)))
})
print(aln.file)
#print(nrow(f))
print(nrow(f[(f$nSeqCDR1==""|f$nSeqFR2==""|f$nSeqCDR2==""|f$nSeqFR3==""|f$nSeqCDR3==""),]))

f[(f$nSeqCDR1==""|f$nSeqFR2==""|f$nSeqCDR2==""|f$nSeqFR3==""|f$nSeqCDR3==""),"cdr13"] = ""
f[f$cdr13=="","mut.n"] = NA

f = f[,c("cloneId","nSeqCDR3","cdr13","mut.n")]
dt = data.table(f)

res = dt[,.N,by=list(cloneId,nSeqCDR3,cdr13,mut.n)]
res = res[order(cloneId,-N),]
res2 = res[,.SD[1],by="cloneId"]

print((res2 %>% filter(is.na(mut.n)) %>% nrow)/(res2 %>% nrow))

#outname = gsub("\\.aln\\.txt",".nMut.txt",aln.file)
#write.table(res2[,c("nSeqCDR3","mut.n")],outname,quote = F,sep="\t",row.names=F)
