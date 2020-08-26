library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
fixer_path = args[2]

#dir = "fastq_final_results/"
#fixer_path = "scripts_final/aln_fixer.R"
source(fixer_path)
files = list.files(dir,pattern = ".seed1.clones.txt")

lapply(files,function(file.name,dir){
  name = gsub('\\.clones\\.txt',"",file.name)
  aln.file = paste(name,".aln.txt",sep="")
  print(file.name)
  print(aln.file)
  
  f = fread(paste(dir,file.name,sep=""))
  print(nrow(f))
 # aln.file = "test10.aln.txt"
  
  aln = fread(paste(dir,aln.file,sep=""))
  aln = aln %>% filter(cloneId>=0,!(nSeqCDR1==""|nSeqFR2==""|nSeqCDR2==""|nSeqFR3==""|nSeqCDR3==""))
  
  aln.fixed = fix.mutations(aln)
  
  aln.fixed = aln.fixed %>% mutate(cdr13 = paste(nSeqCDR1,nSeqFR2,nSeqCDR2,nSeqFR3,nSeqCDR3,sep=""),
                                   nMutations.total = paste(nMutationsInCDR1RelativeVRegion,nMutationsInFR2RelativeVRegion,
                                                            nMutationsInCDR2RelativeVRegion,nMutationsInFR3RelativeVRegion,sep=";"))
  aln.fixed = merge(aln.fixed,(f %>% select(cloneId,nSeqCDR3)),by=c("cloneId","nSeqCDR3"))

  aln.s = aln.fixed %>% group_by(cloneId,bestVGene,aaSeqCDR3,nSeqCDR3,cdr13,nMutations.total) %>% 
    summarize(n = n()) %>% arrange(cloneId,desc(n)) %>%
    group_by(cloneId) %>% slice(1) %>% mutate(mut.n = sapply(nMutations.total,function(x){
      grex = gregexpr("[A-Z0-9]+",x)[[1]]
      return(ifelse(grex[1]==-1,0,length(grex)))
      })
    )
  print(nrow(aln.s))
  
  f.filtered = f %>% filter(cloneId %in% aln.s$cloneId)
  if (nrow(f.filtered)!= nrow(aln.s)) print("clones n problem")
  test = merge(f,aln.s,by="cloneId")
  if (nrow(f.filtered)!= nrow(test)) print("test cloneId - nSeqCDR3 failed")
  
  print(nrow(f.filtered)/nrow(f))
  count.sum = sum(f.filtered$cloneCount)
  #print(count.sum)
  f.filtered = f.filtered %>% mutate(cloneFraction = cloneCount/count.sum)
  
  out.clones.file = paste(dir,name,".clones.withCdr13.txt",sep="")
  targets.file = paste(dir,name,".cdr13.txt",sep="")
  write.table(f.filtered,out.clones.file,quote = F,sep="\t",row.names=F,na = "")
  write.table(aln.s,targets.file,quote = F,sep="\t",row.names=F)
},dir=dir)
