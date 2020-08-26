library("ggplot2",quietly = TRUE, warn.conflicts = FALSE)
library("data.table",quietly = TRUE, warn.conflicts = FALSE)
library("dplyr",quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
sample.file = args[2]
out.file = args[3]

#dir = "fastq_final_results"
#sample.file = "sample_info_2.2_corr.txt"
#out.file = "metadata_corr.convert.txt"

files = list.files(dir, pattern="seed1.clones.withCdr13.decont.txt")
full = paste(dir,"/",files,sep="")
sample = sapply(files,function(x){gsub("\\.IG[HKL].*","",x)})

chain = sapply(files,function(x){
  if (grepl("IGH",x)) return("IGH")
  if (grepl("IGK",x)) return("IGK")
  if (grepl("IGL",x)) return("IGL")
  })
sample_chain = paste(sample,chain,sep="_")
df=data.frame(file.name=full,sample.id=sample_chain,sample=sample,chain=chain)

s = fread(sample.file)
m = merge(df,s,by=c("sample","chain"))
#m = m %>% filter(!(sample %in% c("lom-aug-3","lom-aug-4")))
m = m %>% select(file.name,sample.id,chain,cell_type,subject_id,status,status2,batch)
write.table(m,out.file,quote = F,sep="\t",row.names=F)
print("Making metadata finished!")
