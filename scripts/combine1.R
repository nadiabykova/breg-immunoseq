library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]

#dir = "IBH_25102019/orig/"
dirs = list.dirs(dir,full.names=F,recursive = F)

cols = c()
for (d in dirs){
  f = fread(paste(dir,"/",d,"/",d,"_all_tests.txt",sep=""))
  print(colnames(f))
  cols = union(cols,colnames(f))
}
total = data.frame()
data = data.frame()
dname = "level"
i = 1
while (dname %in% cols){
  dname = paste("level",i,sep="")
  i = i + 1
}
for (d in dirs){
  #total
  f = fread(paste(dir,"/",d,"/",d,"_all_tests.txt",sep=""))
  na.cols = setdiff(cols,colnames(f))
  if (length(na.cols)>0) f[,na.cols]=NA
  f[,dname] = d
  total = rbind(total,f)
  
  #data
  f2 = fread(paste(dir,"/",d,"/",d,"_data.txt",sep=""))
  f2[,dname] = d
  data = rbind(data,f2)
}
last_dir =  sub('/$','',dir)
last_dir = gsub('.*/','',last_dir)

lev_cols = c(colnames(total)[grepl("level",colnames(total))])

# data
data = data %>% select(lev_cols,everything())
write.table(data,paste(dir,last_dir,"_data.txt",sep=""),
            quote = F,sep="\t",row.names=F)

# total
total = total %>% select(lev_cols,test_name,test_type,subset,test_group1,test_group2,chain,everything()) %>%
  arrange(get(lev_cols),test_name,test_type,subset,test_group1,test_group2,chain)
write.table(total,paste(dir,last_dir,"_all_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)
write.table(total %>% filter(adjusted<0.05),paste(dir,last_dir,"_significant_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)
total$total.adj.FDR = p.adjust(total$pvalue,method="fdr")
total$total.adj.holm = p.adjust(total$pvalue,method="holm")
write.table(total %>% filter(total.adj.FDR<0.05),paste(dir,last_dir,"_overall_significant_tests_fdr.txt",sep=""),
            quote = F,sep="\t",row.names=F)
write.table(total %>% filter(total.adj.holm<0.05),paste(dir,last_dir,"_overall_significant_tests_holm.txt",sep=""),
            quote = F,sep="\t",row.names=F)

