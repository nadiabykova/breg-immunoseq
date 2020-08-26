library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(gridExtra,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
stat_source = args[1]
dir = args[2]
outdir = args[3]

#dir = "fastq_2.2_corr_vdj/orig/"
print("Start clonality script..")
source(stat_source)
#file2 = paste(dir,"stat.diversity.strict.resampled.txt",sep="")
file2 = paste(dir,"stat.diversity.strict.exact.txt",sep="")
f2 = fread(file2)
f2 = f2 %>% filter(subject_id!="MS3")

## write data
res.data = f2 %>% mutate(sample = gsub("_IG[HKL]$","",sample_id)) %>%
  melt(measure.var="normalizedShannonWienerIndex_mean") %>%
  select(sample,chain,variable,value)
write.table(res.data,paste(outdir,"clonality_data.txt",sep=""),quote = F,sep="\t",row.names=F)

## pvalue ##
# cell_type

ct.cols = c("chain")
f2.test = f2 %>% rename(value = normalizedShannonWienerIndex_mean)
ct.res = wilcoxon_general_test(ct.cols,f2.test,"cell_type","Breg","Bcell")

f2.test.dc = dcast(f2.test,chain+status+status2+batch+subject_id~cell_type,value.var="value")
ct.res.pair = wilcoxon_pair_test(ct.cols,f2.test.dc,"cell_type","Breg","Bcell")

#ct.res.adj = get.adjusted(c("test_name","test_type"),ct.res.all)
#ct.res.adj = ct.res.adj %>% select(test_name,test_type,test_group1,test_group2,cell_type,
#                                  chain,estimate,pvalue,adjusted,adj.method) %>%
#  arrange(test_type,test_name,cell_type,chain)

# status
stat.cols = c("chain","cell_type")
stat.res = wilcoxon_general_test(stat.cols,f2.test,"status","MS","HD")
#stat.res.adj = get.adjusted(c("test_name"),stat.res)
#stat.res.adj = stat.res.adj %>% select(test_name,test_type,test_group1,test_group2,cell_type,
#                                       chain,estimate,pvalue,adjusted,adj.method) %>%
#  arrange(test_type,test_name,cell_type,chain)

# pvalues: status2
stat2.res = wilcox_general_comb_test(stat.cols,f2.test,"status2")
#stat2.res.adj = get.adjusted(c("test_name"),stat2.res)
#stat2.res.adj = stat2.res.adj %>% select(test_name,test_type,test_group1,test_group2,cell_type,
#                                       chain,estimate,pvalue,adjusted,adj.method) %>%
#  arrange(test_type,test_name,cell_type,chain)

#table totale
ct.res.all = rbind(ct.res,ct.res.pair)
ct.res.all$cell_type = "all"
table.totale = rbind(ct.res.all,stat.res,stat2.res)
table.totale = table.totale %>% rename(subset = cell_type)
table.totale = get.adjusted(c("test_name","test_type"),table.totale)
table.totale = table.totale %>% arrange(test_name,test_type,subset,chain)

## write
write.table(table.totale,paste(outdir,"clonality_all_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)
write.table(table.totale %>% filter(adjusted<0.05),paste(outdir,"clonality_significant_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)

## write data


### pictures ####
levels = unique(f2$chain) %>% sort
labels = sapply(levels, function(x){
  p = table.totale %>% filter(test_name=="cell_type",chain==x,test_type=="Wilcoxon")
  return(paste(x,"\nWilk p.val=",format(round(p$pvalue, 3), nsmall = 3),
         "\n(Holm adj: ",format(round(p$adjusted, 3), nsmall = 3),")",sep=""))
})
labels.pair = sapply(levels, function(x){
  p = table.totale %>% filter(test_name=="cell_type",chain==x,test_type=="Wilcoxon paired")
  return(paste(x,"\npair Wilk p.val=",format(round(p$pvalue, 3), nsmall = 3),
               "\n(Holm adj: ",format(round(p$adjusted, 3), nsmall = 3),")",sep=""))
})
f2$chain.f <- factor(f2$chain,levels =levels,labels = labels)
f2$chain.f.pair <- factor(f2$chain,levels = levels,labels = labels.pair)

co = ggplot(f2,aes(cell_type,normalizedShannonWienerIndex_mean,color=cell_type))+
  geom_boxplot()+geom_point()+
  facet_wrap(~chain)

c1 = ggplot(f2,aes(cell_type,normalizedShannonWienerIndex_mean,color=cell_type))+
  geom_boxplot()+geom_point()+
  facet_wrap(~chain.f)
c2 = ggplot(f2,aes(cell_type,normalizedShannonWienerIndex_mean,color=as.factor(batch)))+
  geom_line(aes(group=subject_id))+geom_point()+
  facet_wrap(~chain.f.pair)
g_cell_type = grid.arrange(c1,c2,nrow=2)

# status
levels = unique(f2$chain) %>% sort
labels = sapply(levels, function(x){
  p.bcell = table.totale %>% filter(test_name=="status",chain==x,subset == "Bcell")
  p.breg = table.totale %>% filter(test_name=="status",chain==x,subset == "Breg")
  return(paste(x,"\nWilk.bcell p.val=",format(round(p.bcell$pvalue, 3), nsmall = 3),
               "\n(Holm adj: ",format(round(p.bcell$adjusted, 3), nsmall = 3),")",
               "\nWilk.breg p.val=",format(round(p.breg$pvalue, 3), nsmall = 3),
               "\n(Holm adj: ",format(round(p.breg$adjusted, 3), nsmall = 3),")",
               sep=""))
})
f2$chain.fs <- factor(f2$chain,levels =levels,labels = labels)

g_status = ggplot(f2,aes(cell_type,normalizedShannonWienerIndex_mean,color=status))+
  geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
  facet_wrap(~chain.fs)

# status2
g_status2 = ggplot(f2,aes(cell_type,normalizedShannonWienerIndex_mean,color=status2))+
  geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
  facet_wrap(~chain)

#
ggsave(paste(outdir,"clonality_otchet2019.pdf",sep=""),co,height=3,width=8)

ggsave(paste(outdir,"clonality_cell_type.pdf",sep=""),g_cell_type,height=8,width=8)
ggsave(paste(outdir,"clonality_status.pdf",sep=""),g_status,height=5,width=8)
ggsave(paste(outdir,"clonality_status2.pdf",sep=""),g_status2,height=5,width=8)

