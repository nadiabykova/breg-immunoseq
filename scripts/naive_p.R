init_data = function(file.name,dir,dir2,file.pattern){
  f = fread(paste(dir,file.name,sep=""))
  prename = gsub(file.pattern,'',file.name)
  name = gsub('_.*','',prename)
  chain = gsub('.*_','',prename)
  chain = gsub('\\.txt','',chain)
  #print(name)
  mut.file = list.files(dir2,pattern = paste(name,".*",chain,".seed1.cdr13.txt",sep=""))
  #print(mut.file)
  m = fread(paste(dir2,mut.file[1],sep=""))
  if (nrow(m[!is.na(m$mut.n),])!=nrow(m)) print("NAs in mutations file!")
  mf = merge(f,m,by.x="cdr3nt",by.y="nSeqCDR3")
  mf = mf %>% select(count,freq,cdr3nt,cdr3aa,v,d,j,VEnd,DStart,DEnd,JStart,mut.n)
  if (nrow(mf)!=nrow(f)) print("Hmm, lost ccdr13 clonotypes")
  #print(nrow(mf))
  if (nrow(mf)<10) {
    print("SMALL file!")
  }
  mf$name=name
  mf$chain=chain
  return(mf)
}

one_df = function(mf.init){
  data = mf.init %>% select(name,chain) %>% unique
  mf.res = data.frame()
  for (i in c(1:nrow(data))){
    mf = mf.init %>% merge(data[i,],by=c("name","chain"))
    naive.table = lapply(c(0:30),function(x){
      naive.n = mf %>% filter(mut.n<=x) %>% nrow
      naive.p = naive.n/nrow(mf)
      mf.mem = mf %>% filter(mut.n>x)
      ml = mean(mf.mem$mut.n)
      naive.sum.freq = sum((mf %>% filter(mut.n<=x))$freq)
      total.sum.freq = sum(mf$freq)
      naive.p.freq = naive.sum.freq/total.sum.freq
      res = data.frame(naive.t = x, naive.p = naive.p, naive.p.freq = naive.p.freq, 
                       read.count=sum(mf$count))
      return(res)
    }) %>% rbindlist()
    naive.table$name=unique(mf$name)
    naive.table$chain=unique(mf$ch)
    mf.res = rbind(mf.res,naive.table)
  }
  return(mf.res)
}

library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
stat_source = args[1]
dir = args[2]
dir2 = args[3]
sample.file = args[4]
outdir = args[5]

source(stat_source)
# after change
#dir = "fastq_final_vdj/orig/"
#dir2 = "fastq_final_results/"
#sample.file = "sample_info_2.2_corr.txt"

print("Start naive_p script..")
# init data
file.pattern = 'vdj\\.corr\\.func\\.csd\\.down\\.'
files = list.files(dir,pattern = "vdj.corr.func.csd.down.[Ll]")
print("Loading files..")
init.df = lapply(files,init_data,dir=dir,dir2=dir2,file.pattern) %>% rbindlist()
res.df = one_df(init.df)

# add samples
samples = fread(sample.file)
samples$name = gsub('_.*','',samples$sample)
res.df.s = merge(res.df,samples,by=c("name","chain"))
res.df.s = res.df.s %>% mutate(naive.p.freq = ifelse(subject_id=="MS3",NA,naive.p.freq))
#write.table(res.df.s %>% filter(naive.t<=3), paste(outdir,"naive_p_data3.txt",sep=""),quote = F,sep="\t",row.names=F)

## melt and dcast
res.m = melt(res.df.s,id.vars = c("sample","chain","naive.t","status","status2","batch","subject_id","cell_type"),
             measure.vars = c("naive.p","naive.p.freq"),variable.name = "w")
res.dc = dcast(res.m,chain+naive.t+status+status2+batch+w+subject_id~cell_type,value.var = "value")
res.dc = res.dc %>% mutate(Breg.norm=Breg/Bcell)
res.m2 = melt(res.dc,id.vars=c("chain","naive.t","status","status2","batch","w","subject_id"),
              variable.name="cell_type")

### write naive_p table
res.data = res.m %>% filter(naive.t<=3) %>% mutate(variable = paste(w,naive.t,sep=".t")) %>%
  select(sample,chain,variable,value)
write.table(res.data, paste(outdir,"naive_p_data.txt",sep=""),quote = F,sep="\t",row.names=F)
res.data2 = res.m %>% filter(naive.t<=3)
write.table(res.data2, paste(outdir,"naive_p_data2.txt",sep=""),quote = F,sep="\t",row.names=F)

######### P-VALUES ############
# pvalues: cell type
ct.cols = c("chain","naive.t","w")
ct.res = wilcoxon_general_test(ct.cols,res.m,"cell_type","Breg","Bcell")
ct.res.pair = wilcoxon_pair_test(ct.cols,res.dc,"cell_type","Breg","Bcell")

# pvalues: status 
stat.cols = c("chain","w","naive.t","cell_type")
stat.res = wilcoxon_general_test(stat.cols,res.m2,"status","MS","HD")
# pvalues: status2
stat2.res = wilcox_general_comb_test(stat.cols,res.m2,"status2")

# table totale
ct.res.all = rbind(ct.res,ct.res.pair)
ct.res.all$cell_type = "all"
table.totale = rbind(ct.res.all,stat.res,stat2.res)
table.totale = get.adjusted(c("test_name","test_type","w","naive.t"),table.totale)
table.totale = table.totale %>% rename(subset = cell_type)
table.totale = table.totale %>% arrange(test_name,test_type,w,subset,naive.t,chain)
table.totale = table.totale %>% filter(naive.t<=3)

write.table(table.totale,paste(outdir,"naive_p_all_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)
write.table(table.totale %>% filter(adjusted<0.05),paste(outdir,"naive_p_significant_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)

###### pics #### 
#color_subj = c("Evtushkina","Ovchinnikova","Chekanova","Lukmanov","additional","additional2")
#res.df.s = res.df.s %>% mutate(color = ifelse(subject_id %in% color_subj,subject_id,"others"))

w.cols = c("naive.p","naive.p.freq")
for (wi in w.cols){
  res.m.w = res.m %>% filter(w==wi)
  res.m2.w = res.m2 %>% filter(w==wi,naive.t==3)
  g.ct3 = ggplot(res.m.w %>% filter(naive.t==3),aes(cell_type,value,color=cell_type))+
    geom_boxplot()+geom_point()+
    facet_grid(~chain)+
    labs(title = paste("Fraction of naive clones",sep=""))
  
  g.ct = ggplot(res.m.w %>% filter(naive.t<=3),aes(cell_type,value,color=status))+
    geom_point()+geom_line(aes(group=subject_id))+
    facet_grid(chain~naive.t)+
    labs(title = paste("Fraction of naive clones (",wi,")",sep=""))
  g.ct_cum = ggplot(res.m.w,aes(naive.t,value,color=status,linetype=cell_type))+
    geom_line(aes(group=sample))+
    facet_grid(~chain)+
    labs(title = paste("Cumulative curve of naive clones (",wi,")",sep=""))
  g.stat = ggplot(res.m.w %>% filter(naive.t<=3),aes(cell_type,value,color=status))+
    geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
    facet_grid(chain~naive.t)+
    labs(title = paste("Naive colnes by status (",wi,")",sep=""))
  g.stat2 = ggplot(res.m.w %>% filter(naive.t<=3),aes(cell_type,value,color=status2))+
    geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
    facet_grid(chain~naive.t)+
    labs(title = paste("Naive colnes by status2  (",wi,")",sep=""))
  g.stat.withNorm.t3 = ggplot(res.m2.w,aes(status,value,color=status))+
    geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
    facet_grid(cell_type~chain,scales="free")+
    labs(title = paste("Naive clones by status (",wi,")",
                       "\nwith Breg.norm, naive.t = 3",sep=""))
  g.stat2.withNorm.t3 = ggplot(res.m2.w,aes(status2,value,color=status2))+
    geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
    facet_grid(cell_type~chain,scales="free")+
    labs(title = paste("Naive clones by status2 (",wi,")",
                       "\nwith Breg.norm, naive.t=3",sep=""))
  
  ggsave(paste(outdir,"naive_p_otchet2019_",wi,".pdf",sep=""),g.ct3,height=3,width=8)
  ggsave(paste(outdir,"naive_p_cell_type_",wi,".pdf",sep=""),g.ct,height=8,width=10)
  ggsave(paste(outdir,"naive_p_cell_type_CUM_",wi,".pdf",sep=""),g.ct_cum,height=6,width=10)
  ggsave(paste(outdir,"naive_p_status_",wi,".pdf",sep=""),g.stat,height=8,width=10)
  ggsave(paste(outdir,"naive_p_status2_",wi,".pdf",sep=""),g.stat2,height=8,width=10)
  ggsave(paste(outdir,"naive_p_status_withNorm_t3_",wi,".pdf",sep=""),g.stat.withNorm.t3,height=8,width=10)
  ggsave(paste(outdir,"naive_p_status3_withNorm_t3_",wi,".pdf",sep=""),g.stat2.withNorm.t3,height=8,width=10)
}

