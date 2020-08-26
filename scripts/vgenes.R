edgeR_general_test = function(cols,df,test.column,g1,g2){
  df = as.data.frame(df)
  cols_data = df %>% select(cols) %>% unique
  res = data.frame()
  for (i in c(1:nrow(cols_data))){
    i_line = cols_data %>% slice(i)
    df.i = merge(df,i_line,by=cols)
    res.i = eR(df.i,test.column,g1,g2)
    res = rbind(res,cbind(i_line,res.i))
  }
  return(res)
}

eR = function(df.vgene,test.column,g1,g2){
  variable = c(test.column)
  group1 = unique((df.vgene %>% filter(get(variable)==g2))$sample)
  group2 = unique((df.vgene %>% filter(get(variable)==g1))$sample)
  
  samples_n = length(unique(df.vgene$sample))
  bad_v = df.vgene %>%
    group_by(v) %>% summarise(c = length(unique(sample))) %>% filter(c<samples_n)
  #df.vgene.2 = df.vgene %>% filter(!(v %in% bad_v$v))
  df.vgene.2 = df.vgene
  temp = dcast(df.vgene.2,v ~ sample,value.var ="count" )
  temp[is.na(temp)]=0
  x = cbind(temp[,group1],temp[,group2])
  row.names(x)=temp[,1]
  group <- factor(c(rep(1,length(group1)),rep(2,length(group2))))
  y <- DGEList(counts=x,group=group)
  y <- calcNormFactors(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design,tagwise = TRUE)
  #print(y$common.dispersion)
  #print(summary(y$tagwise.dispersion))
  #plot(y$tagwise.dispersion)
  
  et <- exactTest(y,dispersion = "tagwise")
  #et <- exactTest(y, pair=c("1","2"),dispersion = 0.04)
  #fit <- glmQLFit(y,design)
  #qlf <- glmQLFTest(fit,coef=2)
  edge.res = as.data.frame(topTags(et, n = nrow(temp)))
  m = merge(df.vgene %>% select(v) %>% unique,
            edge.res %>% select(logFC,PValue),by.x ="v",by.y = "row.names",all=T)
  m$test_type = "edgeR"
  m$test_name = test.column
  m$test_group1 = g1
  m$test_group2 = g2
  return(m)
}

library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(edgeR,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
stat_source = args[1]
dir = args[2]
sample.file = args[3]
outdir = args[4]
#dir = "fastq_final_vdj/orig/"
#sample.file = "sample_info_2.2_corr.txt"
#stat_source="scripts_final/stat_test.util.R"


print("Start vgenes script..")
source(stat_source)
cap = "Significance: * - FDR<0.05, ** - FDR<0.01, *** - FDR<0.001"

file.pattern = 'vdj\\.corr\\.func\\.csd\\.down\\.'
files = list.files(dir,pattern = "vdj.corr.func.csd.down.[Ll]")
res.df = lapply(files,function(file.name,dir,file.pattern){
  f = fread(paste(dir,file.name,sep=""))
  prename = gsub(file.pattern,'',file.name)
  name = gsub('_.*','',prename)
  chain = gsub('.*_','',prename)
  chain = gsub('\\.txt','',chain)
  f$name=name
  f$chain=chain
  return(f)
},dir=dir,file.pattern) %>% rbindlist
samples = fread(sample.file)
samples$name = gsub('_.*','',samples$sample)
res.df.s = merge(res.df,samples,by=c("name","chain"))

df.vgene = res.df.s %>% group_by(sample,chain,subject_id,cell_type,status,status2,batch,v) %>%
  summarise(byReads = sum(count),byClones = n())
df.vgene = melt(df.vgene,measure.vars = c("byReads","byClones"),variable.name="w",value.name="count")
df.sums = df.vgene %>% group_by(sample,w,chain) %>% summarise(total.count = sum(count))
df.vgene = df.vgene %>% merge(df.sums,by = c("sample","w","chain"))
df.vgene = df.vgene %>% mutate(p = count/total.count)

test = dcast(df.vgene,chain+subject_id+status+status2+batch+v+w~cell_type,value.var ="p")
test = test %>% mutate(Breg = ifelse(is.na(Breg),0,Breg),
                       Bcell = ifelse(is.na(Bcell),0,Bcell))
test = test %>% mutate(Breg.norm = Breg/Bcell,
                       Breg.norm2 = Breg-Bcell)
test.m = melt(test,measure.vars = c("Bcell","Breg","Breg.norm","Breg.norm2"),
              variable.name = "cell_type")

## write data
res.data = df.vgene %>% mutate(variable = paste(w,v,sep=".")) %>% rename(value = p) %>%
  select(sample,chain,variable,value)
write.table(res.data,paste(outdir,"vgenes_data.txt",sep=""),quote = F,sep="\t",row.names=F)

##### wilcoxon
# status
stat.cols = c("chain","cell_type","v","w")
stat.res = wilcoxon_general_test(stat.cols,test.m,"status","MS","HD")

# cell type, paired
ct.cols = c("chain","v","w")
ct.res = wilcoxon_general_test(ct.cols,test.m,"cell_type","Breg","Bcell")
ct.res.pair = wilcoxon_pair_test(ct.cols,test,"cell_type","Breg","Bcell")
ct.res.all = rbind(ct.res,ct.res.pair)
ct.res.all$cell_type="all"

## edge R 
ct.cols2 = c("chain","w")
m.ct = edgeR_general_test(ct.cols2,df.vgene,"cell_type","Breg","Bcell")
m.ct = m.ct %>% mutate(cell_type="all")

m.batch = edgeR_general_test(ct.cols2,df.vgene,"batch",2,1)
m.batch = m.batch %>% mutate(cell_type="all")

stat.cols2 = c("chain","cell_type","w")
m.status = edgeR_general_test(stat.cols2,df.vgene,"status","MS","HD")

df.vgene.MS = df.vgene %>% filter(status == "MS")
m.status2 = edgeR_general_test(stat.cols2,df.vgene.MS,"status2","aggressive","benign")

### all pvalues table
edgeR.total = rbind(m.ct,m.status,m.status2,m.batch)
edgeR.total = edgeR.total %>% rename(estimate = logFC,pvalue=PValue)

wilcox.total = rbind(stat.res,ct.res.all)
table.totale = rbind(edgeR.total,wilcox.total)
table.totale = get.adjusted(c("test_type","test_name","cell_type","chain","w"),table.totale,"fdr")
table.totale = table.totale %>%
  mutate(v2 = ifelse((is.na(adjusted)|adjusted>0.05),v,
                     ifelse(adjusted>0.01,paste("*",v),
                            ifelse(adjusted>0.001,paste("**",v),paste("***",v)))))
# write tables
table.totale = table.totale %>% rename(subset = cell_type)
table.totale = table.totale %>% arrange(test_name,test_type,test_group1,test_group2,subset,chain,w,adjusted)

write.table(table.totale,paste(outdir,"vgenes_all_tests.txt",sep=""),quote = F,sep="\t",row.names=F)
write.table(table.totale %>% filter(adjusted<0.05),paste(outdir,"vgenes_significant_tests.txt",sep=""),quote = F,sep="\t",row.names=F)

############## PICS ###########

##### pics
## per subject pics

ws = c("byReads","byClones")
for (wi in ws){
  chains = c("IGH","IGK","IGL")
  for (ch in chains){
    df.vgene.chain = df.vgene %>% filter(chain==ch,w==wi)
    test.ch = test %>% filter(chain==ch,w==wi)
    vgenes.mean.p = df.vgene.chain %>% group_by(v) %>% summarise(mean.p = mean(p)) %>%
      arrange(desc(mean.p))
    top.vgenes = vgenes.mean.p[1:7,]$v
    test.ch = test.ch %>% 
      mutate(v2 = ifelse(v %in% top.vgenes,v,"other"),
             v.type = ifelse(v %in% top.vgenes,1,2))
    gv1 = ggplot(test.ch)+geom_point(aes(Bcell,Breg,color=v2,shape=factor(v.type)))+
      scale_shape_manual(values=c(19,1),guide="none")+
      facet_wrap(~subject_id+status)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste(outdir,"vgenes_per_subject_",wi,".",ch,".pdf",sep=""),gv1,height=12,width=12)
  }
  ## pics for Yasha
  norm.types = c("Breg.norm","Breg.norm2")
  test.w = test %>% filter(w==wi)  
  wilcox.status = table.totale %>% filter(test_type!="edgeR",test_name=="status",w==wi)
  for (t in norm.types){
    res.adj = wilcox.status %>% filter(subset==t) %>% arrange(chain,adjusted)
    g.ya = ggplot(test.w,aes(factor(v,levels = res.adj$v,
                                    labels = res.adj$v2),
                             get(t),color=status))+
      geom_boxplot()+facet_wrap(~chain,scales = "free",nrow = 3)+
      geom_point(position=position_dodge(width=0.75))+
      xlab("V gene")+
      theme_minimal()+labs(title = paste("MS vs HD: status wilcoxon test",t),
                           caption = cap)+
      theme(axis.text.x = element_text(angle = 90))
    if (t=="Breg.norm") g.ya = g.ya +scale_y_log10(limits=c(1e-02,100))
    ggsave(paste(outdir,"vgenes_status_",t,"_",wi,".pdf",sep=""),g.ya,height=12,width=18)
  }
}
## pics edgeR
edgeR.total = table.totale %>% filter(test_type=="edgeR")
cols = c("test_name","subset","test_group1","test_group2","w")
g.data = edgeR.total[,cols] %>% unique
df.vgene$batch = as.character(df.vgene$batch)
for(j in c(1:nrow(g.data))){
  data.j = g.data %>% slice(j)
  test_name = data.j[,1]
  subset = data.j[,2]
  wi =  data.j[,5]
  df.j = merge(edgeR.total,data.j,by=cols)
  df.j = df.j %>% arrange(chain,adjusted)
  df.j.s = df.j %>% group_by(chain) %>% arrange(adjusted) %>%
    slice(c(1:20)) %>% ungroup
  v20 = df.vgene %>% filter(v %in% df.j.s$v,w==wi)
  gv20 = ggplot(v20,aes(factor(v,levels = df.j.s$v,labels=df.j.s$v2),p,color=get(test_name)))+
    geom_boxplot()+
    geom_point(position=position_dodge(width=0.75))+
    facet_wrap(~chain,scales="free",nrow = 3)+scale_y_log10(limits=c(1e-05,1))+
    xlab("V gene (top 20, sorted by edgeR pvalue for difference)")+
    theme_minimal()+labs(title = paste(data.j[,3]," vs. ",data.j[,4],
                                       " (",test_name," test)",
                                       "\nsubset: ",subset,sep=""),
                         caption = cap)+
    theme(axis.text.x = element_text(angle = 90))
  ggsave(paste(outdir,"vgenes_",test_name,"_",subset,"_",wi,".top20.pdf",sep=""),gv20,height=12,width=8)
}