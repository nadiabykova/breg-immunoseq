library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(gridExtra,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
stat_source = args[1]
dir = args[2]
outdir = args[3]
pics_dir = "pics/"

print("Start aa_features script..")
source(stat_source)
#outdir = "IBH_18102019"
#dir = "fastq_2.2_corr_vdj/orig/"
#stat_source = "stat_test.util.R"

regions = c("cdr3-full","cdr3-center-5")
aas = fread(paste(dir,"aas.cdr3aa.stat.wt.norm.txt",sep=""))
aas = aas %>% select(-..filter..,-q25,-q75)
aas2 = melt(aas,id.vars = c("sample_id","chain","cell_type","subject_id","status",
"status2","batch","region","property"),
            measure.vars = c("mean")) %>% 
  mutate(variable = paste(property,variable,sep="_"))
aas2$prop_type = "aa_feature"
aas2 = aas2 %>% filter(property!="count")

basic = fread(paste(dir,"basic.basicstats.txt",sep=""))
basic = basic %>% select(-..filter..,-nc_diversity,-nc_frequency)
basic.melt = melt(basic,id.vars = c("sample_id","chain","cell_type",
                                    "subject_id","status","status2","batch"))
basic.melt$region = "cdr3-full"
basic.melt$property = basic.melt$variable
basic.melt$prop_type = "basic"

df = rbind(aas2,basic.melt)
df = df %>% filter(region %in% regions)
df = rename(df, property = "variable",property_orig = "property")
df = df %>% select(-property_orig)

# Sidorov
df = df %>% filter(subject_id!="MS3")

## write data
res.data = df %>% mutate(sample = gsub("_IG[HKL]$","",sample_id),
  variable = paste(region,property,sep=".")) %>%
  select(sample,chain,variable,value)
write.table(res.data,paste(outdir,"aa_features_data.txt",sep=""),quote = F,sep="\t",row.names=F)

# status t-test
dc = dcast(df,chain+subject_id+status+status2+batch+region+prop_type+property~cell_type,
           value.var="value")
dc$Breg.norm=dc$Breg-dc$Bcell

m = melt(dc,measure.vars = c("Bcell","Breg","Breg.norm"),variable.name="cell_type")
stat.cols = c("chain","prop_type","property","region","cell_type")
prop.df = t_test_general_test(stat.cols,m,"status","MS","HD")

## status2: Anova test

prop.df2 = anova_general_test(stat.cols,m,"status2")

table.totale = rbind(prop.df,prop.df2)
table.totale = get.adjusted(c("test_name","test_type"),table.totale,"fdr")
table.totale = table.totale %>% rename(subset = cell_type)
table.totale = table.totale %>% arrange(test_name,test_type,subset,prop_type,property,region,chain)

table.selected = table.totale %>% filter(pvalue<0.05) %>% arrange(test_type,test_name,chain,property,region,subset)
#print(prop.df.selected)
table.selected2 = table.selected %>% group_by(test_type,test_name,chain,prop_type,property,region) %>%
  filter(length(unique(subset))>1,"Breg.norm" %in% unique(subset)) %>% ungroup %>%
  arrange(test_type,test_name,chain,property,region,subset)

#### write ####
write.table(table.totale,paste(outdir,"aa_features_all_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)
write.table(table.totale %>% filter(adjusted<0.05),paste(outdir,"aa_features_significant_tests.txt",sep=""),
            quote = F,sep="\t",row.names=F)

write.table(table.selected,paste(outdir,"aa_features_selected.txt",sep=""),
            quote = F,sep="\t",row.names=F)
write.table(table.selected2,paste(outdir,"aa_features_selected2.txt",sep=""),
            quote = F,sep="\t",row.names=F)

# pics
dir.create(paste(outdir,pics_dir,sep=""))
for (p in unique(m$property)){
  d = table.totale %>% filter(test_name=="status",property==p,pvalue<0.05)
  d = d %>% mutate(label = paste(chain,region,subset,
                                 format(round(pvalue, 3), nsmall = 3),
                                 format(round(adjusted, 3), nsmall = 3)))
  signif_lines = NULL
  if (nrow(d)>0) signif_lines = paste(d$label,collapse ="\n")
  m.p = m %>% filter(property==p,cell_type!="Breg.norm")
  type = unique(m.p$prop_type)
  g = ggplot(m.p,aes(cell_type,value,color=status))+
    geom_boxplot()+geom_point(position=position_dodge(width=0.75))+
    facet_wrap(region~chain,scales="free")+
    labs(title = paste(p,"(",type,")"),
         subtitle = signif_lines,
         caption = "Property values for each sample was calculated by VDJTools")
  #
  ggsave(paste(outdir,pics_dir,p,".pdf",sep=""),g,height=8,width=8)
}

