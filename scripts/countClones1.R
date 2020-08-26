library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
sample_info_file=args[1]
dir = args[2]
vdj_dir = args[3]
outdir = args[4]

get_value = function(text,lines){
  line = lines[grep(text,lines)]
  #print(line)
  if (length(line)==0) return(0)
  e = regexpr("[0-9]+",line)
  return(as.numeric(substr(line,e,(e+attr(e,"match.length")-1))))
}
getData = function(name,chain,dir,vdj_dir){
  # read1
  print(name)
  print(chain)
  print(dir)
  a1 = readLines(paste(dir,name,".align.report",sep=""))
  total.reads = get_value("Total sequencing reads:",a1)
  chain.reads.IGH = get_value("IGH chains:",a1)
  chain.reads.IGK = get_value("IGK chains:",a1)
  chain.reads.IGL = get_value("IGL chains:",a1)
  total.aligned = chain.reads.IGH+chain.reads.IGK+chain.reads.IGL
  chain.reads = get_value(paste(chain,"chains:"),a1)

  a2 = readLines(paste(dir,name,".",chain,".align.report",sep=""))
  chain.reads2 = get_value("Total sequencing reads:",a2)
  chain.reads3 = get_value(paste(chain,"chains:"),a2)
  
  ass = readLines(paste(dir,name,".",chain,".assemble.report",sep=""))
  clones.n1 = get_value("Final clonotype count:",ass)
  clones.reads1 = get_value("Reads used in clonotypes, percent of total:",ass)
  
  clones = fread(paste(dir,name,".",chain,".clones.txt",sep="")) 
  clones.n2 = nrow(clones)
  clones.reads2 = sum(clones$cloneCount)
  

  #check1
  reads1 = c(chain.reads,chain.reads2,chain.reads3)
  aln.reads = unique(reads1)
  if(length(aln.reads)>1) print(paste("Aln reads:",reads1))
  #reads2 = c(clones.reads1,clones.reads2,clones.reads3.1,clones.reads3.2,clones.reads4)
  reads2 = c(clones.reads1,clones.reads2)
  clone.reads = unique(reads2)
  if(length(clone.reads)>1) print(paste("Clone reads:",reads2))
  clones1 = c(clones.n1,clones.n2)
  clone.n = unique(clones1)
  if(length(clone.n)>1) print(paste("Clone n:",clones1))

  res1 = data.frame(sample=name,chain=chain,
                   total.reads=total.reads,total.aligned=total.aligned,
                   aln.reads=reads1[length(reads1)],
                   clones.reads=reads2[length(reads2)],clones.n=clones1[length(clones1)])
  
  # read2
  a1_2 = readLines(paste(dir,name,".",chain,".seed1.align.report",sep=""))
  chain.reads.seed1 = get_value("Total sequencing reads:",a1_2)
  chain.reads.seed2 = get_value(paste(chain,"chains:"),a1_2)
  
  ass_2 = readLines(paste(dir,name,".",chain,".seed1.assemble.report",sep=""))
  clones.seed.n1 = get_value("Final clonotype count:",ass_2)
  clones.seed.reads1 = get_value("Reads used in clonotypes, percent of total:",ass_2)
  
  clones.seed = fread(paste(dir,name,".",chain,".seed1.clones.txt",sep="")) 
  clones.seed.n2 = nrow(clones.seed)
  clones.seed.reads2 = sum(clones.seed$cloneCount)
  
  clones.cdr13 = fread(paste(dir,name,".",chain,".seed1.clones.withCdr13.txt",sep="")) 
  clones.cdr13.n = nrow(clones.cdr13)
  clones.cdr13.reads = sum(clones.cdr13$cloneCount)
  
  clones.decont = fread(paste(dir,name,".",chain,".seed1.clones.withCdr13.decont.txt",sep="")) 
  clones.decont.n = nrow(clones.decont)
  clones.decont.reads = sum(clones.decont$cloneCount)
  
  # check 2
  reads2_2 = c(clones.seed.reads1,clones.seed.reads2)
  clone.seed.reads = unique(reads2_2)
  if(length(clone.seed.reads)>1) print(paste("Clone reads 2:",reads2_2))
  clones1_2 = c(clones.seed.n1,clones.seed.n2)
  clone.seed.n = unique(clones1_2)
  if(length(clone.seed.n)>1) print(paste("Clone n 2:",clones1_2))
  
  res2 = data.frame(sampled.total.reads=chain.reads.seed1,sampled.aln.reads=chain.reads.seed2,
                    sampled.clones.reads=reads2_2[length(reads2_2)],
                    sampled.clones.n=clones1_2[length(clones1_2)],
                    sampled.clones.reads.withCdr13 = clones.cdr13.reads,
                    sampled.clones.n.withCdr13=clones.cdr13.n,
                    sampled.clones.reads.decont = clones.decont.reads,
                    sampled.clones.n.decont=clones.decont.n)
  
  # vdj data
  vdj0 = fread(paste(vdj_dir,"vdj.",name,"_",chain,".txt",sep="")) 
  vdj0.clones.n = nrow(vdj0)
  vdj0.reads = sum(vdj0$count)
  
  vdj1 = fread(paste(vdj_dir,"vdj.corr.",name,"_",chain,".txt",sep="")) 
  vdj1.clones.n = nrow(vdj1)
  vdj1.reads = sum(vdj1$count)
  
  vdj2 = fread(paste(vdj_dir,"vdj.corr.func.",name,"_",chain,".txt",sep="")) 
  vdj2.clones.n = nrow(vdj2)
  vdj2.reads = sum(vdj2$count)
  
  vdj3 = fread(paste(vdj_dir,"vdj.corr.func.csd.",name,"_",chain,".txt",sep="")) 
  vdj3.clones.n = nrow(vdj3)
  vdj3.reads = sum(vdj3$count)
  
  vdj4 = fread(paste(vdj_dir,"orig/vdj.corr.func.csd.down.",name,"_",chain,".txt",sep="")) 
  vdj4.clones.n = nrow(vdj4)
  vdj4.reads = sum(vdj4$count)
  
  res.vdj = data.frame(init.reads = vdj0.reads,init.clones = vdj0.clones.n,
                       corrected.reads = vdj1.reads,functional.reads = vdj2.reads,
                       crosssample.reads = vdj3.reads,
                       downsample.reads = vdj4.reads,
                       corrected.clones = vdj1.clones.n,functional.clones = vdj2.clones.n,
                       crosssample.clones = vdj3.clones.n,
                       downsample.clones = vdj4.clones.n)
  
  return(cbind(res1,res2,res.vdj))
}

sample_info = fread(sample_info_file) %>% select(sample,chain,cell_type,subject_id,batch) %>% unique
#sample_info = sample_info %>% filter(sample=="Lom8_S29_L001")
res = mapply(getData,name=sample_info$sample,chain=sample_info$chain,
       MoreArgs=list(dir=dir,vdj_dir=vdj_dir),SIMPLIFY = F) %>% rbindlist

res = merge(res,sample_info,by=c("sample","chain"))
cellN = fread("sample_cell_numbers.txt")
res = merge(res,cellN,by="subject_id")

#####
write.table(res,paste(outdir,"read_stat_data.txt",sep=""),
            quote = F,sep="\t",row.names=F)

# LINES
total.lines.reads = c("aln.reads","clones.reads","sampled.total.reads","sampled.aln.reads",
                "sampled.clones.reads","sampled.clones.reads.withCdr13",
                "sampled.clones.reads.decont","corrected.reads",
                "functional.reads","crosssample.reads","downsample.reads")
total.lines.clones = c("clones.n","sampled.clones.n","sampled.clones.n.withCdr13",
                       "sampled.clones.n.decont",
                       "corrected.clones","functional.clones",
                       "crosssample.clones","downsample.clones")
l.total = melt(res,id.vars = c("sample","chain","cell_type","batch"),
               measure.vars = total.lines.reads)
l.total.clones = melt(res,id.vars = c("sample","chain","cell_type","batch"),
               measure.vars = total.lines.clones)
tr = ggplot(l.total,aes(variable,value,color=chain))+
  geom_line(aes(group=sample))+
  facet_wrap(~chain)+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

tc = ggplot(l.total.clones,aes(variable,value,color=chain))+
  geom_line(aes(group=sample))+
  facet_wrap(~chain)+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## P's
#mixcr
res$aligned.p =  res$total.aligned/res$total.reads
res$clones.reads.p =  res$clones.reads/res$aln.reads
#sampling
res$clones.p =  res$sampled.clones.n/res$clones.n
res$reads.p =  res$sampled.clones.reads/res$clones.reads
#cdr13
res$sampled.cdr13.p.clones =  res$sampled.clones.n.withCdr13/res$sampled.clones.n
res$sampled.cdr13.p.reads =  res$sampled.clones.reads.withCdr13/res$sampled.clones.reads
#decontamination
res$decont.p.clones =  res$sampled.clones.n.decont/res$sampled.clones.n.withCdr13
res$decont.p.reads =  res$sampled.clones.reads.decont/res$sampled.clones.reads.withCdr13
#corr
res$vdj.corr.p.clones =  res$corrected.clones/res$init.clones
res$vdj.corr.p.reads =  res$corrected.reads/res$init.reads
#funct
res$vdj.func.p.clones =  res$functional.clones/res$corrected.clones
res$vdj.func.p.reads =  res$functional.reads/res$corrected.reads
#cross-sample
res$vdj.crosssample.p.clones =  res$crosssample.clones/res$functional.clones
res$vdj.crosssample.p.reads =  res$crosssample.reads/res$functional.reads
#downsample
res$vdj.downsample.p.clones =  res$downsample.clones/res$crosssample.clones
res$vdj.downsample.p.reads =  res$downsample.reads/res$crosssample.reads

ps1 = c("aligned.p","clones.reads.p",
        "sampled.cdr13.p.clones","sampled.cdr13.p.reads",
        "clones.p","reads.p",
        "decont.p.clones","decont.p.reads")
mp1 = melt(res,id.vars = c("sample","chain","cell_type","batch"),
         measure.vars = ps1)

ps.vdj = c("vdj.corr.p.clones","vdj.corr.p.reads",
           "vdj.func.p.clones","vdj.func.p.reads",
           "vdj.crosssample.p.clones","vdj.crosssample.p.reads",
           "vdj.downsample.p.clones","vdj.downsample.p.reads")
mp.vdj = melt(res,id.vars = c("sample","chain","cell_type","batch"),
           measure.vars = ps.vdj)

### mixcr align and sample
p1 = ggplot(mp1,aes(chain,value,color=as.factor(batch)))+
  geom_point(position=position_dodge(width=0.75))+
  facet_wrap(~variable,scales="free",nrow=2)
## vdj
p.vdj = ggplot(mp.vdj,aes(chain,value,color=as.factor(batch)))+
  geom_point(position=position_dodge(width=0.75))+
  facet_wrap(~variable,scales="free",nrow=2)

##### write #####

ggsave(paste(outdir,"read_stat_all_reads.pdf",sep=""),tr,height=5,width=8)
ggsave(paste(outdir,"read_stat_all_clones.pdf",sep=""),tc,height=5,width=8)
ggsave(paste(outdir,"read_stat_mixcr_p.pdf",sep=""),p1,height=5,width=8)
ggsave(paste(outdir,"read_stat_vdj_p.pdf",sep=""),p.vdj,height=5,width=8)


ggplot(res,aes(clones.reads,clones.n,color=cell_type))+
  geom_point()+facet_wrap(~chain,scales = "free")

ggplot(res %>% filter(cell_type=="Breg"),aes(cell_number,clones.n,color=cell_type))+
  geom_point()+facet_wrap(~chain,scales = "free")

chains = unique(res$chain)
cr = sapply(chains,function(ch){
  print(ch)
  res.ch = res %>% filter(chain==ch)
  print("reads/clones test")
  rc.corr = cor.test(res.ch$clones.n,res.ch$clones.reads)$p.value
  print(rc.corr)
  print("cell_n/clones test")
  res.ch = res.ch %>% filter(cell_type=="Breg")
  cc.corr = cor.test(res.ch$clones.n,res.ch$cell_number)$p.value
  print(cc.corr)
})


