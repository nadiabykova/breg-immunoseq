library(rjson)
library(jsonlite)
imgt.lib = fromJSON("imgt.201822-5.sv4.json") %>% as.data.frame %>% filter(taxonId==9606)
imgt.genes = imgt.lib$genes[[1]]
imgt.seqs = imgt.lib$sequenceFragments[[1]]
getRegionBoundary = function(gene,reg){
  if (reg=="CDR1"){
    b1="CDR1Begin"
    b2="FR2Begin"
  } else if (reg=="FR2"){
    b1="FR2Begin"
    b2="CDR2Begin"
  } else if (reg=="CDR2"){
    b1="CDR2Begin"
    b2="FR3Begin"
  } else if (reg=="FR3"){
    b1="FR3Begin"
    b2="CDR3Begin"
  }
  res = (imgt.genes %>% filter(name==gene) %>% select(anchorPoints))[1,1] %>%
    select(one_of(c(b1,b2)))
  colnames(res) = c("begin","next.begin")
  return(res)
}
getRegionGermline = function(gene,reg){
  z = (imgt.genes %>% filter(name==gene))$baseSequence
  res = imgt.seqs %>% filter(uri==z)
  from = (res$range)$from
  bonds = getRegionBoundary(gene,reg)
  seq = res$sequence
  return(substr(seq,bonds[1,1]+1-from,bonds[1,2]-from))
}
mut.reg = function(mut.line,begin,next.begin){
  #mut.line= "SG71T,SA89T,SG123T,SA141G"
  #mut.line = "I75A,I75C,ST77A,ST78A,SC82T,ST84G,ST86C,DC88,SA90T,DC92,ST93A,SG96T,I99C,I99C"
  mut.a = strsplit(mut.line,",")[[1]]
  pos = sapply(mut.a, function(x){as.numeric(gsub("[SDIATGC]","",x))})
  res = mut.a[which(pos>=begin&pos<next.begin)]
  edge.insert1 = paste("I",begin,sep="")
  edge.insert2 = paste("I",next.begin,sep="")
  e1 = grep(edge.insert1,res)
  if (length(e1)>0) res = res[-e1]
  e2 = grep(edge.insert2,mut.a)
  return(paste(c(res, mut.a[e2]),collapse=","))
}
parseMut = function(aln.line){
  #aln.line = "98|297|316|2|201|SG101TSA103TSG123TSA141GSA146GSC152GSA153GST164AST167GSA171GSA178TSA182TSG194ASG196ASG215AST222GSC236GSC239TSG245ASG259CSG266ASG278ASA283TSC284TSG290T|645.0"
  x = strsplit(aln.line,"\\|")[[1]][6]
  x.expr = gregexpr("[SDI][ATGC]*[0-9]+[ATGC]*",x)[[1]]
  x.array = mapply(substr,start=x.expr,stop=x.expr+attr(x.expr,"match.length")-1,
                   MoreArgs = list(x=x))
  return(paste(x.array,collapse = ","))
}
targets.parseMut = function(targets.line){
  x = strsplit(targets.line,",")[[1]]
  if(length(x)==1) x[2]=""
  res = sapply(x,parseMut)
  return(paste(res,collapse = ";"))
}
anchor.seq = function(seq,anchor.line,region){
  if (region=="CDR1"){
    start.ind=6
  } else if (region=="FR2"){
    start.ind=7
  } else if (region=="CDR2"){
    start.ind=8
  } else if (region=="FR3"){
    start.ind=9
  }
  #anchor.line=":::::69:93:144:168:279:11:300::::::::::"
  idx = strsplit(anchor.line,":")[[1]]
  start = idx[start.ind]
  end = idx[start.ind+1]
  if (start!=""&end!=""){
    return(substr(seq,(as.numeric(start)+1),as.numeric(end)))
  } else return(NA) 
}
check.mut = function(mut,seq,germline,begin,next.begin){
  l = nchar(germline)
  if (next.begin-begin!=l) {
    print("germline length fail")
  }
  mut.a = strsplit(mut,",")[[1]]
  shift=0
  res = strsplit(germline,"")[[1]]
  for (m in mut.a){
    type = substr(m,1,1)
    if (type=="S"){
      from = substr(m,2,2)
      pos = as.numeric(substr(m,3,nchar(m)-1))
      to = substr(m,nchar(m),nchar(m))
      pos.local = pos-begin+shift+1
      if (res[pos.local]!=from) {
        print(m)
        print("S not good")
      }
      res[pos.local]=to
    } else if (type=="D"){
      from = substr(m,2,2)
      pos = as.numeric(substr(m,3,nchar(m)))
      pos.local = pos-begin+shift+1
      if (res[pos.local]!=from) {
        print(m)
        print("D not good")
      }
      res = res[-pos.local]
      shift = shift - 1
    } else{
      if (type!="I") print("wft")
      pos = as.numeric(substr(m,2,nchar(m)-1))
      to = substr(m,nchar(m),nchar(m))
      pos.local = pos-begin+shift+1
      if (pos.local>length(res)){
        res = c(res[1:(pos.local-1)],to)
      } else {
        res = c(res[1:(pos.local-1)],to,res[pos.local:length(res)]) 
      }
      shift = shift + 1
    }
  }
  res.seq = paste(res,collapse = "")
  if (seq==res.seq) out = TRUE
  else {
    out = FALSE
  }
  return(out)
}

fix.mutations = function(f){
  n1 = nrow(f)
  f = f %>% filter(!(nSeqCDR1==""|nSeqFR2==""|nSeqCDR2==""|nSeqFR3==""|nSeqCDR3==""))
  
  if (nrow(f)!=n1) print("Not all readsd has all regions!!!")
  #print(nrow(f))
  #print(table(f$numberOfTargets))
  
  fix = f %>% select(readId,numberOfTargets,bestVHit,nMutationsInCDR1RelativeVRegion,nSeqCDR1,
                     nMutationsInFR2RelativeVRegion,nSeqFR2,
                     nMutationsInCDR2RelativeVRegion,nSeqCDR2,
                     nMutationsInFR3RelativeVRegion,nSeqFR3)
  m.fix = melt(fix,id.vars = c("readId","bestVHit","numberOfTargets"),variable.name="col.name",value.name="data") %>%
    mutate(type = ifelse(grepl("nSeq",col.name),"seq","mutations"),
           region = gsub("nMutationsIn|RelativeVRegion|nSeq","",col.name)) %>%
    select(-col.name)
  dc.fix = dcast(m.fix,readId+numberOfTargets+bestVHit+region~type,value.var="data")
  
  genes = dc.fix %>% select(bestVHit,region) %>% unique
  germlines = genes %>% mutate(germline=mapply(getRegionGermline,gene=genes$bestVHit,reg=genes$region))
  germline.data = cbind(germlines,mapply(getRegionBoundary,gene=germlines$bestVHit,reg=germlines$region,SIMPLIFY = F) %>% rbindlist)
  
  dc.fix2 = merge(dc.fix,germline.data,by=c("bestVHit","region")) %>% 
    mutate(mut.match = mapply(check.mut,mut=mutations,seq = seq,
                              germline = germline,begin = begin,next.begin = next.begin))
  
  # check that problems only when reads N = 2
  table(dc.fix2 %>% select(region,mut.match,numberOfTargets))
  test1 = dc.fix2 %>% filter(numberOfTargets==1) %>% select(mut.match) %>% unique
  if (nrow(test1)>1|test1[1,1]!=TRUE) print("Test 1 not passed! Wrong mutations when reads N = 1")
  
  # fix mutations in 2
  f2 = f %>% filter(numberOfTargets==2) %>% select(readId,bestVAlignment) %>%
    mutate(read.mutations.total = sapply(bestVAlignment,targets.parseMut)) %>% 
    mutate(read1.mutations.total = sapply(read.mutations.total,function(x){strsplit(x,";")[[1]][1]}),
           read2.mutations.total = sapply(read.mutations.total,function(x){strsplit(x,";")[[1]][2]}))
  
  dc.fix3 = merge(dc.fix2,f2 %>% select(readId,read1.mutations.total,read2.mutations.total),by="readId") %>%
    mutate(read1.mutations = mapply(mut.reg, mut.line = read1.mutations.total, begin = begin, next.begin = next.begin),
           read2.mutations = mapply(mut.reg, mut.line = read2.mutations.total, begin = begin, next.begin = next.begin))
  
  dc.fix3 = dc.fix3 %>% mutate(read1.mut.match = mapply(check.mut,mut = read1.mutations, seq = seq,
                                                        germline = germline,begin = begin, next.begin = next.begin),
                               read2.mut.match = mapply(check.mut,mut = read2.mutations, seq = seq,
                                                        germline = germline, begin = begin, next.begin = next.begin))
  # check that ALL not matched are matched in read #2
  table(dc.fix3 %>% select(read2.mut.match,region,mut.match),useNA="ifany")
  test2 = dc.fix3 %>% filter(mut.match==F) %>% select(read2.mut.match,read1.mut.match) %>% table(useNA="ifany")
  if (nrow(test2)>1|rownames(test2)!=TRUE|colnames(test2)!=FALSE) print("Test 2 not passed! Unmatched cases have no matches in read 2")
  
  # finally, fix
  dc.fix3 = dc.fix3 %>% mutate(new.mutations = ifelse(mut.match,mutations,read2.mutations))                         
  dc.fix3 = dc.fix3 %>% mutate(new.mut.match = mapply(check.mut, mut = new.mutations, seq = seq,
                                                      germline = germline, begin = begin, next.begin = next.begin))
  table(dc.fix3 %>% select(region,new.mut.match))
  test3 = dc.fix3 %>% select(new.mut.match) %>% unique
  if (nrow(test3)>1|test3[1,1]!=TRUE) print("Test 3 not passed! Wrong mutations after fix :(")
  
  # convert back
  
  dc.fixed = rbind(dc.fix3 %>% select(readId,region,new.mutations,seq) %>% rename(mutations = new.mutations),
                   dc.fix2 %>% filter(numberOfTargets==1) %>% select(readId,region,mutations,seq))
  
  m.fixed = melt(dc.fixed,id.vars = c("readId","region"),variable.name="type",value.name="data") %>%
    mutate(name = ifelse(type=="seq",paste("nSeq",region,sep=""),paste("nMutationsIn",region,"RelativeVRegion",sep=""))) %>%
    select(-type,-region)
  fixed = dcast(m.fixed,readId~name,value.var="data")
  
  f.fixed = f %>% select(-nMutationsInCDR1RelativeVRegion,-nMutationsInCDR2RelativeVRegion,
                         -nMutationsInFR2RelativeVRegion,-nMutationsInFR3RelativeVRegion) %>%
    merge(fixed %>% select(readId,nMutationsInCDR1RelativeVRegion,nMutationsInCDR2RelativeVRegion,
                           nMutationsInFR2RelativeVRegion,nMutationsInFR3RelativeVRegion),by="readId")
  
  return(f.fixed)
}
