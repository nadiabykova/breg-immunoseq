args = commandArgs(TRUE)
dir = args[1]
chain = args[2]
sample_info = args[3]
seed = as.integer(args[4])
outdir = args[5]

#dir = "/Users/koshka/Documents/fav_new/final"
#chain = "IGH"
#seed=1
#sample_info = "/Users/koshka/Documents/fav_new/final/sample_info_1.txt"
#tun = ".cdr3"

tun = ""
si = read.csv(sample_info,sep="\t",header=T)
si = si[si$chain == chain,]
sizes = mapply(function(name,chain){
  i1 = paste(dir,"/",name,".",chain,tun,".reads_R1.fastq",sep="")
  i2 = paste(dir,"/",name,".",chain,tun,".reads_R2.fastq",sep="")
  conn <- file(i1,open="r")
  lines1 <-readLines(conn)
  close(conn)
  conn <- file(i2,open="r")
  lines2 <-readLines(conn)
  close(conn)
  if (length(lines1)!=length(lines2)) print(paste("Unequal fastq files:",i))
  size = length(lines1)/4
  print(paste(name,size))
  return(size)
},name = si$sample,chain = si$chain)
min_size = min(sizes)

df = data.frame(name = si$sample,size = sizes, sample_size = min_size)
write.table(df,file = paste(outdir,"/",chain,tun,".seed",seed,".sample_fasta.report",sep=""),
            row.names = F, sep="\t",quote=F)

# sample
mapply(function(name,chain){
  i1 = paste(dir,"/",name,".",chain,tun,".reads_R1.fastq",sep="")
  i2 = paste(dir,"/",name,".",chain,tun,".reads_R2.fastq",sep="")
  conn <- file(i1,open="r")
  lines1 <-readLines(conn)
  close(conn)
  conn <- file(i2,open="r")
  lines2 <-readLines(conn)
  close(conn)
  ids = c(1:(length(lines1)/4))
  set.seed(seed)
  ids_sampled = sample(ids,min_size)
  new_line_ids = as.vector(mapply(seq,4*ids_sampled-3,4*ids_sampled))
  sampled_lines_1 = lines1[new_line_ids]
  sampled_lines_2 = lines2[new_line_ids]
  i1_s = paste(outdir,"/",name,".",chain,tun,".sampled_seed",seed,"_R1.fastq",sep="")
  i2_s = paste(outdir,"/",name,".",chain,tun,".sampled_seed",seed,"_R2.fastq",sep="")
  unlink(i1_s)
  conn <- file(i1_s,open="w")
  writeLines(sampled_lines_1,conn)
  close(conn)
  unlink(i2_s)
  conn <- file(i2_s,open="w")
  writeLines(sampled_lines_2,conn)
  close(conn)
},name = si$sample,chain = si$chain)

