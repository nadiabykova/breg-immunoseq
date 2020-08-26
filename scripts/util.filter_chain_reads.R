args = commandArgs(TRUE)
id_file = args[1]
fasta_file = args[2]
new_fasta_file = args[3]

ids.df = read.csv(id_file,header=T,sep="\t")
ids = ids.df$readId
new_line_ids = as.vector(mapply(seq,4*ids+1,4*ids+4))

conn <- file(fasta_file,open="r")
lines <-readLines(conn)
lines_IGH = lines[new_line_ids]
close(conn)

unlink(new_fasta_file)
conn <- file(new_fasta_file,open="w")
writeLines(lines_IGH,conn)
close(conn)