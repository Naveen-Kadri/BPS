infiles = snakemake@input[["infiles"]]
outfile=snakemake@output[["outfile"]]

chr <- 0
for (myfile in infiles){
    chr  <- chr +1
    inf <- read.table (myfile,head=T)
    inf$chr  <- chr
    if (chr ==1) {
        allchr  <- inf
    }else {
        allchr  <- rbind (allchr, inf)
    }
}

write.table (allchr, outfile, col.names=T,row.names=F,quote=F,sep="\t")
