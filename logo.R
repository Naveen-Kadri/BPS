library ('ggplot2')
library ('cowplot')
library ('ggseqlogo')



in_file=snakemake@input[["in_file"]]
plot_file=snakemake@output[["plot_file"]]

##in_file <- '/cluster/work/pausch/naveen/CNS/BPS/BPSII/formatted_bpp_res.txt'
##plot_file <- 'BP_logo.pdf'


##pdf (plot_file, height=4,width=20)
pdf (plot_file, height=6,width=8)

inf  <- read.table (in_file ,stringsAsFactors=F)
colnames (inf)  <-  c ('myid', 'bps', 'bps_pos','sc_bps', 'sc_ppt', 'sc', 'zsc', 'zsc_ppt', 'zsc')
##par (mar = c (6,4,4,4), las=1,cex.lab=1.5, cex.axis=1.5)
##ggseqlogo (inf$bps, method='bits', facet='wrap')
## nbases  <-  10
## toadd  <-  4
## mymat  <- matrix (nrow=nrow(inf), ncol=nbases+toadd)
## bases  <- c('A', 'C', 'G', 'T')

## for (k in 1:(nbases+toadd)  ){
##     set.seed (k)
##     mymat [,k] <-  sample (bases, nrow (inf), replace=T)
##     cat (k,"\n")
## }

## mypaste <- function (vec) {
##     paste(vec, collapse="")
## }

## dummy1  <- apply (mymat[, 1:nbases],1,mypaste)
## dummy2  <- apply (mymat[, 1:(nbases+toadd)],1,mypaste)







##myseq  <- paste (dummy1, inf$bps, dummy2,sep="")
##length(myseq)
msize <- 35
##ggplot() + theme_bw() +  geom_logo (myseq) +    theme(text = element_text(size=msize), axis.title.x=element_blank() ,  axis.text.x=element_blank() ,  axis.ticks.x=element_blank() )


##without grids
myplot <- ggplot () + geom_logo (inf$bps)
myplot + theme_bw() + theme(panel.grid.major = element_blank () , panel.grid.minor = element_blank () ,axis.text.x=element_blank (), axis.ticks.x=element_blank(), text=element_text (size=msize)) 
dev.off()
mytab <- table (substr (inf$bps,6,6))
100*(mytab /sum (mytab))
