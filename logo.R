library ('ggplot2')
library ('cowplot')
library ('ggseqlogo')

in_file=snakemake@input[["in_file"]]
plot_file=snakemake@output[["plot_file"]]

pdf (plot_file, height=6,width=8)

inf  <- read.table (in_file ,stringsAsFactors=F)
colnames (inf)  <-  c ('myid', 'bps', 'bps_pos','sc_bps', 'sc_ppt', 'sc', 'zsc', 'zsc_ppt', 'zsc')

msize <- 35

##without grids
myplot <- ggplot () + geom_logo (inf$bps)
myplot + theme_bw() + theme(panel.grid.major = element_blank () , panel.grid.minor = element_blank () ,axis.text.x=element_blank (), axis.ticks.x=element_blank(), text=element_text (size=msize)) 
dev.off()
mytab <- table (substr (inf$bps,6,6))
100*(mytab /sum (mytab))
