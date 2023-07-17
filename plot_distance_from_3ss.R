##in_file  <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN2/MUTATION/NEW/allinfo_mutation_count.txt"
in_file  <- snakemake@input[["in_file"]]
plot_file  <- snakemake@output[["plot_file"]]


##in_file  <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN5/MUTATION/NEW/allinfo_mutation_count.txt"
##plot_file  <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN5/NEW/distance_from_3ss.pdf"

inf <- read.table (in_file, head=T)
##pdf ('distance_from_3ss.pdf', height=6,width=8)
pdf (plot_file, height=6,width=8)

par (mar = c (6,6,4,4),cex.lab=1.5, cex.axis=1.5)
myseq <- seq (min(inf$bp_pos), max(inf$bp_pos), by=1)
all <- hist (inf$bp_pos, breaks=myseq,plot=F)
all$counts  <- 100*all$counts/sum (all$counts)

ylims  <- c (0,7)
ylabs <- "Variants per 100 bp"
ylabs  <- '%'
plot (all ,las=1,main="", xlab="Distance from 3' SS", ylab=ylabs, col='gray',ylim=ylims, border='whitesmoke')

can <- inf [inf$type=="canonical",]
can <- hist (can$bp_pos, breaks=myseq,plot=F)
can$counts  <- 100*can$counts/sum (can$counts)

noncan <- inf [inf$type=="noncanonical",]
noncan <- hist (noncan$bp_pos, breaks=myseq, plot=F)
noncan$counts  <- 100*noncan$counts/sum (noncan$counts)
#cols  <- c('indianred2', 'darkseagreen')
#cols  <- c('#5A76A5', '#A5895A')
cols  <- c('darkseagreen', 'indianred2')
par (new =T )
plot (can, col=cols[1],xlab='',ylab='',main='',ylim=ylims,las=1, border="whitesmoke")

par (new =T )
plot (noncan, col=cols [2],xlab='',ylab='',main='',ylim=ylims,las=1, border="whitesmoke")

par (new=T)
new <- can
new$counts [new$counts > noncan$counts] =0 
plot (new, col=cols[1],xlab='',ylab='',main='',ylim=ylims,las=1, border="whitesmoke")

#legend ("topright", col=cols, legend=c("canonical BP", "Non canonical BP"),pch=15, bty="n",cex=1.5)
legend ("topright", col=cols, legend=c("TNA", "Non TNA"),pch=15, bty="n",cex=1.5)

dev.off ()



