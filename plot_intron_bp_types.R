##in_file  <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN5/MUTATION/NEW/allinfo_mutation_count.txt"
##plot_file <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/scPPT/RUN5/NEW/intron_bp_types.pdf"

in_file  <- snakemake@input[["in_file"]]
plot_file  <- snakemake@output[["plot_file"]]

inf <- read.table (in_file, head=T,stringsAsFactors=F)

spl  <- strsplit (inf$bps.y, "" )
bptype <- paste(sapply (spl, "[", 4) , sapply (spl, "[", 6) , sep="N")
inttype  <- paste(inf$SS5, inf$SS3, sep="-")
comp  <- data.frame(int=inttype, bp=bptype)
head (comp)
mynum <- sort (table (comp$int), decreasing=T);head (mynum)



##look from the point of view of the BPs

inttab  <- sort (table (comp$int),decreasing=T)
toplot  <-  4
for (k in 1:toplot) {
    myint<- names (inttab [k])
    mytab  <- comp [comp$int == myint ,]
    n <- nrow (mytab)
    res  <- data.frame (head (100*sort (table (mytab$bp),decreasing=T) / n, toplot))
    colnames (res)  <-  c('bp', myint)
    if (k ==1) {
        allres  <- res
    }else {
        allres  <- merge (allres, res, by='bp', all=T)
        
    }
}

allres [is.na (allres)]   <- 0


toplot <- 4
toadd  <-  0.7
xlims <- c(0.8, toplot+toadd)
ylims  <- c(0,110)
xlabs <- 'Intron types'
ylabs <- '%'
##pdf ('Intron_bp_type.pdf',height=6, width=8)
pdf (plot_file, height=6, width=8)
mymat  <- matrix (1:2,ncol=1) ;mymat
layout (mymat, heights=c(0.85,0.15))
par (mar = c (4,6,2,2),cex.lab=1.5,cex.axis=1.5,las=1)
plot (1, type= "n", xlim=xlims, ylim=ylims, las=1, xlab=xlabs,ylab=ylabs,xaxt='n',yaxt='n')
rect (1:toplot, 0,  1:toplot + toadd, 100, col='lightgray', border='lightgray')
yats  <- seq (0,100,,length=5)
axis (side=2,at=yats,yats,las=1,cex.axis=1.5)
percentages  <- round (100* as.numeric (inttab) /sum (inttab),1)
##text (1:toplot + toadd /2, 105,  paste0(percentages [1:toplot], "%")  ,cex=1.5)
text (1:toplot + toadd /2, 105,  format(inttab[1:toplot], big.mark=",")  ,cex=1.2)


axis (side=1,1:toplot + toadd/2,names (inttab) [1:toplot] ,cex.axis=1.5)


#cols  <- c ('indianred2', 'darkseagreen', 'steelblue', 'violet', 'black' )
#cols  <- c ('darkseagreen', 'indianred2', 'steelblue', 'violet', 'black' )
cols  <- c('#5A76A5','#A5895A', 'indianred2', 'violet', 'black')

border_col <- 'orange'



prevres  <- rep (0, toplot)
for (k in 1:nrow (allres)) {
    myres  <- as.numeric (allres [k, -c(1)])
    rect (1:toplot, prevres, 1:toplot + toadd, myres+prevres, col=cols [k], lwd=0.5, border='whitesmoke')
    prevres  <- prevres + myres
}

#par (mar=c(0,6,0,4))

par (mar = c (0,6,0,2))
plot (1, type='n',xlab='',ylab='',axes=F)
legend ('center', horiz=T, col=cols, pch=15,cex=1, legend=(allres$bp), bty="o", title="BP types",pt.cex=1.5)

dev.off ()


