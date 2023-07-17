info_file <- snakemake@input[["info_file"]]
mutation_file<- snakemake@input[["mutation_file"]]
res_file <-  snakemake@input[["bpp_res"]]
out_file  <- snakemake@output[["out_file"]]
myrange  <- snakemake@params[["myrange"]]
myrange  <- seq (myrange [1], myrange [2]-1)


#info_file ="/cluster/work/pausch/naveen/CNS/BPS/BPSII/seq_info.txt"
## mutation_file  <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/test.txt"
## res_file  <- "/cluster/work/pausch/naveen/CNS/BPS/BPSII/bpp_res.txt"
## ## out_file <-  "todel.txt"




res  <- read.table (res_file, stringsAsFactors=F)
cnames  <- c("id", "bps", "bp_pos", "sc_bps", "sc_ppt", "sc", "zsc_bps", "zsc_ppt", "zsc")
colnames (res)  <- cnames
head (res)

info  <- read.table (info_file, head=T, stringsAsFactors=F)
info$intron_nr  <-  as.numeric (substr (info$intron,7,nchar(info$intron)))
info$id  <- paste (info$transcript,info$intron_nr, info$chr, info$strand, info$start, info$end, sep="_")
info$id  <-  paste0(">", info$id)
allinfo  <- merge (res, info, by="id")




mut_res  <- read.table (mutation_file, stringsAsFactors=F)
colnames (mut_res)  <- c ('cord', 'id', 'bps', 'type', paste0("dist",myrange))
fres  <- merge (mut_res, allinfo, by="id")
write.table (fres, out_file, col.names=T,row.names=F,quote=F,sep="\t")



## dim (fres)








## nrow(fres)
## nrow (mut_res)

## ##maxi <- max (apply (inf [, 4:ncol (inf)] ,2, sum) / nrow (inf))
## muts  <- apply( inf[, 5:ncol (inf)],2,sum)
## data.frame (muts)
## plot (-51:20,muts, lty=1, lwd=3, type="h", col="green")




## can  <- inf [inf$type == "canonical",]
## x  <- apply( can[, 5:ncol (inf)],2,sum)
## plot (-51:20,x, lty=1, lwd=3, type="h", col="green")

## ncan  <- inf [inf$type == "noncanonical",]
## y  <- apply( ncan[, 5:ncol (inf)],2,sum)
## plot (-51:20,y, lty=1, lwd=3, type="h", col="green")



## dist_cols <- grep ('dist', colnames (allinfo))


## canss  <- allinfo [allinfo$SS5== "GT" & allinfo$SS3=="AG",]
## sum (canss$type=="canonical") /nrow (canss)
## sum (allinfo$type =="canonical")/ nrow (allinfo)






## z  <- apply(canss[, dist_cols],2,sum)
## mycol  <- "#DA8D25"
## plot (-51:20,z, lty=1, type="h", col=mycol,lwd=10)
## abline (v = 4, col='red')

## mydata  <- allinfo[allinfo$strand == "+",]
## z  <- apply(mydata[, dist_cols],2,sum)
## mycol  <- "#DA8D25"
## plot (-51:20,z, lty=1, type="h", col=mycol,lwd=10)
## abline (v = 4, col='red')




## old  <- read.table ("/cluster/work/pausch/naveen/CNS/BPS/ALL/mutation_count.txt")
## colnames (old)  <- c ('cord', 'bps', 'type', paste0("dist",-51:20))
## old <- apply (old [, 4:ncol (inf)],2,sum)
## plot(-51:20,old, lty=1, lwd=3, type="b", col='red')


## abline (v = c (4,6), col="red")
## x  <- data.frame (x)


## x  <- inf [, c("cord", "bps", "type", "dist4", "dist6") ]
## xcan <- x [x$type =="canonical",]
## xnoncan <- x [x$type =="noncanonical",]

## 100*sum (xcan$dist6 ==1) / nrow (xcan)
## 100*sum (xnoncan$dist6 ==1) / nrow (xnoncan)


## sum (xcan$dist6 ==1) 
## sum (xnoncan$dist6 ==1)

## 141/689

## 100*sum (xcan$dist4 ==1) / nrow (xcan)
## 100*sum (xnoncan$dist4 ==1) / nrow (xnoncan)




## can <- inf [inf$type =='canonical',]
## percan <- 100*round (nrow (can)/nrow (inf),4)
## noncan  <- inf [inf$type =='noncanonical',]


## ylims  <- c (0,maxi)

## myx <-  -51:20
## n=length (myx)
## pdf (plot_file, height=12,width=20)
## par (mar = c (6,6,4,4),cex.lab=1.5, cex.axis=1.5)
## plot (1, type='n', xlim=c(0,n), ylim=ylims,xlab='',ylab='',las=1,xaxt="n")
## #axis (side=1, 1:n, labels=myx)
## mycols  <- rep ('indianred2', n)
## mycols [myx<=0] ="darkseagreen"
## mycols [myx>7] ="darkseagreen"

## rect(1:n,  0, 1:n+0.80,apply (inf[,4:ncol(inf)],2,sum)/nrow(inf), col=mycols, border=mycols)

## mycols  <- rep ('firebrick1', n)
## mycols [myx<=0] ="darkseagreen4"
## mycols [myx>7] ="darkseagreen4"

## rect(1:n,  0, 1:n+0.80,apply (can[,4:ncol(can)],2,sum)/nrow(can), col=mycols, border=mycols)
## legend ('topright',legend= paste0(percan, "%"),cex=1,bty="n")
## genome_wide <- nvar/genome_length
## abline (h=genome_wide, col='red', lty=3)
## dev.off ()

