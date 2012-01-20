#!/usr/bin/env Rscript

# This script takes the results of a differential expression analysis and a mapping to peak data
# and looks to see if genes annotated with peaks are significantly enriched in the differentially
# exchanged set
# results go in the same directory as the mapped.peaks.filename
# we're calling sig diff anything with a p value of < 0.01

options(stringsToFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)


xpn.filename <- args[1]
mapped.peaks.filename <- args[2]


#for testing
xpn.filename <- "/mnt/astro_xpn/limma_rd.csv"
mapped.peaks.filename <- "/mnt/data/Astro_chip_xpn_merge.csv"

#xpn.filename<-"/mnt/rory_data/ESC/expression/limma_rd.csv"
#mapped.peaks.filename<-"/mnt/data/ESC_chip_xpn_merge.csv"

#xpn.filename<-"/mnt/rory_data/NS5/expression/limma_rd.csv"
#mapped.peaks.filename<-"/mnt/data/NSC_chip_xpn_merge.csv"




xpn <- read.csv(xpn.filename)
nearest <- read.csv(mapped.peaks.filename)

geneid <- gsub(".*\\((.+)\\)", "\\1",xpn[,"Proportion_Ensembl_transcripts"])
xpn <- cbind(xpn, geneid)

# remove anything we can't map to Ensembl as we have no information about its proximity to peaks
rm <- which(xpn[,"geneid"]=="")
if(length(rm)!=0) xpn <- xpn[-1*rm,]


# Where a probe maps across an exon-exon boundary, a split genome position is given by the ReMOAT
# annotation and we have put these in 2 separate entries in the table.
# Just remove the duplicated entries for counting purposes.
xpn <- xpn[!duplicated(xpn[,"IlluminaID"]),]

# duplicates in the nearest file are because of multiple peaks mapping to the same nearest gene
# the question we're asking here is "are the \de set of genes enriched in genes that are nearest to a peak",
# so we only want to count each gene once, so just delete dups for counting here too
nearest <- nearest[!duplicated(nearest[,"IlluminaID"]),]

not.nearest <- xpn[which(!xpn[,"IlluminaID"] %in% nearest[,"IlluminaID"]),]

n.nearest <- nrow(nearest)
n.not.nearest <- nrow(not.nearest)


cut.off <- 0.01
n.nearest.de <- length(which(nearest[,"P.Value"] <= cut.off))
n.nearest.not.de <- n.nearest-n.nearest.de

n.not.nearest.de <- length(which(not.nearest[,"P.Value"] <= cut.off))
n.not.nearest.not.de <- n.not.nearest-n.not.nearest.de


tab <- cbind(c(n.nearest.de, n.not.nearest.de),
      c(n.nearest.not.de, n.not.nearest.not.de))


cs <- chisq.test(tab)

outfile <- sub("\\.csv","_chisq.txt",mapped.peaks.filename)
write(cs$method, outfile)
write(cs$statistic, outfile, append=T)
write(cs$p.value, outfile, append=T)

# comparative boxplots of fold changes:
outfile <- sub("\\.csv","_fc_boxplot.png",mapped.peaks.filename)
bmp(outfile)
boxplot(nearest[,"logFC"], not.nearest[,"logFC"], names=c("with peak","without peak"), main="log2 FC")
dev.off()


