library(rtracklayer)
library(diffHic)
gr_obj <- import("/home/gdstantonlab/mxw010/Data/SPICE-C/data/bins.bed")
temp <- read.table("/home/gdstantonlab/mxw010/Data/SPICE-C/data/chrom.sizes")
#remove y chr
seqlengths(gr_obj) <- temp[,2]
hs.param <- pairParam(gr_obj)
library(edgeR)

# working directory:
# /home/gdstantonlab/mxw010/Data/SPICE-C/data/diffHiC

input <- c('RH4_rep1.h5','RH4_rep2.h5', 'RH30_rep1.h5', 'RH30_rep2.h5', 'CTR_rep1.h5', 'CTR_rep2.h5', 'RD_rep1.h5','RD_rep2.h5')

finder <- domainDirections(input, hs.param, width=50000, span=20)

all.counts <- cbind(assay(finder, "up"), assay(finder, "down"))
totals <- totalCounts(input, hs.param)


ydom <- DGEList(all.counts, lib.size=rep(totals, 2))


Cell <- factor(rep(c("RH4", 'RH30', "CTR", "RD"),each=2))
Sample <- factor(seq_along(input))
Cell <- rep(Cell, 2)
Sample <- rep(Sample, 2)
Condition=factor(rep(c('FP','FN'), each=4))
Condition=rep(Condition,2)
Direction <- rep(c("Up", "Down"), each=length(input))
design <- model.matrix(~0 + Cell + Direction:Condition)

design <- design[,!grepl("DirectionDown", colnames(design))]
colnames(design) <- sub("DirectionUp:", "", colnames(design))


ab <- aveLogCPM(ydom)
keep <- ab > 0
ydom <- ydom[keep,]
summary(keep)

cur.regions <- rowRanges(finder)[keep,]

ydom <- estimateDisp(ydom, design)
plotBCV(ydom)

fitdom <- glmQLFit(ydom, design, robust=TRUE)
plotQLDisp(fitdom)

con <- makeContrasts(ConditionFP - ConditionFN, levels=design)
resdom <- glmQLFTest(fitdom, contrast=con)
topTags(resdom)

#LRT

resdom <- glmLRT(fitdom)

output <- data.frame(as.data.frame(cur.regions)[,1:3],
FP=fitdom$coefficients[,"ConditionFP"]/log(2),
FN=fitdom$coefficients[,"ConditionFN"]/log(2),resdom$table)
output$FDR <- p.adjust(resdom$table$PValue, method="BH")
o <- order(output$PValue)
output <- output[o,]
length(which(abs(output$logFC) > log(1.5,base=2) & output$FDR < 0.01))



is.sig <- output$FDR <= 0.05 & abs(output$logFC) > log(1.5,base=2)
change.type <- abs(output$FP) > abs(output$FN)

# summary(change.type[is.sig])
#    Mode   FALSE    TRUE
# logical     251     376

length(which(output$FDR < 0.01))

xchr <- as.character(seqnames(region))
x.min <- max(1L, start(region))
x.max <- min(seqlengths(fragments)[[xchr]], end(region))

write.table(output[output$FDR<0.01,],'differential_TAD_50kb.bed',sep="\t", quote=F, row.names=F)