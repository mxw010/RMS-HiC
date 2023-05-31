args = commandArgs(trailingOnly=TRUE)
chr <- args[1]
res <- args[2]
window <- unlist(strsplit(args[3:length(args)],","))
cell <- args[4]

#setwd(paste("result_", res,"/",sep=""))
#print(paste("working directory has changed to ", getwd(), sep=""))
library(TopDom)
x <- readHiC(paste("/home/gdstantonlab/mxw010/Data/SPICE-C/", cell,"/TopDom/data_", res, "/",chr,sep=""))
chr=unique(x$bins$chr)
#result <- TopDom(x, window.size=5, outFile = chr)
for (size in window) {
    print(paste("running window size ", size, sep=""))
    TopDom(x, window.size=size, outFile = paste("result_", res,"/", chr, "_", size, sep=""))
}
