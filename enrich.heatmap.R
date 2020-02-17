#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library('gplots')
library('RColorBrewer')
library('zoo')
map2color<-function(x,pal,limits=NULL){
	if(is.null(limits)) limits=range(x)
	pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

matr <- 1:2001
enrich.files <- args
enrich <- read.table(enrich.files[1], sep = "\t", header = F)
enrich.order <- order(rowMeans(enrich[,503:1503], na.rm = T), decreasing = T)
enname <- paste(enrich[,1], enrich[,2], sep = ":")
nameorder <- enname[enrich.order]
for (i in 1:length(enrich.files)) {
	enrich <- read.table(enrich.files[i], sep = "\t", header = F)
	names <- paste(enrich[,1], enrich[,2], sep = ":")
	ordername <- unlist(sapply(nameorder, function(x){grep(x, names)}))
	hmatrix <- as.matrix(enrich[, matr + 2])
	hmatrix[hmatrix > 8] <- 8
	hmatrix[is.na(hmatrix)] <- 0
	
	hmatrix <- hmatrix[ordername,]
	png(filename = paste(enrich.files[i], "heatmap.png", sep = "."), width = 632, height = 1308)
	par(mar=c(0,0,0,0))
	hmcol <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))
	plot(0,0, ylim = c(-nrow(hmatrix), 0 ), xlim = range(matr), type = "n", xaxs = "i", yaxs = "i")
	sapply(1:nrow(hmatrix), function(x){points(1:ncol(hmatrix), rep(-x, ncol(hmatrix)), pch = ".", col = map2color(hmatrix[x,], hmcol, limit = c(-5,5)))})
	dev.off()
}
