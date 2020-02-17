## usage: Rscript qqenrich.R ref.chip ref.input sample.chip sample.input spike.chip spike.input

### modify these parameters accordingly
## chromosome lists of your sample
achrs <- c("MullerB", "MullerE") #autosomes
chrs <- c("Contig_Y1", "Contig_Y2", "MullerA", "MullerAD","MullerB", "MullerC", "MullerE", "MullerF") # all chromosomes


## chromosome lists of your spike in
spk.achrs <- c("2L", "2R", "3L", "3R") # autosomes
spk.chrs <- c("X", "2L", "2R", "3L", "3R", "4") # all chromosomes

## number of window to bin from input for the output. For example if your input window size is 1kb, and you want 10kb windows, the number will be 10
win <- 25
win.spk <- win
###

# miranda chip-seq with melanogaster spike in
options(scipen=999)
library(zoo)

arg = commandArgs(trailingOnly=TRUE)

ref_chp <- arg[1] # reference chip sample
ref_inp <- arg[2] # reference input sample
sam_chp <- arg[3] # actual chip sample
sam_inp <- arg[4] # actual input sample
spk_chp <- arg[5] # spikein chip sample
spk_inp <- arg[6] # spikein input sample

chp_name <- tools::file_path_sans_ext(basename(sam_chp))

pdf(file = paste(c(chp_name, "pdf"), collapse= "."), width = 11.5, height = 8)
par(mar = c(4.5,4.5,4.5,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), widths=c(1,2), heights=c(1,1))


chp <- read.table(ref_chp, sep = "\t", header = F)
inp <- read.table(ref_inp, sep = "\t", header = F)
inps <- c()
chps <- c()
for (c in spk.achrs) {
	inp.t <- inp[inp[,1] == c, ]
	chp.t <- chp[chp[,1] == c, ]
	inp.v <- rollapply(inp.t[,3], win.spk, by = win.spk, mean, partial = T, align = "left")
	chp.v <- rollapply(chp.t[,3], win.spk, by = win.spk, mean, partial = T, align = "left")
	inps <- c(inps, inp.v)
	chps <- c(chps, chp.v)
}
std_inps_med <- median(inps)
std_chps_med <- median(chps)

inps <- c()
chps <- c()
for (c in spk.chrs) {
	inp.t <- inp[inp[,1] == c, ]
	chp.t <- chp[chp[,1] == c, ]
	inp.v <- rollapply(inp.t[,3], win.spk, by = win.spk, mean, partial = T, align = "left")
	chp.v <- rollapply(chp.t[,3], win.spk, by = win.spk, mean, partial = T, align = "left")
	inps <- c(inps, inp.v)
	chps <- c(chps, chp.v)
}
std_en <- (chps/std_chps_med+0.01)/(inps/std_inps_med+0.01)
qr <- quantile(log2(std_en), probs= seq(0,1, by = 0.001))
names(qr) <- as.numeric(sub("%", "", names(qr)))/100


## analyze spike in samples
chp <- read.table(spk_chp, sep = "\t", header = F)
inp <- read.table(spk_inp, sep = "\t", header = F)
chps <- c()
inps <- c()
for (c in spk.achrs) {
	inp.t <- inp[inp[,1] == c, ]
	chp.t <- chp[chp[,1] == c, ]
	inp.v <- rollapply(inp.t[,3], win *2, by = win *2, mean, partial = T, align = "left")
	chp.v <- rollapply(chp.t[,3], win *2, by = win *2, mean, partial = T, align = "left")
	inps <- c(inps, inp.v)
	chps <- c(chps, chp.v)
}
inps_med <- median(inps)
chps_med <- median(chps)

chps <- c()
inps <- c()
for (c in spk.chrs) {
	inp.t <- inp[inp[,1] == c, ]
	chp.t <- chp[chp[,1] == c, ]
	inp.v <- rollapply(inp.t[,3], win, by = win, mean, partial = T, align = "left")
	chp.v <- rollapply(chp.t[,3], win, by = win, mean, partial = T, align = "left")
	inps <- c(inps, inp.v)
	chps <- c(chps, chp.v)
}
lci <- log2((chps/chps_med+0.01)/(inps/inps_med+0.01))
q1 <- quantile(lci, probs= seq(0,1, by = 0.001))
names(q1) <- as.numeric(sub("%", "", names(q1)))/100
qq <- q1-qr # qr is the quantiles in the reference chip

## linear model fit

names(qq) <- round(q1, digits = 4)
l2fc <- as.numeric(names(qq))
lmycept <- sapply(6:(length(qq)-5), function(x) {lm(as.numeric(qq[(x-2):(x+2)]) ~ l2fc[(x-2):(x+2)])$coefficients[1] })
lmslope <- sapply(6:(length(qq)-5), function(x) {lm(as.numeric(qq[(x-2):(x+2)]) ~ l2fc[(x-2):(x+2)])$coefficients[2] })

l2fc <- as.numeric(names(qq)[6:(length(qq)-5)])
lmycept <- c(lm(as.numeric(qq[1:5]) ~ as.numeric(names(qq[1:5])))$coefficients[1], lmycept,
						 lm(as.numeric(qq[(length(qq)-4):length(qq)]) ~ as.numeric(names(qq[(length(qq)-4):length(qq)])))$coefficients[1])
lmslope <- c(lm(as.numeric(qq[1:5]) ~ as.numeric(names(qq[1:5])))$coefficients[2], lmslope,
						 lm(as.numeric(qq[(length(qq)-4):length(qq)]) ~ as.numeric(names(qq[(length(qq)-4):length(qq)])))$coefficients[2])
lmycept <- c(lmycept[1], lmycept)
lmslope <- c(lmslope[1], lmslope)
l2fc <- c(-100, as.numeric(qq[5]), l2fc, 100)
qqdex <- sapply(lci, function(x) {tail(which(l2fc < x), n = 1) + 1})
tr <- as.numeric(lci*lmslope[qqdex] + lmycept[qqdex])

plot(names(qq), qq, pch = 20, las = 1, cex.axis = 0.8, cex.axis = 0.8,
		 ylab = "adjustment", xlab = "spike in log2(enrichment)", main = "QQ adjustment")
sapply(1:length(lmycept), function(x){abline(a=lmycept[x], b=lmslope[x], col = rgb(0,0,1,0.05))})

lci.norm <- lci - tr
plot(lci, pch = 20, cex = 0.4, ylim = c(-3,6), cex.axis = 0.8, cex.lab = 0.8, col = rgb(0,0,1,0.5), las = 1, ylab = "log2(enrichment)",
		 main = paste(c(basename(spk_chp), "spike enrichment"), collapse = " "))
points(lci.norm, pch = 20, cex = 0.4, col = rgb(1,0,0,0.5))


## analyze sample chips
chp <- read.table(sam_chp, sep = "\t", header = F)
inp <- read.table(sam_inp, sep = "\t", header = F)

inps <- c()
chps <- c()

for (c in achrs) {
	inp.t <- inp[inp[,1] == c, ]
	chp.t <- chp[chp[,1] == c, ]
	inp.v <- rollapply(inp.t[,3], win*2, by = win*2, mean, partial = T, align = "left")
	chp.v <- rollapply(chp.t[,3], win*2, by = win*2, mean, partial = T, align = "left")
	inps <- c(inps, inp.v)
	chps <- c(chps, chp.v)
}
inps_med <- median(inps)
chps_med <- median(chps)

inps <- c()
chps <- c()
contigs <- c(0)
chr.v <- c()
start.v <- c()
for (c in chrs) {
	inp.t <- inp[inp[,1] == c, ]
	chp.t <- chp[chp[,1] == c, ]
	inp.v <- rollapply(inp.t[,3], win, by = win, mean, partial = T, align = "left")
	chp.v <- rollapply(chp.t[,3], win, by = win, mean, partial = T, align = "left")
	inps <- c(inps, inp.v)
	chps <- c(chps, chp.v)
	contigs <- c(contigs, contigs[length(contigs)] + length(inp.v))
	chr.v <- c(chr.v, rep(c, length(chp.v)))
	start.v <- c(start.v, seq(1, by = win*1000, length.out = length(chp.v)))
}
zero <- which(chps == 0 & inps == 0)
lci <- log2((chps/chps_med+0.01)/(inps/inps_med+0.01))

qqdex <- sapply(lci, function(x) {tail(which(l2fc < x), n = 1) + 1})
tr <- as.numeric(lci*lmslope[qqdex] + lmycept[qqdex])

plot(lci, pch = 20, cex = 0.4, ylim = c(-2,5), cex.axis = 0.8, cex.lab = 0.8, col = rgb(0,0,1,0.5),
		 main = paste(c(basename(sam_chp), "spike enrichment"), collapse = " "), xaxt='n',
		 xlab = "Chr", ylab="log2(enrichment)", las = 1)
points(lci-tr, pch = 20, cex = 0.4, col = rgb(1,0,0,0.5))
abline(v=contigs, lty = 2)
axis(1, at = contigs[-length(contigs)] + (contigs[-1] - contigs[-length(contigs)])/2, labels = chrs)

write.table(file = paste(c(chp_name, "enrich"), collapse= "."), matrix(c(chr.v, start.v, lci-tr), byrow = F, ncol = 3), col.names = F, row.names = F, sep = "\t", quote = F)

dev.off()
