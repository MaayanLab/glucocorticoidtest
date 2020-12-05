library("limma")
library("voom")
library("edgeR")
library("preprocessCore")
library("Rtsne")

setwd("C:/Users/Rubberchicken/Documents/GitHub/glucocorticoidtest")


counts = read.table("data/GSE120783counts.txt", sep="\t", stringsAsFactors=F)
cols = unlist(counts[1,-c(1,2)])
ens = unlist(counts[-1,1])
genes = unlist(counts[-1,2])

counts = counts[-1, -c(1,2)]
counts = data.matrix(counts)

qexp = normalize.quantiles(log2(counts+1))
rownames(qexp) = genes
colnames(qexp) = cols

rownames(counts) = genes
colnames(counts) = cols

# fix weird ordering
cols2 = cols
cols2[34] = cols[33]
cols2[33] = cols[34]
qexp = qexp[,cols2]
counts = counts[,cols2]

info = read.table("data/sampleinfo.tsv", skip=1, sep="\t", stringsAsFactors=F)
rownames(info) = info[,2]
info = info[,-1]
info = info[colnames(qexp),]

colo = rep(1,34)
colo[which(info[,5] == "batch2")] = 2
colo[which(info[,5] == "batch3")] = 3

pcc = rep(20,24)
pcc[which(info[,2] == "gc")] = 17

pp = prcomp(t(qexp))
plot(pp$x[,1], pp$x[,2], col=colo, pch=pcc, xlab="PC1", ylab="PC2", cex=1.3)
legend("topleft", legend=c("Batch 1", "Batch 2", "Batch 3", "control", "GC"), col=c(1,2,3,1,1), pch=c(20,20,20,20,17))

# since the samples are matched there is a lot more power in a matched analysis

controls = qexp[,which(info[,2] == "control")]
gc = qexp[,which(info[,2] == "gc")]

tts = list()
ttp = list()

for(i in 1:nrow(qexp)){
    tv = t.test(controls[i,], gc[i,], paired = T)
    tts[[length(tts)+1]] = tv$statistic
    ttp[[length(ttp)+1]] = tv$p.value
}

tts = unlist(tts)
tts[is.na(tts)] = 0
ttp = unlist(ttp)
ttp[is.na(ttp)] = 1
apv = p.adjust(ttp, method="fdr")

ww1 = which(tts > 0 & apv < 0.1)
ww2 = which(tts < 0 & apv < 0.1)

oo = rev(order(tts))

genes[oo][1:10]
genes[rev(oo)][1:10]

plot(tts[oo])

# the results point to very similar results

fbatch = as.factor(info[,5])
ftreatment = as.factor(info[,2])
fgender = as.factor(info[,4])
fethnicity = as.factor(info[,7])
findividual = as.factor(info[,6])

# this model is not matched by subject. This is possible somehow too and will make the results better
design <- model.matrix(~ftreatment)

dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

top = topTable(fit, n=Inf)

write.table(top, file="limma_unmatched_gc.tsv", sep="\t", quotes=F)

