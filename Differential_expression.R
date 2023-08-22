library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)

########## filter undetected / lowly expressed genes.
x <- as.factor(rep(c("HC", "Learning"), each=4))
names(x) <- colnames(tximport.counts_learning)
filter <- apply(tximport.counts_learning, 1, function(x) length(x[which(x>10)])>5)
filtered <- as.matrix(tximport.counts_learning)[filter,]

########## EDASeq normalization.
uq <- betweenLaneNormalization(filtered, which="upper")
plotRLE(uq, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq, cex=1, cex.axis=1, cex.lab=1)

groups <- matrix(data=c(1:4, 5:8), nrow=2, byrow=TRUE)

########## edgeR with RUVs normalization.
controls <- rownames(uq)
s <- RUVs(uq, controls, k=1, groups)
plotRLE(s$normalizedCounts, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(s$normalizedCounts, cex=1, cex.axis=1, cex.lab=1)

#designRUV
design <- model.matrix(~x + s$W)
design

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
topRsFC <- topTags(lrt, n=Inf)$table
write.table(topRsFC, file="DGEanalysis_RUV.csv", sep = ",")
