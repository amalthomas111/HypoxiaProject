#libraries
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript chipcounts.R inputfile outputname ", call.=FALSE)
}
suppressPackageStartupMessages({
  library(gplots)
  library(RColorBrewer)
  #library(MASS)

})
inputfile = args[1]
outputname = args[2]
#inputfile="07_countspromter.tsv"
#outputname="t"
df = read.table(inputfile,sep = "\t",header = T,row.names = 1)
head(df)
nrow(df)
df.filt <- df[rowSums(df)>=10,]
head(df.filt)
nrow(df.filt)
write.table(df.filt,file = paste0(outputname,"_nonzerorows10.tsv"),sep = "\t",quote = F)
df.log = log2(df.filt+1)
corr.matrix <- cor(do.call("cbind",df.log))
corr.matrix
write.table(corr.matrix, file=paste0(outputname,"_correlation.tsv"),sep="\t")
colfunc <- colorRampPalette(c("grey98",brewer.pal(6,"Blues")))
#heatmap
pdf(file = paste0(outputname,"_heatmap.pdf"))
heatmap.2(corr.matrix, col=colfunc(10),margins=c(15,15),symm=F,
          trace="none",distfun=function (x) as.dist(1-x),hclust=function(x) hclust(x,method="ward.D"))
dev.off()
subread_counts<- df.filt[,1:12]

X<- subread_counts
sv = svd(t(X))
U = sv$u
V = sv$v
D = sv$d

Z<- t(X) %*% V

cols = as.numeric(as.factor(colnames(X)))
pdf(file = paste0(outputname,"_pca.pdf"))
plot(Z[,1], Z[,2], type ="n", xlab="PC1", ylab="PC2")
text(Z[,1], Z[,2], colnames(X), col=cols)
dev.off()

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- dist(t(df.filt)) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf(file=paste0(outputname,"_MDS.pdf"))
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="Metric	MDS",	type="n")
text(x, y, labels = colnames(df.filt), cex=.7)
dev.off()

# d <- dist(t(df.filt)) # euclidean distances between the rows
# fit <- isoMDS(d, k=2) # k is the number of dim
# x <- fit$points[,1]
# y <- fit$points[,2]
# pdf(file=paste0(outputname,"_nonmetricMDS.pdf"))
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
#      main="Nonmetric MDS", type="n")
# text(x, y, labels = colnames(df.filt), cex=.7)
# dev.off()

