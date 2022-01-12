#libraries
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript chipcounts.R inputfile outputname", call.=FALSE)
}

suppressPackageStartupMessages({
library(genefilter)
library(gplots)
library(rafalib)
library(RColorBrewer)
})
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

inputfile = args[1]
outputname = args[2]
df = read.table(inputfile)#,sep = "\t",header = T,row.names = 1)
head(df)
df_new= as.matrix(df)
diag(df_new)=NA
#type = c("parent","parent","parent","parent","parent","parent",
#                  "Kidney-mets","Lung-mets","Lung-mets","Lung-mets", "Lung-mets","Lynd-mets",
#                  "Ovary-mets","Ovary-mets","Ovary-mets","Ovary-mets")
type = c("68.1","68.1","68.2","68.2","68.3",
"68.3","68.3","68.4","68.4","mut.1","mut.1","mut.2","mut.2","mut.3","mut.3","wt.1","wt.1","wt.2","wt.2")
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(type)]
head(cbind(colnames(df_new),cols))
pdf(file=paste0(outputname,"_jaccard.pdf"))
heatmap.2(as.matrix(df_new), trace="none", na.rm=TRUE,labCol=FALSE, RowSideColors=cols,col=hmcol,margins=c(6,12),symm=F,,Rowv=TRUE, Colv=TRUE)
dev.off()
