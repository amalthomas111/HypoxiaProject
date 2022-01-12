#libraries
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript chipcounts.R inputfile outputname ", call.=FALSE)
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
#inputfile = "allsamples_jaccard.tsv"
#outputname = "allsamples"
df = read.table(inputfile)#,sep = "\t",header = T,row.names = 1)
head(df)
type = c("Brx07","Brx07","Brx07","Brx07","Brx07","Brx07","Brx07","Brx07","Brx07","Brx07","Brx07","Brx07",
         "Brx42","Brx42","Brx42","Brx42","Brx42",
         "Brx50","Brx50","Brx50","Brx50","Brx50",
"Brx68","Brx68","Brx68","Brx68","Brx68"
)

metstat = c("parent","parent","Kidneymets","Lungmets","Lungmets","Lungmets","Lungmets",
            "Lyndmets","Ovarymets","Ovarymets","Ovarymets","Ovarymets",
             "Brainmets","Brainmets","parent","parent","Ovarymets",
             "Bonemets","Brainmets","Brainmets","parent","parent",
            "Bonemets","Bonemets","parent","parent","Lungmets")
celltype <- palette(brewer.pal(8, "Dark2"))[as.fumeric(type)]
parent_mets = palette(brewer.pal(8, "Dark2"))[as.fumeric(metstat)]
#head(cbind(colnames(df),cols),n=20)

df_new= as.matrix(df)
diag(df_new)=NA

pdf(file=paste0(outputname,"_jaccard.pdf"))
heatmap.2(as.matrix(df_new), trace="none", ColSideColors=celltype, col=hmcol,margins=c(15,15),symm=F)
dev.off()
#heatmap.2(as.matrix(df), trace="none", ColSideColors=c(cols,cols1), col=hmcol,margins=c(15,15),symm=F)

library("gplots")
library("devtools")
source("/home/master/git/heatmap3/heatmap.3.R")


clab = cbind(celltype,parent_mets)
pdf(file=paste0(outputname,"_heatmap3.pdf"))
heatmap.3(as.matrix(df_new),Rowv=TRUE, Colv=TRUE, ColSideColors=clab,col=hmcol,margins=c(6,12),symm=F
          )
dev.off()
