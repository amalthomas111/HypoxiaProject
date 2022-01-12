args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript GC.R <counts> <gc table>", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input file not found!\
Usage: Rscript peakannot.R <bedfile1> <bedfile2>", call.=FALSE)
}
suppressPackageStartupMessages({
library(geneplotter)
})
counts=args[1]
gc = args[2]
#counts="hg19_genome.ctcf.tsv"
#gc="hg19.genome.10kb.gc" 

df_count=read.table(counts,sep="\t",header = T,row.names = 1)#,col.names = 1,row.names = 1)
head(df_count)


df_gc=read.table(gc,sep="\t",row.names = 1)
names(df_gc)="gc"
head(df_gc)
count=1
for(i in 1:ncol(df_count)){
  png(filename = paste0(count,".smooth.",colnames(df_count)[i],".png"))
  smoothScatter(df_gc$gc,df_count[,i],xlim=c(0.2,max(df_gc$gc)),ylim=c(1,300),main = colnames(df_count)[i],pch = '.',ylab="Read Counts",xlab="GC content")
  dev.off()
  count=count+1
  }

for(i in 1:ncol(df_count)){
  png(filename = paste0(count,".dot.",colnames(df_count)[i],".png"))
  plot(df_gc$gc,df_count[,i],xlim=c(0.2,max(df_gc$gc)),ylim=c(1,300),main = colnames(df_count)[i],pch = '.',ylab="Read Counts",xlab="GC content")
  dev.off()
  count=count+1
}


####################
#log plot
####################


count=1
for(i in 1:ncol(df_count)){
  png(filename = paste0(count,".log.smooth.",colnames(df_count)[i],".png"))
  smoothScatter(df_gc$gc,log2(df_count[,i]+1),xlim=c(0.2,max(df_gc$gc)),ylim=c(1,15),main = colnames(df_count)[i],pch = '.',ylab=expression(paste("log"[2]," Read counts")),xlab="GC content")
  dev.off()
  count=count+1
}

for(i in 1:ncol(df_count)){
  png(filename = paste0(count,".log.dot.",colnames(df_count)[i],".png"))
  plot(df_gc$gc,log2(df_count[,i]+1),xlim=c(0.2,max(df_gc$gc)),ylim=c(1,15),main = colnames(df_count)[i],pch = '.',ylab=expression(paste("log"[2]," Read counts")),xlab="GC content")
  dev.off()
  count=count+1
}

