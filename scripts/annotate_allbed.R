suppressPackageStartupMessages({
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(clusterProfiler)
})
files <- list.files(path=getwd(), pattern="*.txt$", full.names=T, recursive=FALSE)
files
names(files) = gsub("_peaks.bed","",basename(files))
names(files) = gsub(".R0","",names(files))
names(files) = gsub("_peaks.narrowPeak.blacklistcleared","",names(files))


mainDir=getwd()
subDir="annotation"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
#taglist <-  lapply(files, getTagMatrix, windows=promoter)
#plotAvgProf(taglist, xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")

count=1
for (i in files) {
  peakAnno <- annotatePeak(i, tssRegion=c(-3000, 3000), 
                           TxDb=txdb, annoDb="org.Hs.eg.db")
  output=paste0("annotation/annotate_",names(files[count]),".tsv")
  count=count+1
  write.table(as.data.frame(peakAnno),file = output,sep="\t",quote = F, row.names = F)
}
