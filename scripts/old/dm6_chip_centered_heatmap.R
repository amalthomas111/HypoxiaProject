library(EnrichedHeatmap) # for making heatmap
library(rtracklayer)  # for reading bigwig and bed files
library(GenomicRanges)
library(circlize)
library(dplyr)
library(tidyverse)

asi_bed = import("AsiSi.bed",format="BED")
asi_bed.diff = import("Comparison2_phH2av_TRT-vs-UNT_peakintersectn_DESEQ2_up_10kb.0.01.bed",format = "BED")

#read bigwig file

ph.TRT.R1 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/deeptools_bw/phH2Av_Treat_R1.bw")
ph.TRT.R2 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/deeptools_bw/phH2Av_Treat_R2.bw")
h3k9.TRT.R1 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/deeptools_bw/H3K9me3_Treat_R1.bw")
h3k9.TRT.R2 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/deeptools_bw/H3K9me3_Treat_R2.bw")
ph.UNT.R1= import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/deeptools_bw/phH2Av_UNT_R1.bw")
ph.UNT.R2= import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/deeptools_bw/phH2Av_UNT_R2.bw")


ph.TRT.R1 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/my_bw_tracks/phH2Av_Treat_R1.bw")
ph.TRT.R2 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/my_bw_tracks/phH2Av_Treat_R2.bw")
h3k9.TRT.R1 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/my_bw_tracks/H3K9me3_Treat_R1.bw")
h3k9.TRT.R2 = import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/my_bw_tracks/H3K9me3_Treat_R2.bw")
ph.UNT.R1= import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/my_bw_tracks/phH2Av_UNT_R1.bw")
ph.UNT.R2= import.bw("/media/master/hdd/WORK/chipseq/private/DSB_Laeiti/bowtie1/my_bw_tracks/phH2Av_UNT_R2.bw")


cutsite.10Kb = resize(asi_bed, width = 10000, fix = "center")
cutsite.center = resize(cutsite.10Kb,width = 1, fix = "center")

cutsite.10Kb.diff = resize(asi_bed.diff, width = 10000, fix = "center")
cutsite.center.diff = resize(cutsite.10Kb.diff,width = 1, fix = "center")

ph.TRT.R1.mat = normalizeToMatrix(ph.TRT.R1,cutsite.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)
ph.TRT.R2.mat = normalizeToMatrix(ph.TRT.R2 ,cutsite.center, value_column = "score", mean_mode="w0", w=100, extend = 5000)
ph.UNT.R1.mat = normalizeToMatrix(ph.UNT.R1,cutsite.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)
ph.UNT.R2.mat = normalizeToMatrix(ph.UNT.R2,cutsite.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)

quantile(ph.TRT.R1.mat , probs = c(0.005, 0.5,0.95))
quantile(ph.TRT.R2.mat, probs = c(0.005, 0.5,0.95))
quantile(ph.UNT.R1.mat, probs = c(0.005, 0.5,0.95))
quantile(ph.UNT.R2.mat, probs = c(0.05,0.5,0.95))

ph.TRT.R1.diff.mat = normalizeToMatrix(ph.TRT.R1,cutsite.center.diff, value_column = "score",mean_mode="w0", w=100, extend = 5000)
ph.TRT.R2.diff.mat = normalizeToMatrix(ph.TRT.R2 ,cutsite.center.diff, value_column = "score", mean_mode="w0", w=100, extend = 5000)
ph.UNT.R1.diff.mat = normalizeToMatrix(ph.UNT.R1,cutsite.center.diff, value_column = "score",mean_mode="w0", w=100, extend = 5000)
ph.UNT.R2.diff.mat = normalizeToMatrix(ph.UNT.R2,cutsite.center.diff, value_column = "score",mean_mode="w0", w=100, extend = 5000)

quantile(ph.TRT.R1.diff.mat , probs = c(0.005, 0.5,0.95))
quantile(ph.TRT.R2.diff.mat, probs = c(0.005, 0.5,0.95))
quantile(ph.UNT.R1.diff.mat, probs = c(0.005, 0.5,0.95))
quantile(ph.UNT.R2.diff.mat, probs = c(0.05,0.5,0.95))

## from the quantile, I choose the color mapping range
col_fun <- circlize::colorRamp2(c(0, 350), c("white", "red"))
col_fun_ER.nogata.upregion<- circlize::colorRamp2(c(0, 350), c("white", "red"))
col_fun_GATA.upregion<- circlize::colorRamp2(c(0, 350), c("white", "red"))

col_fun.diff <- circlize::colorRamp2(c(0, 15), c("white", "red"))

EnrichedHeatmap(ph.TRT.R1.mat, axis_name_rot = 0, name = "ph.TRT.R1",
                column_title = "ph.TRT.R1", use_raster = TRUE, col = col_fun,
                top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                top_annotation_height = unit(2, "cm")) +
  EnrichedHeatmap(ph.TRT.R2.mat, axis_name_rot = 0, name = "ph.TRT.R2", 
                  column_title = "ph.TRT.R2", use_raster = TRUE, col = col_fun,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                  top_annotation_height = unit(2, "cm")) + 
  EnrichedHeatmap(ph.UNT.R1.mat, axis_name_rot = 0, name = "ph.UNT.R1", 
                  column_title = "ph.UNT.R1", use_raster = TRUE, col = col_fun,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                  top_annotation_height = unit(2, "cm")) +
  EnrichedHeatmap(ph.UNT.R2.mat, axis_name_rot = 0, name = "ph.UNT.R2", 
                  column_title = "ph.UNT.R2", use_raster = TRUE, col = col_fun,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                  top_annotation_height = unit(2, "cm"))


EnrichedHeatmap(ph.TRT.R1.diff.mat, axis_name_rot = 0, name = "ph.TRT.R1",
                column_title = "ph.TRT.R1", use_raster = TRUE, col = col_fun.diff,
                top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                top_annotation_height = unit(2, "cm")) +
  EnrichedHeatmap(ph.TRT.R2.diff.mat, axis_name_rot = 0, name = "ph.TRT.R2", 
                  column_title = "ph.TRT.R2", use_raster = TRUE, col = col_fun.diff,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                  top_annotation_height = unit(2, "cm")) + 
  EnrichedHeatmap(ph.UNT.R1.diff.mat, axis_name_rot = 0, name = "ph.UNT.R1", 
                  column_title = "ph.UNT.R1", use_raster = TRUE, col = col_fun.diff,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                  top_annotation_height = unit(2, "cm")) +
  EnrichedHeatmap(ph.UNT.R2.diff.mat, axis_name_rot = 0, name = "ph.UNT.R2", 
                  column_title = "ph.UNT.R2", use_raster = TRUE, col = col_fun.diff,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()),
                  top_annotation_height = unit(2, "cm"))


ER.withgata.upregion_mean = data.frame(avg = colMeans(ER.withgata.upregion.mat), 
                              CI_lower = apply(ER.withgata.upregion.mat, 2, Hmisc::smean.cl.normal)[2,],
                              CI_upper = apply(ER.withgata.upregion.mat, 2, Hmisc::smean.cl.normal)[3,]) %>% 
  mutate(factor = "ER control", pos = colnames(ER.withgata.upregion.mat)) 

ER.nogata.upregion_mean = data.frame(avg = colMeans(ER.nogata.upregion.mat), 
                                     CI_lower = apply(ER.nogata.upregion.mat, 2, Hmisc::smean.cl.normal)[2,],
                                     CI_upper = apply(ER.nogata.upregion.mat, 2, Hmisc::smean.cl.normal)[3,]) %>%
  mutate(factor = "ER siGATA3", pos = colnames(ER.nogata.upregion.mat))
GATA.upgregion_mean = data.frame(avg = colMeans(GATA.upregion.mat), 
                                 CI_lower = apply(GATA.upregion.mat, 2, Hmisc::smean.cl.normal)[2,],
                                 CI_upper = apply(GATA.upregion.mat, 2, Hmisc::smean.cl.normal)[3,]) %>%
  mutate(factor = "GATA3", pos = colnames(GATA.upregion.mat))


ph.TRT.R1.diff.mat_mean = data.frame(avg = colMeans(ph.TRT.R1.diff.mat), 
                                     CI_lower = apply(ph.TRT.R1.diff.mat, 2, Hmisc::smean.cl.normal)[2,],
                                     CI_upper = apply(ph.TRT.R1.diff.mat, 2, Hmisc::smean.cl.normal)[3,]) %>% 
  mutate(factor = "ph.TRT.R1", pos = colnames(ph.TRT.R1.diff.mat)) 
ph.TRT.R2.diff.mat_mean = data.frame(avg = colMeans(ph.TRT.R2.diff.mat), 
                                     CI_lower = apply(ph.TRT.R2.diff.mat, 2, Hmisc::smean.cl.normal)[2,],
                                     CI_upper = apply(ph.TRT.R2.diff.mat, 2, Hmisc::smean.cl.normal)[3,]) %>% 
  mutate(factor = "ph.TRT.R2", pos = colnames(ph.TRT.R2.diff.mat)) 
ph.UNT.R1.diff.mat_mean = data.frame(avg = colMeans(ph.UNT.R1.diff.mat), 
                                     CI_lower = apply(ph.UNT.R1.diff.mat, 2, Hmisc::smean.cl.normal)[2,],
                                     CI_upper = apply(ph.UNT.R1.diff.mat, 2, Hmisc::smean.cl.normal)[3,]) %>% 
  mutate(factor = "ph.UNT.R1", pos = colnames(ph.UNT.R1.diff.mat)) 
ph.UNT.R2.diff.mat_mean = data.frame(avg = colMeans(ph.UNT.R2.diff.mat), 
                                     CI_lower = apply(ph.UNT.R2.diff.mat, 2, Hmisc::smean.cl.normal)[2,],
                                     CI_upper = apply(ph.UNT.R2.diff.mat, 2, Hmisc::smean.cl.normal)[3,]) %>% 
  mutate(factor = "ph.UNT.R2", pos = colnames(ph.UNT.R2.diff.mat)) 

combine_both = bind_rows(ph.TRT.R1.diff.mat_mean,ph.TRT.R2.diff.mat_mean,ph.UNT.R1.diff.mat_mean,ph.UNT.R2.diff.mat_mean)
combine_both$pos<- factor(combine_both$pos, levels= ph.TRT.R2.diff.mat_mean$pos)

ggplot(combine_both, aes(x = pos,y = avg, group = factor)) + geom_line(aes(color = factor)) + 
  theme_bw(base_size = 14) +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(breaks = c("u1", "d1", "d50"), labels =c ("-5kb", "Peak center", "5kb")) +
  xlab(NULL) + 
  ylab("RPKM")+
  ggtitle("ChIP-seq signal") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "avgplot_diffupregions.0.01.10kb.pdf")
## take some touch up to finalize the figure
ggplot(combine_both, aes(x = pos,y = avg, group = factor)) + geom_line(aes(color = factor)) + 
  geom_ribbon(aes(ymin= CI_lower, ymax=CI_upper), alpha=0.2, col = "#8B7E66") +
  theme_bw(base_size = 14) +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(breaks = c("u1", "d1", "d50"), labels =c ("-5kb", "center", "5kb")) +
  xlab(NULL) + 
  ylab("RPKM")+
  ggtitle("ChIP-seq signal")+ theme(plot.title = element_text(hjust = 0.5))
