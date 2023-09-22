setwd(R"(C:\Dessertation\metabric_analysis)")
list.files()
rm(list=ls())
options(warn=-1)
sw=suppressPackageStartupMessages
sw(library(limma))
sw(library(pheatmap))
sw(library(EnhancedVolcano))
sw(library(dplyr))
sw(library(tidyverse))
sw(library(illuminaHumanv3.db))
sw(library(RColorBrewer))
#Enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GOplot)
DEGHR=read.table('DEGSselectedinHR+.txt', sep='\t',header = 1)
as.data.frame(DEGHR)
GO_results<-read.table('HR+Enrichment3top.txt', sep='\t',fill=TRUE,header = 1)
head(GO_results)
GO_sub<-subset(GO_results,GO_results$adj_pval  <= 0.05 )

head(DEGHR)
dim(DEGHR)
dim(GO_results)

circ <- circle_dat(GO_results, DEGHR)

# Generate the plot with custom colors
plots <- GOBar(subset(circ,category == ' Biological Process'))


ggsave("HR+EnrichmentanalysisBarplot.png", plots, width = 49, height =40,dpi = 500)

g=GOBubble(circ, title = 'Bubble plot', colour = c("red", "green", "blue", "yellow", "purple", "orange", "black", "white", "gray", "pink", "brown", "cyan", "magenta", 'darkred', 'gold'), display = ' Biological Process', labels = 15)  

ggsave("HR+EnrichmentanalysisBubbleplot.png", g, width = 49, height =40,dpi = 500)

 



q=GOBubble(circ,labels = 3)

ggsave("HR+EnrichmentanalysisBubbleplot3only.png", q, width = 49, height =40,dpi = 500)

L=GOBubble(circ, title = 'Bubble plot', colour = c("pink", "brown", "cyan"), display = 'multiple', labels = 3) 

ggsave("HR+EnrichmentanalysisBubbleplot3onlyseperated.png", L, width = 49, height =40,dpi = 500)

c=GOCircle(circ)
ggsave("HR+EnrichmentanalysisCircle3top.png", c, width = 15, height =15,dpi = 500)


DEGHER2=read.table('DEGSselectedinHER2+.txt', sep='\t',header = 1)
GO_results2<-read.table('HER2+enrichmentTOP2.txt', sep='\t',fill=TRUE,header = 1)
circ2 <- circle_dat(GO_results2, DEGHER2)
plots2 <- GOBar(subset(circ2,category == ' Biological Process'))

ggsave("HER2+EnrichmentanalysisBarplot.png", plots2, width = 20, height =20,dpi = 500)
e=GOBubble(circ2,labels = 2)

ggsave("HER2+EnrichmentanalysisBubbleplot2only.png", e, width =20 , height =20,dpi = 500)

y=GOBubble(circ2, title = 'Bubble plot', colour = c( "red", "green"), display = 'multiple', labels = 2) 

ggsave("HER2+EnrichmentanalysisBubbleplot2onlyseperated.png", y, width = 20, height =20,dpi = 500)
w=GOCircle(circ2)
ggsave("HER2+EnrichmentanalysisCircle2top.png", w, width = 15, height =15,dpi = 500)

DEGTNBC=read.table('DEGSselectedinTNBC.txt', sep='\t',header = 1)
GO_results3<-read.table('TNBCEnrichmentTOP3.txt', sep='\t',fill=TRUE,header = 1)
circ3 <- circle_dat(GO_results3, DEGTNBC)
plots3 <- GOBar(subset(circ3,category == ' Biological Process'))

ggsave("TNBCEnrichmentanalysisBarplot.png", plots3, width = 25, height =25,dpi = 500)
t=GOBubble(circ3,labels = 3)

ggsave("TNBCEnrichmentanalysisBubbleplot3only.png", t, width = 25, height =20,dpi = 500)

q=GOBubble(circ3, title = 'Bubble plot', colour = c( "purple", "orange"), display = 'multiple', labels = 3) 

ggsave("TNBCEnrichmentanalysisBubbleplot3onlyseperated.png", q, width = 15, height =15,dpi = 500)
o=GOCircle(circ3)
ggsave("TNBCEnrichmentanalysisCircle3top.png", o, width = 15, height =15,dpi = 500)
