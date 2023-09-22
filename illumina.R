setwd(R"(C:\Dessertation\metabric_analysis)")

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
sw(library(dplyr))
sw(library(ggpubr))

list.files()

df = read.table('illumina_microarraywithoutentrenz.txt', sep='\t',header = 1)

row_names <- df[, 1]
row.names(df) <- make.unique(row_names)
sub_df <- df[, -c(1)]
targets = read.table('cleaned_clinical_patient.txt', sep='\t',header = 1)
colnames(sub_df)
head(targets)
df_sub=read.table('cleaned_dataagin.txt', sep='\t',header = 1)
head(df_sub)
row.names(df_sub)<-row.names(sub_df)
write.table(df_sub, "NewlysortedDF.txt", sep='\t')
conditions<- paste(targets$Groups,targets$INFERRED_MENOPAUSAL_STATE)
head(conditions)
conditions <- factor(conditions, levels=unique(conditions))
design <- model.matrix(~0+ conditions)
colnames(design) <- levels(conditions)
unique(conditions)
Newname=make.names(colnames(design))
colnames(design)<- Newname
#fit the design to the expression data
fit = lmFit(df_sub, design)
# Use the valid names in the contrast matrix
contrast.matrix <- makeContrasts(HRstatus=HR..Pre-HR..Post,
                                 HER2status=HER2..Pre-HER2..Post,
                                 TNBCstatus=TNBC.Pre-TNBC.Post,
                                 levels = design)
#Fit contrast matrix
fit.cont = contrasts.fit(fit, contrast.matrix)
#Fit ebayes
fit.cont = eBayes(fit.cont)


res=topTable(fit.cont, coef = "HRstatus", adjust.method = "BH", number = nrow(df_sub), sort.by = 'p')
head(res)
write.table(res, "HRstatusDEGS.txt", sep='\t')
p=EnhancedVolcano(res,
                  lab = row.names(res),
                  x = 'logFC',
                  y = 'P.Value',
                  title = 'HR+', subtitle='',
                  pCutoff = 0.05,
                  legendLabels = c("NS", expression(Log[2] ~ FC), "pvalue", expression(pvalue ~ and ~ log[2] ~ FC)),
                  FCcutoff = 0.5)
png("HR+StatusDEGS.volcano.png", units="in", width=10, height=10, res=900)
p
dev.off()

resher=topTable(fit.cont, coef = "HER2status", adjust.method = "BH", number = nrow(df_sub), sort.by = 'p')
head(resher)
write.table(resher, "HER2statusDEGS.txt", sep='\t')
q=EnhancedVolcano(resher,
                  lab = row.names(resher),
                  x = 'logFC',
                  y = 'P.Value',
                  title = 'HR+', subtitle='',
                  pCutoff = 0.05,
                  legendLabels = c("NS", expression(Log[2] ~ FC), "pvalue", expression(pvalue ~ and ~ log[2] ~ FC)),
                  FCcutoff = 0.5)
png("HER2StatusDEGS.volcano.png", units="in", width=10, height=10, res=900)
q
dev.off()

restnbc=topTable(fit.cont, coef = "TNBCstatus", adjust.method = "BH", number = nrow(df_sub), sort.by = 'p')
head(restnbc)
write.table(restnbc, "TNBCstatusDEGS.txt", sep='\t')
k=EnhancedVolcano(restnbc,
                  lab = row.names(restnbc),
                  x = 'logFC',
                  y = 'P.Value',
                  title = 'HR+', subtitle='',
                  pCutoff = 0.05,
                  legendLabels = c("NS", expression(Log[2] ~ FC), "pvalue", expression(pvalue ~ and ~ log[2] ~ FC)),
                  FCcutoff = 0.5)
png("TNBCStatusDEGS.volcano.png", units="in", width=10, height=10, res=900)
k
dev.off()




customgenes = c("PTK2",	"PTK2B",	"PTK7",	"PTK6",	"TP53",	"BRCA1",	"BRCA2",	"ERCC1",	"ERCC2",	"ERC1",	"ERC2"	,"EGFR",	"STAT1",	"STAT2",	"AKT1",	"AKT2",	"AKT3",	"SRC")

group_df = data.frame(Groups=rep(c("HERPost", "HRPost",
                                   "TNBCPost", "HERPre",
                                   "HRPre", "TNBCPre"), c(169,1180,207,78,233,113)))
rownames(group_df) <- colnames(df_sub)

a=pheatmap(df_sub[rownames(df_sub) %in% customgenes,], scale="row",
           cluster_cols=FALSE,
           annotation_col=group_df,
            show_colnames=F,
           colorRampPalette(c("navy", "white", "red"))(75),gaps_col=cumsum(c(169,1180,207,78,233,113)))
png(file="heatmap_selected_genes.png", width=14, height=14, units="in", res=300)
a
dev.off()











