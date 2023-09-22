setwd(R"(C:\Dessertation\metabric_analysis\Ana;ysisExtra)")
rm(list=ls())
options(warn=-1)
sw = suppressPackageStartupMessages
sw(library(dplyr))
sw(library(data.table))
sw(library(survival))
sw(library(survminer))
sw(library(limma))
sw(library(ggpubr))
sw(library(ggplot2))
sw(library(imogimap))
sw(library(immunedeconv))
sw(library(corrplot))
df = read.delim("ptk2_mta1_data.txt", sep='\t')
dim(df)
head(df,2)
q1 = quantile(df$PTK2,0.33)
q2 = quantile(df$PTK2,0.66)
df$PTK2Group = with(df,
                    ifelse(PTK2 <= q1, 'Low',
                           ifelse(PTK2 > q1 & PTK2 < q2, 'Moderate',
                                  ifelse(PTK2 >= q2, 'High', 'unknown'))))
tnbc = filter(df, StatusGroup == "TNBC")
hrpos = filter(df, StatusGroup == "HRPOS")
Her2=filter(df,StatusGroup=="HER2Enriched")
dim(tnbc)
dim(hrpos)
dim(Her2)
q3 = quantile(tnbc$PTK2,0.33)
q4 = quantile(tnbc$PTK2,0.66)
tnbc$PTK2Group = with(tnbc,
                      ifelse(PTK2 <= q3, 'Low',
                             ifelse(PTK2 > q3 & PTK2 < q4, 'Moderate',
                                    ifelse(PTK2 >= q4, 'High', 'unknown'))))
q5 = quantile(hrpos$PTK2,0.33)
q6 = quantile(hrpos$PTK2,0.66)
hrpos$PTK2Group = with(hrpos,
                      ifelse(PTK2 <= q5, 'Low',
                             ifelse(PTK2 > q5 & PTK2 < q6, 'Moderate',
                                    ifelse(PTK2 >= q6, 'High', 'unknown'))))
q7 = quantile(Her2$PTK2,0.33)
q8 = quantile(Her2$PTK2,0.66)
Her2$PTK2Group = with(Her2,
                      ifelse(PTK2 <= q7, 'Low',
                             ifelse(PTK2 > q7 & PTK2 < q8, 'Moderate',
                                    ifelse(PTK2 >= q8, 'High', 'unknown'))))
write.table(tnbc, "tnbc_metabric.fak.signature.txt", sep='\t')
write.table(hrpos, "hrpos_metabric.fak.signature.txt", sep='\t')
write.table(Her2, "Her2_metabric.fak.signature.txt", sep='\t')
centr = theme_grey() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
wdata = filter(df, PTK2Group == "High" | PTK2Group == "Low")
wdata= arrange(wdata, PTK2Group)
wdata_M = filter(wdata, OS_MONTHS <= 120)
SAHR = filter(hrpos, PTK2Group == "High" | PTK2Group == "Low")
SAHR = arrange(SAHR, PTK2Group)
SAHR_M = filter(SAHR, OS_MONTHS <= 120)
SAHER = filter(Her2, PTK2Group == "High" | PTK2Group == "Low")
SAHER = arrange(SAHER, PTK2Group)
SAHER_M = filter(SAHER, OS_MONTHS <= 120)
SATNBC = filter(tnbc, PTK2Group == "High" | PTK2Group == "Low")
SATNBC = arrange(SATNBC, PTK2Group)
SATNBC_M = filter(SATNBC, OS_MONTHS <= 120)
#Survival analysis without cutoff
fitdata1 = survfit(Surv(OS_MONTHS, OS_STATUS)~PTK2Group, data=wdata_M)
fitdata2 = survfit(Surv(RFS_MONTHS, RFS_STATUS)~PTK2Group, data=wdata_M)

fitHR1 = survfit(Surv(OS_MONTHS, OS_STATUS)~PTK2Group, data=SAHR_M)
fitHR2 = survfit(Surv(RFS_MONTHS, RFS_STATUS)~PTK2Group, data=SAHR_M)

fitHER1 = survfit(Surv(OS_MONTHS, OS_STATUS)~PTK2Group, data=SAHER_M)
fitHER2 = survfit(Surv(RFS_MONTHS, RFS_STATUS)~PTK2Group, data=SAHER_M)

fitTNBC1 = survfit(Surv(OS_MONTHS, OS_STATUS)~PTK2Group, data=SATNBC_M)
fitTNBC2 = survfit(Surv(RFS_MONTHS, RFS_STATUS)~PTK2Group, data=SATNBC_M)
splots <- list()
#Survival for GR
gr=ggsurvplot(fitdata1, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
              xlab = "Overall survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
              ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
              legend.labs=c("High", "Low"), font.legend = c(18), legend.title="All Tumors", title="FAK", font.title=c(18, "bold.italic"),
              size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr$table = gr$table + labs(x = NULL, y = NULL)
splots[[1]]=gr

gr2=ggsurvplot(fitdata2, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Relapse free survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="All Tumors", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr2$table = gr2$table + labs(x = NULL, y = NULL)
splots[[2]] = gr2

gr3=ggsurvplot(fitHR1, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Overall survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="HR+", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr3$table = gr3$table + labs(x = NULL, y = NULL)
splots[[3]]=gr3

gr4=ggsurvplot(fitHR2, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Relapse free survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="HR+", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr4$table = gr4$table + labs(x = NULL, y = NULL)
splots[[4]] = gr4

gr5=ggsurvplot(fitTNBC1, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Overall survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="TNBC", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr5$table = gr5$table + labs(x = NULL, y = NULL)
splots[[5]]=gr5

gr6=ggsurvplot(fitTNBC2, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Relapse free survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="TNBC", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr6$table = gr6$table + labs(x = NULL, y = NULL)
splots[[6]] = gr6
gr7=ggsurvplot(fitHER1, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Overall survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="HER2+", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr7$table = gr7$table + labs(x = NULL, y = NULL)
splots[[7]]=gr7

gr8=ggsurvplot(fitHER2, pval=TRUE, risk.table=TRUE, risk.table.col = "strata", conf.int = TRUE, conf.int.style = "step",
               xlab = "Relapse free survival (Months)", break.x.by=15, font.x = c(20, face = "bold"), font.tickslab = c(16, "bold"),
               ylab="Survival Probability", font.y = c(20, face = "bold"), tables.y.text = FALSE,
               legend.labs=c("High", "Low"), font.legend = c(18), legend.title="HER2+", title="FAK", font.title=c(18, "bold.italic"),
               size=1,  ggtheme = centr, surv.median.line = "hv", palette = c("red", "navy"))
gr8$table = gr8$table + labs(x = NULL, y = NULL)
splots[[8]] = gr8
m=arrange_ggsurvplots(splots, print = TRUE, ncol = 4, nrow = 2, risk.table.height = 0.25)
png(file="FAK_survival.png", width=30, height=17, units="in", res=500)
m
dev.off()

p1=ggboxplot(wdata, x = "PTK2Group", y = "Proliferation", xlab="All Tumors", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p2=ggboxplot(SAHR, x = "PTK2Group", y = "Proliferation", xlab="HR+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p3=ggboxplot(SATNBC, x = "PTK2Group", y = "Proliferation", xlab="TNBC", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p4=ggboxplot(SAHER, x = "PTK2Group", y = "Proliferation", xlab="HER2+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p5=ggboxplot(wdata, x = "PTK2Group", y = "Breast.Stemness", xlab="All Tumors", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p6=ggboxplot(SAHR, x = "PTK2Group", y = "Breast.Stemness", xlab="HR+", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p7=ggboxplot(SATNBC, x = "PTK2Group", y = "Breast.Stemness", xlab="TNBC", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p8=ggboxplot(SAHER, x = "PTK2Group", y = "Breast.Stemness", xlab="HER2+", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('GR.Proliferation.Stemness.png', width=25, height=12, units="in", dpi=300)

pa=ggboxplot(wdata, x = "PTK2Group", y = "CYT", xlab="All Tumors", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pb=ggboxplot(SAHR, x = "PTK2Group", y = "CYT", xlab="HR+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pc=ggboxplot(SATNBC, x = "PTK2Group", y = "CYT", xlab="TNBC", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pd=ggboxplot(SAHER, x = "PTK2Group", y = "CYT", xlab="HER2+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pe=ggboxplot(wdata, x = "PTK2Group", y = "EMT", xlab="All Tumors", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pf=ggboxplot(SAHR, x = "PTK2Group", y = "EMT", xlab="HR+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pg=ggboxplot(SATNBC, x = "PTK2Group", y = "EMT", xlab="TNBC", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ph=ggboxplot(SAHER, x = "PTK2Group", y = "EMT", xlab="HER2+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pa,pb,pc,pd,pe,pf,pg,ph, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('EMT.CytolyticActivity.png', width=25, height=12, units="in", dpi=300)

pi=ggboxplot(wdata, x = "PTK2Group", y = "CYT", xlab="All Tumors", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "PTK2Group", y = "CYT", xlab="HR+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "PTK2Group", y = "CYT", xlab="TNBC", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "PTK2Group", y = "CYT", xlab="HER2+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pm=ggboxplot(wdata, x = "PTK2Group", y = "EMT", xlab="All Tumors", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pn=ggboxplot(SAHR, x = "PTK2Group", y = "EMT", xlab="HR+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

po=ggboxplot(SATNBC, x = "PTK2Group", y = "EMT", xlab="TNBC", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pp=ggboxplot(SAHER, x = "PTK2Group", y = "EMT", xlab="HER2+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl,pm,pn,po,pp, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('Angioscore_ifng_score.png', width=25, height=12, units="in", dpi=300)


bind = arrange(df,PATIENT_ID)
write.table(bind,"datafullreadytoanalyse.txt",sep='\t')


pi=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="All Tumors", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="HR+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="TNBC", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="HER2+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pm=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="All Tumors", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pn=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="HR+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

po=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="TNBC", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pp=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="HER2+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl,pm,pn,po,pp, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('Angioscore_ifng_score.png', width=25, height=12, units="in", dpi=300)


p1=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="All Tumors", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p2=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="HR+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p3=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="TNBC", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p4=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="HER2+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p5=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="All Tumors", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p6=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="HR+", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p7=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="TNBC", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p8=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="HER2+", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('GR.Proliferation.Stemness.png', width=25, height=12, units="in", dpi=300)



pi=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="All Tumors", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="HR+", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="TNBC", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="HER2+", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pm=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="All Tumors", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pn=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="HR+", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

po=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="TNBC", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pp=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="HER2+", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl,pm,pn,po,pp, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('Angioscore_ifng_scorebasedonMANO.png', width=25, height=12, units="in", dpi=300)



pi=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="All Tumors", ylab= "PTK2expression", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c(" #00008B", "#FF8C00"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="HR+", ylab= "PTK2expression", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c(" #00008B", "#FF8C00"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="TNBC", ylab= "PTK2expression", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c(" #00008B", "#FF8C00"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="HER2+", ylab= "PTK2expression", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c(" #00008B", "#FF8C00"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl, labels = c("A", "B", "C", "D"), ncol = 4, nrow = 1, common.legend = T)
ggsave('PTK2expressionbasedonMANO.png', width=25, height=12, units="in", dpi=300)