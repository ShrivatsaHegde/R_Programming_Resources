setwd(R"(C:\Dessertation\metabric_analysis\Ana;ysisExtra\Shrivatsa)")
#set working directory
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
library(imogimap)
library("survival")
library("survminer")
#All the packages needed
list.files()
#list out the files in the folder

df = read.table('illumina_microarraywithoutentrenz.txt', sep='\t',header = 1)
# the whole dataset with genes as rows and patients as colnames
row_names <- df[, 1]
#set the row names as gene names
row.names(df) <- make.unique(row_names)
# unique for making any genes which are not readable to readabel
sub_df <- df[, -c(1)]
# remove the first coloumn which contains the gene names
# fuction for EMT scoring
get_emt_score=function(pdf){
    pdf <- as.matrix(pdf)
    missing_EMT <- EMT_gene_list[-which(EMT_gene_list$genes %in% rownames(pdf)),]$gene
    if(length(missing_EMT)>0){
      warning(length(missing_EMT)," missing EMT signature genes:  \n  ",
              lapply(missing_EMT, function(x)paste0(x,"  ")),"\nCheck EMT_gene_list for all signature genes.\n")
    }
    
    sub_list <- EMT_gene_list[!(EMT_gene_list$genes %in% missing_EMT),]
    pdf_sub <- pdf[rownames(pdf) %in% sub_list$genes,]
    if(nrow(pdf_sub)==0){
      warning("No signature genes found!")
      pscore <-tibble("score"=numeric())
      
    }else{
      pdf_sub <- log2(pdf_sub+1)
      pdf_sub <- as.data.frame( t(scale(t(pdf_sub), center = T, scale = T)))
      pdf_sub <- sweep(pdf_sub, 1, sub_list$sign, "*")
      
      pscore <-as.data.frame(colMeans(pdf_sub,na.rm = T))
      colnames(pscore) <- "EMTscore"
      #pscore$EMTscore <- 2^pscore$EMTscore
      pscore$Tumor_Sample_ID <- rownames(pscore)
      pscore <- pscore[,c("Tumor_Sample_ID","EMTscore")]
    }
    return(pscore)
  }
# give the cleaned dataset to this through this function
EMTscore<-get_emt_score(sub_df)

# Function for getting ANgioScore

get_angio_score=function(pdf){
  pdf <- as.matrix(pdf)
  missing_genes <- AG_gene_list[-which(AG_gene_list %in% rownames(pdf))]
  if(length(missing_genes)>0){
    warning(length(missing_genes)," missing Angiogenesis signature genes:  \n  ",
            lapply(missing_genes, function(x)paste0(x,"  ")),"\n Check AG_gene_list for
      all signature genes.\n")
  }
  
  sub_list <- AG_gene_list[!(AG_gene_list %in% missing_genes)]
  pdf_sub <- pdf[rownames(pdf) %in% sub_list,]
  if(nrow(pdf_sub)==0){
    warning("No signature genes found!")
    pscore <-tibble("score"=numeric())
    
  }else{
    pdf_sub <- log2(pdf_sub+1)
    pdf_sub <- as.data.frame( t(scale(t(pdf_sub), center = T, scale = T)))
    
    pscore <-as.data.frame(colMeans(pdf_sub,na.rm = T))
    colnames(pscore) <- "AGscore"
    #pscore$AGscore <- 2^pscore$AGscore
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","AGscore")]
  }
  return(pscore)
}

# give the data to that function
Angioscore<-get_angio_score(sub_df)

# defining the function for getting the IFNG score
get_ifng_score=function(pdf){
  pdf <- as.matrix(pdf)
  IFNG_gene_list <- "IFNG"
  missing_IFNG <- IFNG_gene_list[-which(IFNG_gene_list %in% rownames(pdf))]
  if(length(missing_IFNG)>0){
    warning(length(missing_IFNG)," missing IFNG genes:",
            lapply(missing_IFNG, function(x)paste0(x,"  ")),"\n")
  }
  
  sub_list <- IFNG_gene_list[!(IFNG_gene_list %in% missing_IFNG)]
  pdf_sub <- as.matrix(pdf[rownames(pdf) %in% sub_list,])
  if(nrow(pdf_sub)==0){
    warning("No signature genes found!")
    pscore <-tibble("score"=numeric())
    
  }else{
    pdf_sub <- log2(pdf_sub+1)
    pscore <- as.data.frame( scale((pdf_sub), center = T, scale = T) )
    
    colnames(pscore) <- "IFNGscore"
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","IFNGscore")]
  }
  return(pscore)
}
# pass the data from this fuction 
IFNGscore<- get_ifng_score(sub_df)
# read the file which are having all the clinical data required in a variable
targets = read.table('clinical_patient.txt', sep='\t',header = 1)
# just creating a dataset called bind which will have all the information  about the patients
bind<-cbind(targets, EMTscore$EMTscore)
colnames(bind)[colnames(bind) == "EMTscore$EMTscore"] <- "EMTscore"
# this steps ensure that the EMT score and Angio score and everything will be included in bind
bind<-cbind(bind, Angioscore$AGscore)
colnames(bind)[colnames(bind) == "Angioscore$AGscore"] <- "AGscore"
bind<-cbind(bind, IFNGscore$IFNGscore)
colnames(bind)[colnames(bind) == "IFNGscore$IFNGscore"] <- "IFNGscore"
# then take the data of PTK2 or the gene you are intrested and include that in bind
target_row <- df["PTK2", ]
target_column <- t(target_row)
target_column <- target_column[-1, ]
bind <- cbind(bind, target_column)
colnames(bind)[colnames(bind) == "target_column"] <- "PTK2"#renaming the colounm as the selected gene names
bind$OS_STATUS <- gsub(pattern = ":LIVING", replacement = "", x = bind$OS_STATUS)
bind$OS_STATUS <- gsub(pattern = ":DECEASED", replacement = "", x = bind$OS_STATUS)
bind$RFS_STATUS <- gsub(pattern = ":Not Recurred", replacement = "", x = bind$RFS_STATUS)
bind$RFS_STATUS <- gsub(pattern = ":Recurred", replacement = "", x = bind$RFS_STATUS)
#converting the vales to numeric value to fit them into the survival analysis
bind$OS_STATUS <- as.numeric(bind$OS_STATUS)
bind$RFS_STATUS <- as.numeric(bind$RFS_STATUS)
# Proliferation markers
markers<- sub_df["MKI67",]
Pmark<- sub_df["MCM6",]
Pm<- sub_df["PCNA",]
# transpose of the data taken 
markers_col<- t(markers)
Pmark_col<- t(Pmark)
Pm_col<- t(Pm)
bind <- cbind(bind, markers_col)
colnames(bind)[colnames(bind) == "markers_col"] <- "MKI67"
bind <- cbind(bind, Pmark_col)
colnames(bind)[colnames(bind) == "Pmark_col"] <- "MCM6"
bind <- cbind(bind, Pm_col)
colnames(bind)[colnames(bind) == "Pm_col"] <- "PCNA"
bind$Proliferation_mean <- rowMeans(bind[, c("MKI67", "MCM6", "PCNA")])
# Stemness Markers
Stemness<- sub_df["EPCAM",]
Smarker<-sub_df["CD24",]
marker<-sub_df["BMI1",]
Stem<-sub_df["ITGA6",]
ness<-sub_df["KLF4",]
Sma<-sub_df["CXCR1",]
Mark<-sub_df["SALL4",]
kers<-sub_df["SOX2",]
ke<-sub_df["MYC",]
ne<-sub_df["THY1",]
ers<-sub_df["CXCR4",]
#transpose 
Stemness_col<- t(Stemness)
Smarker_col<-t(Smarker)
marker_col<- t(marker)
Stem_col<- t(Stem)
ness_col<- t(ness)
Sma_col<- t(Sma)
Mark_col<- t(Mark)
kers_col<- t(kers)
ke_col<- t(ke)
ne_col<- t(ne)
ers_col<-t(ers)
#to the main dataframe
bind <- cbind(bind, Stemness_col)
bind <- cbind(bind, marker_col)
bind <- cbind(bind, Stem_col)
bind <- cbind(bind, ness_col)
bind <- cbind(bind, Sma_col)
bind <- cbind(bind, Mark_col)
bind <- cbind(bind, kers_col)
bind <- cbind(bind, Smarker_col)
bind <- cbind(bind, ke_col)
bind <- cbind(bind, ne_col)
bind <- cbind(bind, ers_col)
bind$Stemness_mean <- rowMeans(bind[, c("EPCAM", "CD24", "BMI1" ,"ITGA6", "KLF4","CXCR1","SALL4","SOX2","MYC","THY1","CXCR4")])
extra = read.table('datafullreadytoanalyse.txt', sep='\t',header = 1)
bind <- cbind(bind, extra$CYT)
colnames(bind)[colnames(bind) == "extra$CYT"] <- "CYT"
bind<- cbind(bind,extra$Breast.Stemness)
colnames(bind)[colnames(bind) == "extra$Breast.Stemness"] <- "Breast.Stemness"
group<-read.table('GROUPSandsample.txt', sep='\t',header = 1)
bind<-cbind(bind, group$Groups)
colnames(bind)[colnames(bind) == "group$Groups"] <- "Groups"
bind<-cbind(bind, group$status)#just getting the group and status data . which is nothing but HR+,HER,TNBC data.
colnames(bind)[colnames(bind) == "group$status"] <- "status"
# Grouping Again
HRgroup<-bind[bind$Groups == "HR+", ]
HRgroup<-HRgroup[, -c(1)]
HERgroup<-bind[bind$Groups == "HER2+", ]
HERgroup<-HERgroup[, -c(1)]
TNBCgroup<-bind[bind$Groups == "TNBC", ]
TNBCgroup<-TNBCgroup[, -c(1)]
library(ggpubr)
library(RColorBrewer)
bind<-bind[, -c(1)]
bind$PTK2 <- as.numeric(bind$PTK2)
q1 = quantile(bind$PTK2,0.33)
q2 = quantile(bind$PTK2,0.66)
bind$PTK2Group = with(bind,
                    ifelse(PTK2 <= q1, 'Low',
                           ifelse(PTK2 > q1 & PTK2 < q2, 'Moderate',
                                  ifelse(PTK2 >= q2, 'High', 'unknown'))))
dim(TNBCgroup)
dim(HRgroup)
dim(HERgroup)
table(bind$PTK2Group)
TNBCgroup$PTK2 <- as.numeric(TNBCgroup$PTK2)
HRgroup$PTK2 <- as.numeric(HRgroup$PTK2)
HERgroup$PTK2 <- as.numeric(HERgroup$PTK2)
q3 = quantile(TNBCgroup$PTK2,0.33)
q4 = quantile(TNBCgroup$PTK2,0.66)
TNBCgroup$PTK2Group = with(TNBCgroup,
                      ifelse(PTK2 <= q3, 'Low',
                             ifelse(PTK2 > q3 & PTK2 < q4, 'Moderate',
                                    ifelse(PTK2 >= q4, 'High', 'unknown'))))
q5 = quantile(HRgroup$PTK2,0.33)
q6 = quantile(HRgroup$PTK2,0.66)
HRgroup$PTK2Group = with(HRgroup,
                       ifelse(PTK2 <= q5, 'Low',
                              ifelse(PTK2 > q5 & PTK2 < q6, 'Moderate',
                                     ifelse(PTK2 >= q6, 'High', 'unknown'))))
q7 = quantile(HERgroup$PTK2,0.33)
q8 = quantile(HERgroup$PTK2,0.66)
HERgroup$PTK2Group = with(HERgroup,
                      ifelse(PTK2 <= q7, 'Low',
                             ifelse(PTK2 > q7 & PTK2 < q8, 'Moderate',
                                    ifelse(PTK2 >= q8, 'High', 'unknown'))))

write.table(TNBCgroup, "tnbc_metabric.PTK2.signature.txt", sep='\t')
write.table(HRgroup, "hrpos_metabric.PTK2.signature.txt", sep='\t')
write.table(HERgroup, "Her2_metabric.PTK2.signature.txt", sep='\t')
centr = theme_grey() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
wdata = filter(bind, PTK2Group == "High" | PTK2Group == "Low")
wdata= arrange(wdata, PTK2Group)
wdata_M = filter(wdata, OS_MONTHS <= 120)
SAHR = filter(HRgroup, PTK2Group == "High" | PTK2Group == "Low")
SAHR = arrange(SAHR, PTK2Group)
SAHR_M = filter(SAHR, OS_MONTHS <= 120)
SAHER = filter(HERgroup, PTK2Group == "High" | PTK2Group == "Low")
SAHER = arrange(SAHER, PTK2Group)
SAHER_M = filter(SAHER, OS_MONTHS <= 120)
SATNBC = filter(TNBCgroup, PTK2Group == "High" | PTK2Group == "Low")
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
png(file="PTK2_survival.png", width=30, height=17, units="in", res=500)
m
dev.off()
p1=ggboxplot(wdata, x = "PTK2Group", y = "Proliferation_mean", xlab="All Tumors", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p2=ggboxplot(SAHR, x = "PTK2Group", y = "Proliferation_mean", xlab="HR+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p3=ggboxplot(SATNBC, x = "PTK2Group", y = "Proliferation_mean", xlab="TNBC", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p4=ggboxplot(SAHER, x = "PTK2Group", y = "Proliferation_mean", xlab="HER2+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
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

pe=ggboxplot(wdata, x = "PTK2Group", y = "EMTscore", xlab="All Tumors", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pf=ggboxplot(SAHR, x = "PTK2Group", y = "EMTscore", xlab="HR+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pg=ggboxplot(SATNBC, x = "PTK2Group", y = "EMTscore", xlab="TNBC", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ph=ggboxplot(SAHER, x = "PTK2Group", y = "EMTscore", xlab="HER2+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pa,pb,pc,pd,pe,pf,pg,ph, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('EMT.CytolyticActivity.png', width=25, height=12, units="in", dpi=300)
# Agioscore and IFNG score


pi=ggboxplot(wdata, x = "PTK2Group", y = "AGscore", xlab="All Tumors", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "PTK2Group", y = "AGscore", xlab="HR+", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "PTK2Group", y = "AGscore", xlab="TNBC", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "PTK2Group", y = "AGscore", xlab="HER2+", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pm=ggboxplot(wdata, x = "PTK2Group", y = "IFNGscore", xlab="All Tumors", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pn=ggboxplot(SAHR, x = "PTK2Group", y = "IFNGscore", xlab="HR+", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

po=ggboxplot(SATNBC, x = "PTK2Group", y = "IFNGscore", xlab="TNBC", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pp=ggboxplot(SAHER, x = "PTK2Group", y = "IFNGscore", xlab="HER2+", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="PTK2Group", palette = c("red", "navy"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl,pm,pn,po,pp, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('Angioscore_ifng_score.png', width=25, height=12, units="in", dpi=300)

pi=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="All Tumors", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="HR+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="TNBC", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "CYT", xlab="HER2+", ylab= "Cytolytic Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pm=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="All Tumors", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pn=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="HR+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

po=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="TNBC", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pp=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "EMTscore", xlab="HER2+", ylab= "EMT Activity", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl,pm,pn,po,pp, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('EMTscoreandCYTbasedonMANO.png', width=25, height=12, units="in", dpi=300)


p1=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="All Tumors", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p2=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="HR+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p3=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="TNBC", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p4=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "Proliferation_mean", xlab="HER2+", ylab= "Proliferation", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p5=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="All Tumors", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p6=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="HR+", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p7=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="TNBC", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =  c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

p8=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "Breast.Stemness", xlab="HER2+", ylab= "Stemness", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('Proliferation.StemnessbasedonMANO.png', width=25, height=12, units="in", dpi=300)

pi=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="All Tumors", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="HR+", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="TNBC", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "AGscore", xlab="HER2+", ylab= "Angioscore", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pm=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="All Tumors", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pn=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="HR+", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("blue", "orange"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

po=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="TNBC", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pp=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "IFNGscore", xlab="HER2+", ylab= "Interferon_gamma", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("orange","blue"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl,pm,pn,po,pp, labels = c("A", "B", "C", "D", "E", "F","G","H"), ncol = 4, nrow = 2, common.legend = T)
ggsave('Angioscore_ifng_scorebasedonMANO.png', width=25, height=12, units="in", dpi=300)
#PTK2expression graph
pi=ggboxplot(wdata, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="All Tumors", ylab= "PTK2expression", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("#00008B", "#FF8C00"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pj=ggboxplot(SAHR, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="HR+", ylab= "PTK2expression", legend = "none", width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("#00008B", "#FF8C00"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pk=ggboxplot(SATNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="TNBC", ylab= "PTK2expression", legend = c("Pre","Post"), width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette = c("#FF8C00","#00008B"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

pl=ggboxplot(SAHER, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2", xlab="HER2+", ylab= "PTK2expression", legend = c("Pre","Post"), width =0.5, size=0.8,
             color="INFERRED_MENOPAUSAL_STATE", palette =c("#FF8C00","#00008B"), bxp.errorbar = TRUE) + stat_compare_means(label = "p.format", label.x = 1.5, face="bold", size=6) + 
  theme_grey() + font("xylab", size = 15, face = "bold") + font("xy.text", size = 14, color = "#606060")

ggarrange(pi,pj,pk,pl, labels = c("A", "B","C", "D"), ncol = 4, nrow = 1, common.legend = T)
ggsave('PTK2expressionbasedon2MANO.png', width=20, height=12, units="in", dpi=300)



