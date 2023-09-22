setwd(R"(C:\Dessertation\metabric_analysis)")
list.files()

df_tnbc=read.table('TNBCstatusDEGS.txt', sep='\t',header = 1)

filtered_df <- subset(df_tnbc, df_tnbc$P.Value <= 0.05 & df_tnbc$logFC < 1)

# Sort the filtered data frame based on the 'p values' column
sorted_df <- filtered_df[order(filtered_df$P.Value), ]
write.table(sorted_df,"DEGSselectedinTNBC.txt",sep = '\t',col.names = TRUE ,row.names = TRUE)


df_HR=read.table('HRstatusDEGS.txt', sep='\t',header = 1)
HR_df <- subset(df_HR, df_HR$P.Value <= 0.05 & df_HR$logFC < 1)
write.table(HR_df,"DEGSselectedinHR+.txt",sep = '\t',col.names = TRUE ,row.names = TRUE)



df_HER=read.table('HER2statusDEGS.txt', sep='\t',header = 1)
HER_df <- subset(df_HER, df_HER$P.Value <= 0.05 & df_HER$logFC < 1)
write.table(HER_df,"DEGSselectedinHER2+.txt",sep = '\t',col.names = TRUE ,row.names = TRUE)

