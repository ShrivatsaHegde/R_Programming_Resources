#Calculate 3rd and 1st Quantile of PTK2 expression
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

df = read.table('QuantilePTK2.txt', sep='\t',header = 1)
expression_values <- c(df$PTK2)
third_quantile <- quantile(expression_values, 0.75)
first_quantile <- quantile(expression_values, 0.25)

PTK2_High <- expression_values >= third_quantile
PTK2_Moderate <- expression_values < third_quantile & expression_values > first_quantile
PTK2_Low <- expression_values <= first_quantile



# Create a dataframe
data <- data.frame(PTK2_Expression = expression_values,
                   PTK2_High = PTK2_High,
                   PTK2_Moderate = PTK2_Moderate,
                   PTK2_Low = PTK2_Low)

# Write the dataframe to a CSV file
write.csv(data, file = "QuantilesortedDF.txt",sep='\t' ,row.names = FALSE)

#HR+

dfHR = read.table('HR+Quatilesorted.txt', sep='\t',header = 1)
expression_valuesHR <- c(dfHR$PTK2)
third_quantileHR <- quantile(expression_valuesHR, 0.75)
first_quantileHR <- quantile(expression_valuesHR, 0.25)

PTK2_HighHR <- expression_valuesHR >= third_quantileHR
PTK2_ModerateHR <- expression_valuesHR < third_quantileHR & expression_valuesHR > first_quantileHR
PTK2_LowHR <- expression_valuesHR <= first_quantileHR



# Create a dataframe
dataHR <- data.frame(PTK2_Expression = expression_valuesHR,
                   PTK2_High = PTK2_HighHR,
                   PTK2_Moderate = PTK2_ModerateHR,
                   PTK2_Low = PTK2_LowHR)

# Write the dataframe to a CSV file
write.table(dataHR, file = "QuantilesortedHR+DF.txt", sep='\t' ,row.names = FALSE)



#HER2+

dfHER = read.table('HER2+Quantilesorted.txt', sep='\t',header = 1)
expression_valuesHER <- c(dfHER$PTK2)
third_quantileHER <- quantile(expression_valuesHER, 0.75)
first_quantileHER <- quantile(expression_valuesHER, 0.25)

PTK2_HighHER <- expression_valuesHER >= third_quantileHER
PTK2_ModerateHER <- expression_valuesHER < third_quantileHER & expression_valuesHER > first_quantileHER
PTK2_LowHER <- expression_valuesHER <= first_quantileHER



# Create a dataframe
dataHER <- data.frame(PTK2_Expression = expression_valuesHER,
                     PTK2_High = PTK2_HighHER,
                     PTK2_Moderate = PTK2_ModerateHER,
                     PTK2_Low = PTK2_LowHER)

# Write the dataframe to a CSV file
write.table(dataHER, file = "QuantilesortedHER2+DF.txt", sep='\t' ,row.names = FALSE)



#TNBC

dfTNBC = read.table('TNBCQuantileSorted.txt', sep='\t',header = 1)
expression_valuesTNBC <- c(dfTNBC$PTK2)
third_quantileTNBC <- quantile(expression_valuesTNBC, 0.75)
first_quantileTNBC <- quantile(expression_valuesTNBC, 0.25)

PTK2_HighTNBC <- expression_valuesTNBC >= third_quantileTNBC
PTK2_ModerateTNBC <- expression_valuesTNBC < third_quantileTNBC & expression_valuesTNBC > first_quantileTNBC
PTK2_LowTNBC <- expression_valuesTNBC <= first_quantileTNBC


# Create a dataframe
dataTNBC <- data.frame(PTK2_Expression = expression_valuesTNBC,
                     PTK2_High = PTK2_HighTNBC,
                     PTK2_Moderate = PTK2_ModerateTNBC,
                     PTK2_Low = PTK2_LowTNBC)

# Write the dataframe to a CSV file
write.table(dataTNBC, file = "QuantilesortedTNBCDF.txt", sep='\t' ,row.names = FALSE)

