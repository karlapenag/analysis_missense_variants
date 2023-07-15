library(NMF)  # load NMF package
library(pheatmap)
library(gplots)
library(openxlsx)
library(readxl)
packages <- c("readxl","graphics", "stats", "stringr", "digest", "grid", "grDevices", "gridBase", "colorspace", "RColorBrewer", "foreach", "doParallel", "ggplot2", "reshape2", "codetools", "BiocManager")

lapply(packages, library, character.only = TRUE)
library(dplyr)

epileptic_encephalopathy_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/epileptic_encephalopathy_input.xlsx")
cancer_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/cancer_input.xlsx")
heart_disease_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/heart_disease_input.xlsx")
intellectual_disability_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/intellectual_disability_input.xlsx")
charcot_marie_tooth_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/charcot_marie_tooth_input.xlsx")
hypercholesterolemia_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/hypercholesterolemia_input.xlsx")
marfan_syndrome_data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/marfan_syndrome_input.xlsx")

## RUN FOR EVERY DISEASE ONE BY ONE
data <- epileptic_encephalopathy_data
data <- cancer_data
data <- heart_disease_data
data <- intellectual_disability_data
data <- charcot_marie_tooth_data
data <- hypercholesterolemia_data
data <- marfan_syndrome_data

## RUN ONE BYE ONE

data <- data %>% 
  rename("id" = "...1")

continuous_mat <- as.data.frame(data)

# Assign the 'id' column as the row.names
row.names(continuous_mat) <- continuous_mat$id
continuous_mat$id <- NULL # Remove the 'id' column from the data frame

continuous_mat[continuous_mat == 0] <- continuous_mat[continuous_mat == 0] + 0.0001


#Specify the data types for each column
column_types <- c("numeric","text","numeric", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text",
                  "text", "text", "text", "text", "text",
                  "text", "text", "text", "text", "text",
                  "text", "text", "text", "text", "text",
                  "text", "text", "text", "text", "text",
                  "text", "text", "text", "text", "text",
                  "text", "text")

# Read the Excel file and specify the column types
epileptic_encephalopathy_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/epileptic_encephalopathy_labels.xlsx", col_types = column_types))
cancer_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/cancer_labels.xlsx", col_types = column_types))
intellectual_disability_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/intellectual_disability_labels.xlsx", col_types = column_types))
heart_disease_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/heart_disease_labels.xlsx", col_types = column_types))
charcot_marie_tooth_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/charcot_marie_tooth_labels.xlsx", col_types = column_types))
marfan_syndrome_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/marfan_syndrome_labels.xlsx", col_types = column_types))
hypercholesterolemia_anns <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/diseases/hypercholesterolemia_labels.xlsx", col_types = column_types))

## RUN FOR EVERY DISEASE ONE BY ONE.
annotations <- epileptic_encephalopathy_anns
annotations <- cancer_anns
annotations <- intellectual_disability_anns
annotations <- heart_disease_anns
annotations <- charcot_marie_tooth_anns
annotations <- marfan_syndrome_anns
annotations <- hypercholesterolemia_anns

## RUN ONE BY ONE

annotations <- annotations %>% 
  rename("id" = "...1")

# Assign the 'id' column as the row.names
row.names(annotations) <- annotations$id
annotations$id <- NULL # Remove the 'id' column from the data frame

pathogenicity_df <- data.frame(CHR =annotations$CHR,
                               row.names = row.names(annotations))

## TRANSPOSE
transposed_input = t(continuous_mat)
transposed_labels = t(pathogenicity_df)

estim.r <- nmf(continuous_mat, method='brunet', 2:10, nrun=10, seed=123456)

plot(estim.r)

consensusmap(estim.r)

saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/epileptic_encephalopathy/estimr_results2906_1.rds')
saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/estimr_results2906_1.rds')
saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/intellectual_disability/estimr_results2906_1.rds')
saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/estimr_results2906_1.rds')

k <- 9

# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)

nmf_fit <- nmf(transposed_input, rank=k, method="brunet", nrun=4, seed=123456)

pathogenicity_df$CHR<-as.character(pathogenicity_df$CHR)
custom_chr <- c("#00FF99","#66FFCC","#33CC99","#33FFFF","#33CCCC","#339999","#336666","#006699","#003399","#3333ff",
                "#3333CC","#333399","#333366","#6633CC","#9966FF","#6600FF","#FF00CC","#FF33CC","#990066","#CC0066",
                "#FF3399","#FF9999","#CC9999")

nmf_fit <- readRDS(file = "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/nmf_fit2906_rank9_1.rds")

png("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/coefmap_2906_rank9_1.png", width = 10, height = 10, units = "in", res = 300)
coefmap(nmf_fit)
dev.off()

png("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/basismap_2906_rank9_1.png", width = 10, height = 10, units = "in", res = 300)
basismap(nmf_fit)
dev.off()


saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/epileptic_encephalopathy/nmf_fit2906_rank4_1.rds')
saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/nmf_fit2906_rank4_1.rds')
saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/intellectual_disability/nmf_fit2906_rank7_1.rds')
saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/nmf_fit2906_rank9_1.rds')


# normalize rows of the coef matrix
coef_norm <- t(apply(coef(nmf_fit), 1, function(x) x/sum(x)))

# print normalized matrix
coef_norm

#### Add covariates/ annotations Pathogenic/Benign to variants
nmf_results <- cbind(annotations, t(coef(nmf_fit)))

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", nmf_results)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/epileptic_encephalopathy/nmf_results_coeflabels_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/nmf_results_coeflabels_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/intellectual_disability/nmf_results_coeflabels_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/nmf_results_coeflabels_2906_2.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef_norm)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/epileptic_encephalopathy/coef_norm_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/coef_norm_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/intellectual_disability/coef_norm_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/coef_norm_2906_2.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/epileptic_encephalopathy/coef_fit_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/coef_fit_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/intellectual_disability/coef_fit_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/coef_fit_2906_2.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/epileptic_encephalopathy/basis_fit_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/cancer/basis_fit_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/intellectual_disability/basis_fit_2906_1.xlsx")
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/5/heart_disease/basis_fit_2906_2.xlsx")

#### proving normalization

colSums(basis(nmf_fit))
rowSums(coef_norm)

plot(silhouette(nmf_fit))
si <- silhouette(nmf_fit, what = 'features')
plot(si)
