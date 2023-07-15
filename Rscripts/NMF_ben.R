library(NMF)  # load NMF package
library(pheatmap)
library(gplots)
library(openxlsx)
library(readxl)
packages <- c("readxl","graphics", "stats", "stringr", "digest", "grid", "grDevices", "gridBase", "colorspace", "RColorBrewer", "foreach", "doParallel", "ggplot2", "reshape2", "codetools", "BiocManager")

lapply(packages, library, character.only = TRUE)
library(dplyr)

data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/benign_complete_input.xlsx")

data <- data %>% 
  rename("id" = "...1")

continuous_mat <- as.data.frame(data)

# Assign the 'id' column as the row.names
row.names(continuous_mat) <- continuous_mat$id
continuous_mat$id <- NULL # Remove the 'id' column from the data frame

#Specify the data types for each column
column_types <- c("numeric","text","numeric", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text", "text")

# Read the Excel file and specify the column types
annotations <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/benign_complete_labels.xlsx", col_types = column_types))

annotations <- annotations %>% 
  rename("id" = "...1")

# Assign the 'id' column as the row.names
row.names(annotations) <- annotations$id
annotations$id <- NULL # Remove the 'id' column from the data frame

pathogenicity_df <- data.frame(Pathogenicity = annotations$Pathogenicity,
                               CHR =annotations$CHR ,
                               row.names = row.names(annotations))


transposed_input = t(continuous_mat)
transposed_labels = t(pathogenicity_df)

estim.r <- nmf(continuous_mat, method='brunet', 2:10, nrun=10, seed=123456)

plot(estim.r)

consensusmap(estim.r)

# Plot consensus matrix
plot(consensus_mat)

saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/estimr_2806_1.rds')

k <- 4
# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)
nmf_fit <- nmf(transposed_input, rank=k, method="brunet", nrun=4, seed=123456)

pathogenicity_df$CHR<-as.character(pathogenicity_df$CHR)

saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank4/nmf_fit2806_rank4_1.rds')

custom_colors <- c("benign" = "green","pathogenic" = "red")
custom_chr <- c("#00FF99","#66FFCC","#33CC99","#33FFFF","#33CCCC",
                "#339999","#336666","#006699","#003399","#3333ff",
                "#3333CC","#333399","#333366","#6633CC","#9966FF",
                "#6600FF","#FF00CC","#FF33CC","#990066","#CC0066",
                "#FF3399","#FF9999","#CC9999")

basismap(nmf_fit)

coefmap(nmf_fit,annCol=pathogenicity_df, annColors = list(Pathogenicity=custom_colors,
                                                          CHR = custom_chr))
# normalize rows of the coef matrix
coef_norm <- t(apply(coef(nmf_fit), 1, function(x) x/sum(x)))

#### Add covariates/ annotations Pathogenic/Benign to variants
nmf_results <- cbind(annotations, t(coef(nmf_fit)))

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", nmf_results)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank4/nmf_results_coeflabels_2806_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef_norm)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank4/coef_norm_2806_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank4/coef_fit_2806_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank4/basis_fit_2806_1.xlsx")

#### proving normalization

colSums(basis(nmf_fit))

rowSums(coef_norm)

plot(silhouette(nmf_fit))
si <- silhouette(nmf_fit, what = 'consensus')
plot(si)

k <- 2
# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)
nmf_fit <- nmf(continuous_mat, rank=k, method="brunet", nrun=10, seed=123456)

basismap(nmf_fit)

coefmap(nmf_fit)

saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank2/nmf_fit2806_rank7_1.rds')

# normalize rows of the coef matrix
coef_norm <- t(apply(coef(nmf_fit), 1, function(x) x/sum(x)))

# print normalized matrix
coef_norm

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef_norm)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank2/coef_norm_2806_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank2/coef_fit_2806_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/3/rank2/basis_fit_2806_1.xlsx")

#### proving normalization

colSums(basis(nmf_fit))

rowSums(coef_norm)

plot(silhouette(nmf_fit))
si <- silhouette(nmf_fit, what = 'consensus')
plot(si)