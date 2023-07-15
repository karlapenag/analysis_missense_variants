library(NMF)  # load NMF package
library(pheatmap)
library(gplots)
library(openxlsx)
library(readxl)
packages <- c("readxl","graphics", "stats", "stringr", "digest", "grid", "grDevices", "gridBase", "colorspace", "RColorBrewer", "foreach", "doParallel", "ggplot2", "reshape2", "codetools", "BiocManager")

lapply(packages, library, character.only = TRUE)
library(dplyr)

data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/variants_combined_input.xlsx")

data <- data %>% 
  rename("id" = "...1")

continuous_mat <- as.data.frame(data)

# Assign the 'id' column as the row.names
row.names(continuous_mat) <- continuous_mat$id
continuous_mat$id <- NULL # Remove the 'id' column from the data frame

#Specify the data types for each column
column_types <- c("numeric","text","numeric", "text", "text", "text", 
                  "text", "text", "text","text", "text", 
                  "text", "text", "text", "text", "text", 
                  "text", "text", "text", "text")
# Read the Excel file and specify the column types
annotations <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/variants_combined_labels.xlsx", col_types = column_types))

annotations <- annotations %>% 
  rename("id" = "...1")

# Assign the 'id' column as the row.names
row.names(annotations) <- annotations$id
annotations$id <- NULL # Remove the 'id' column from the data frame

pathogenicity_df <- data.frame(Pathogenicity = annotations$Pathogenicity,
                               row.names = row.names(annotations))

transposed_input = t(continuous_mat)
transposed_labels = t(pathogenicity_df)

estim.r <- nmf(continuous_mat, method='brunet', 2:10, nrun=10, seed=123456)

plot(estim.r)

consensusmap(estim.r)
consens <- consensus(nmf_fit)

saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/estimr_results2806_1.rds')

k <- 2
# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)
nmf_fit <- nmf(continuous_mat, rank=k, method="brunet", nrun=5, seed=123456)
transposed_nmf<- t(nmf_fit)

basismap(transposed_nmf)

custom_colors <- c("pathogenic" = "red","benign" = "green")

coefmap(transposed_nmf, annCol=pathogenicity_df, annColors = list(Pathogenicity=custom_colors), tracks="basis")

saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank2/nmf_fit2806_rank2_1.rds')

k <- 6
# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)
nmf_fit <- nmf(continuous_mat, rank=k, method="brunet", nrun=10, seed=123456)
transposed_nmf<- t(nmf_fit)

basismap(transposed_nmf)

coefmap(transposed_nmf, annCol=pathogenicity_df, annColors = list(Pathogenicity=custom_colors),tracks="basis")

saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank6/nmf_fit2906_rank6_1.rds')

#  without columns with almost 0 std

estim.r <- readRDS(file = "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/estimr_results2806_1.rds")

# normalize rows of the coef matrix
coef_norm <- t(apply(coef(nmf_fit), 1, function(x) x/sum(x)))

#### Add covariates/ annotations Pathogenic/Benign to variants
nmf_results <- cbind(annotations, basis(nmf_fit))

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", nmf_results)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank6/nmf_factors_coeflabels_rank6_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef_norm)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank6/coef_norm_2906_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(transposed_nmf))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank6/coef_fit_2906_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(transposed_nmf))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank6/basis_fit_2906_1.xlsx")

#### proving normalization

colSums(basis(nmf_fit))

rowSums(coef_norm)

plot(silhouette(nmf_fit))
si <- silhouette(transposed_nmf, what = 'features')
plot(si)

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", nmf_results)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank6/nmf_factors_labels_rank6_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef_norm)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank5/coef_norm_2606_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank5/coef_fit_2606_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank5/basis_fit_2606_1.xlsx")

#### proving normalization

colSums(basis(nmf_fit))

rowSums(coef_norm)

nmf_fit <- readRDS(file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/1/rank5/nmf_fit2606_rank5_1.rds')
plot(silhouette(nmf_fit))
si <- silhouette(nmf_fit, what = 'consensus')
plot(si)
