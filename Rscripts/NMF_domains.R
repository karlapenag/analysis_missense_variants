library(NMF)  # load NMF package
library(pheatmap)
library(gplots)
library(openxlsx)
library(readxl)
library(Matrix)
packages <- c("readxl","graphics", "stats", "stringr", "digest", "grid", "grDevices", "gridBase", "colorspace", "RColorBrewer", "foreach", "doParallel", "ggplot2", "reshape2", "codetools", "BiocManager")

lapply(packages, library, character.only = TRUE)
library(dplyr)

data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/domains/ion_transport_input.xlsx")
data <- read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/domains/fibrillin_1_input.xlsx")

data <- data %>% 
  rename("id" = "...1")

continuous_mat <- as.data.frame(data)


# Assign the 'id' column as the row.names
row.names(continuous_mat) <- continuous_mat$id
continuous_mat$id <- NULL # Remove the 'id' column from the data frame

continuous_mat[continuous_mat == 0] <- continuous_mat[continuous_mat == 0] + 0.000001

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
                  "text", "text", "text", "text", "text","text")

# Read the Excel file and specify the column types
annotations <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/domains/ion_transport_labels.xlsx", col_types = column_types))
annotations <- as.data.frame(read_excel("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/domains/fibrillin_1_labels.xlsx", col_types = column_types))

annotations <- annotations %>%
  rename("id" = "...1")

# Assign the 'id' column as the row.names
row.names(annotations) <- annotations$id
annotations$id <- NULL # Remove the 'id' column from the data frame

pathogenicity_df <- data.frame(Pathogenicity = annotations$Pathogenicity,
                               Heart_disease = annotations$heart_disease,
                               Intellectual_disability = annotations$intellectual_disability,
                               Epileptic_encephalopathy = annotations$epileptic_encephalopathy,
                               CHR =annotations$CHR,
                               row.names = row.names(annotations))


transposed_input = t(continuous_mat)
transposed_labels = t(pathogenicity_df)

estim.r <- nmf(continuous_mat, method='brunet', 2:10, nrun=10, seed=123456)

plot(estim.r)

consensus_mat <- consensusmap(estim.r)

# Plot consensus matrix
plot(consensus_mat)

saveRDS(estim.r, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/estimr_results0707_ion_1.rds')

k <- 7
# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)
nmf_fit <- nmf(transposed_input, rank=k, method="brunet", nrun=4, seed=123456)

coef_matrix <- coef(nmf_fit)
saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/nmf_fit2706_rank7_1.rds')
nmf_fit <- readRDS(file = "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/ion/rank4/nmf_fit2706_rank4_1.rds")

pathogenicity_df$Disease_count<-as.numeric(pathogenicity_df$Disease_count)
pathogenicity_df$CHR<-as.character(pathogenicity_df$CHR)

custom_colors <- c("pathogenic" = "red","benign" = "green")
custom_heart <- c("heart_disease" = "coral", 
                  "other" = "white")
custom_marfan <- c("marfan_syndrome" = "brown", 
                   "other" = "white")
custom_epi <- c("epileptic_encephalopathy" = "purple", 
                "other" = "white")
custom_intel <- c("intellectual_disability" = "lightblue", 
                  "other" = "white")

custom_chr <- c("#00FF99","#6633CC","#FF00CC","#FF9999","#3333ff","#66FFCC","#33CC99","#33FFFF","#33CCCC",
                "#339999","#336666","#006699","#003399",
                "#3333CC","#333399","#333366","#9966FF",
                "#6600FF","#FF33CC","#990066","#CC0066",
                "#FF3399","#CC9999")

basismap(nmf_fit)

coefmap(nmf_fit,annCol=pathogenicity_df, annColors = list(Pathogenicity=custom_colors,
                                                          Intellectual_disability = custom_intel,
                                                          Heart_disease = custom_heart,
                                                          Marfan_syndrome = custom_marfan,
                                                          CHR = custom_chr))
png("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/coefmap_0707_rank4_1.png", width = 10, height = 10, units = "in", res = 300)
coefmap(nmf_fit,annCol=pathogenicity_df, annColors = list(Pathogenicity=custom_colors,
                                                          Intellectual_disability = custom_intel,
                                                          Heart_disease = custom_heart,
                                                          Epileptic_encephalopathy = custom_epi,
                                                          CHR = custom_chr))
dev.off()

png("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/basismap_0707_rank7_1.png", width = 10, height = 10, units = "in", res = 300)
basismap(nmf_fit)
dev.off()


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
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/ion/nmf_results_coeflabels_2806_7.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef_norm)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/ion/coef_norm_2706_7.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/ion/coef_fit_2706_7.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/ion/basis_fit_2706_7.xlsx")

#### proving normalization

colSums(basis(nmf_fit))

rowSums(coef_norm)

plot(silhouette(nmf_fit))
si <- silhouette(nmf_fit, what = 'features')
plot(si)
png("/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/silhouette_rank7_1.png", width = 10, height = 10, units = "in", res = 300)
plot(si)
dev.off()
# feature clustering

k <- 5
# run NMF with a specific algorithm (e.g., "brunet" algorithm is default)
nmf_fit <- nmf(continuous_mat, rank=k, method="brunet", nrun=10, seed=123456)

basismap(nmf_fit, annRow=pathogenicity_df, annColors = list(Pathogenicity=custom_colors))
basismap(nmf_fit)

coefmap(nmf_fit)

saveRDS(nmf_fit, file = '/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/rank5/nmf_fit2706_rank5_1.rds')

#### Add covariates/ annotations Pathogenic/Benign to variants
nmf_results <- cbind(annotations, t(coef(nmf_fit)))

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", nmf_results)

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/rank5/nmf_results_coeflabels_2806_1.xlsx")

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
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/rank5/coef_norm_2706_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", coef(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/rank5/coef_fit_2706_1.xlsx")

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, "Data")

# Write data to worksheet
writeData(wb, "Data", basis(nmf_fit))

# Save workbook
saveWorkbook(wb, "/Users/fakias0a/Documents/karla_thesis/normalized/updated_2706/results_NMF/6/rank5/basis_fit_2706_1.xlsx")

#### proving normalization

colSums(basis(nmf_fit))

rowSums(coef_norm)

plot(silhouette(nmf_fit))
si <- silhouette(nmf_fit, what = 'consensus')
plot(si)