options(repos = c(CRAN = "https://cloud.r-project.org"))

myPaths<-.libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs",myPaths)  # add new path
.libPaths(myPaths)  # reassign them
.libPaths()

library(brant, lib.loc = "/n/home_fasse/txia/R/libs")

#remove.packages("MASS", lib = "/n/home_fasse/txia/R/libs")
if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS", lib = "/n/home_fasse/txia/R/libs", repos = "https://cloud.r-project.org")
}

library(glmnet, lib.loc = "/n/home_fasse/txia/R/libs")
library(ordinalNet, lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr, lib.loc = "/n/home_fasse/txia/R/libs")
library(matrixStats, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
#library(MASS, lib.loc = "/n/home_fasse/txia/R/libs")
library(MASS)

# Load the dataset
data1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/combined_proteomics_pheno.csv")
length(data1$f.eid)
data2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/covariate.csv")
length(data2$f.eid)
dim(data1)
dim(data2)

#Proteomics annotation file
proteomics_info<-read.table("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv",header=TRUE, sep=",")


ukb <- data1 |> left_join(data2, by = "f.eid")
#remove ukb missing
ukb <- ukb[complete.cases(ukb$ckm)]
length(ukb$f.eid)
dim(ukb)

#set.seed(123)  # Set seed for reproducibility
#sampled_participants <- sample(1:nrow(ukb), 500, replace = FALSE)
#ukb <- ukb[sampled_participants, ]

all_proteins <- colnames(ukb[, 2:2924])
length(all_proteins)
#print(all_proteins)
#head(all_proteins)

# Convert 'ckm' to an ordered factor if it's ordinal
ukb$ckm <- factor(ukb$ckm, ordered = TRUE)

# Calculate the percentage of missing values for each proteomics variable
missing_percentage <- sapply(ukb[, ..all_proteins], function(x) sum(is.na(x)) / length(x) * 100)
# Print the missing percentage for each proteomics variable
cat("Missing percentage for each proteomics variable:\n")
print(missing_percentage)

# Identify proteomics variables with more than 80% missing ukb
high_missing_proteins <- names(missing_percentage[missing_percentage > 80])

# Print the number of variables with more than 80% missing ukb
cat("Number of proteomics variables with more than 80% missing ukb:", length(high_missing_proteins), "\n")

# Exclude the proteomics variable with more than 80% missing ukb
#all_proteins <- all_proteins[!all_proteins %in% "P48060.x"]


# Mean imputation for proteomics ukb
#for (col in all_proteins) {
#  if (is.numeric(ukb[[col]])) {
#    ukb[[col]][is.na(ukb[[col]])] <- mean(ukb[[col]], na.rm = TRUE)
#  }
#}


# Step 1: Daily Median Correction
#ukb[, (all_proteins) := lapply(.SD, function(protein) {
#  ukb$day_corrected <- ave(protein, ukb$day, FUN = function(x) x - median(x, na.rm = TRUE))
#  return(ukb$day_corrected)
#}), .SDcols = all_proteins]

# Step 2: Median Imputation for Missing Values
ukb[, (all_proteins) := lapply(.SD, function(protein) {
  # Impute missing values with the median of each protein
  protein[is.na(protein)] <- median(protein, na.rm = TRUE)
  return(protein)
}), .SDcols = all_proteins]

# Step 4: Inverse Normal Transformation (INT)
ukb[, (all_proteins) := lapply(.SD, function(protein) {
  # Convert the protein values to numeric if necessary
  protein <- as.numeric(protein)
  
  # Apply inverse normal transformation (INT) as per previous method
  INT_values <- qnorm(rank(protein) / (length(protein) + 1), mean = 0, sd = 1)
  
  return(INT_values)
}), .SDcols = all_proteins]

# Print the first few rows to verify the transformations
print(head(ukb[, ..all_proteins]))

# Calculate the effective number of tests (Meff) using PCA on the final processed proteomics data
# Make sure all transformations are complete before running this code.

# Select only proteomics columns (update column selection if necessary)
proteomics_data <- ukb[, ..all_proteins]

# Calculate the correlation matrix on the processed data
cor_matrix <- cor(proteomics_data, use = "pairwise.complete.obs")

# Obtain eigenvalues from the correlation matrix
eigen_values <- eigen(cor_matrix)$values

# Define the cutoff for variance explained (e.g., 99.5%)
percent_cut <- 0.995
total_variance <- sum(eigen_values)
cumulative_variance <- 0
meff <- 0

# Calculate Meff by summing eigenvalues until reaching the cumulative variance cutoff
for (i in seq_along(eigen_values)) {
  cumulative_variance <- cumulative_variance + eigen_values[i]
  if (cumulative_variance / total_variance >= percent_cut) {
    meff <- i
    break
  }
}

cat("The effective number of tests (Meff) is:", meff, "\n")

