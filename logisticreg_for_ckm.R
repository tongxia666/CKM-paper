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
all_proteins <- all_proteins[!all_proteins %in% "P48060.x"]


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

ukb$sex.y <- as.factor(ukb$sex.y)
ukb$nonwhite <- as.factor(ukb$nonwhite)
ukb$fast_hrs <- as.factor(ukb$fast_hrs)
ukb$multivit <- as.factor(ukb$multivit)
ukb$alcohol_cat<-as.factor(ukb$alcohol_cat)
#length(all_proteins)
#print(all_proteins)
# Check missing values in 'ckm'
cat("Missing values in 'ckm':", sum(is.na(ukb$ckm)), "\n")

# Check missing values in 'age.y'
cat("Missing values in 'age.y':", sum(is.na(ukb$age.y)), "\n")

# Check missing values in 'sex.y'
cat("Missing values in 'sex.y':", sum(is.na(ukb$sex.y)), "\n")

# Check missing values in 'nonwhite'
cat("Missing values in 'nonwhite':", sum(is.na(ukb$nonwhite)), "\n")

# Check missing values in 'fast_hrs'
cat("Missing values in 'fast_hrs':", sum(is.na(ukb$fast_hrs)), "\n")

# Check missing values in 'townsend'
cat("Missing values in 'townsend':", sum(is.na(ukb$townsend)), "\n")

# Check missing values in 'energy_mean'
cat("Missing values in 'energy_mean':", sum(is.na(ukb$energy_mean)), "\n")

# Check missing values in 'multivit'
cat("Missing values in 'multivit':", sum(is.na(ukb$multivit)), "\n")

# Check missing values in 'amed'
cat("Missing values in 'amed':", sum(is.na(ukb$amed)), "\n")

# Check missing values in 'PA'
cat("Missing values in 'PA':", sum(is.na(ukb$PA)), "\n")

# Check missing values in 'alcohol'
cat("Missing values in 'alcohol':", sum(is.na(ukb$alcohol_cat)), "\n")

#ukb<- ukb[!is.na(ukb$ckm) & !is.na(ukb$age.y)  &!is.na(ukb$sex.y) &!is.na(ukb$nonwhite)  &!is.na(ukb$townsend)  &!is.na(ukb$energy_mean)  
# &!is.na(ukb$multivit)  &!is.na(ukb$amed), ]
dim(ukb)

# Convert 'ckm' to a factor if not already done
ukb$ckm <- factor(ukb$ckm, ordered = TRUE)

# Function to create binary comparison for logistic regression
#create_binary_comparison <- function(ckm, reference_category) {
#  return(ifelse(ckm == reference_category, 0, 1))
#}

# Set seed for reproducibility
set.seed(123)

# Shuffle rows
rows <- sample(seq_len(nrow(ukb)))
ukb <- ukb[rows, ]

# Define training and testing split
train_size <- round(0.75 * nrow(ukb))
training <- ukb[1:train_size, ]
testing <- ukb[(train_size + 1):nrow(ukb), ]

# Function to create binary comparison for logistic regression
create_binary_comparison <- function(ckm, include_category, reference_category) {
  return(ifelse(ckm == include_category, 1, ifelse(ckm == reference_category, 0, NA)))
}

# Apply binary comparison only in the training dataset
training$ckm_1_vs_0 <- create_binary_comparison(training$ckm, 1, 0)
training$ckm_2_vs_0 <- create_binary_comparison(training$ckm, 2, 0)
training$ckm_3_vs_0 <- create_binary_comparison(training$ckm, 3, 0)
training$ckm_4_vs_0 <- create_binary_comparison(training$ckm, 4, 0)
training$ckm_4_vs_1 <- create_binary_comparison(training$ckm, 4, 1)
training$ckm_3_vs_1 <- create_binary_comparison(training$ckm, 3, 1)
training$ckm_2_vs_1 <- create_binary_comparison(training$ckm, 2, 1)
training$ckm_4_vs_2 <- create_binary_comparison(training$ckm, 4, 2)
training$ckm_3_vs_2 <- create_binary_comparison(training$ckm, 3, 2)
training$ckm_4_vs_3 <- create_binary_comparison(training$ckm, 4, 3)

# List of binary comparisons and corresponding output file names
binary_comparisons <- list(
  "ckm_1_vs_0" = "logistic_results_1_vs_0_combined.csv",
  "ckm_2_vs_0" = "logistic_results_2_vs_0_combined.csv",
  "ckm_3_vs_0" = "logistic_results_3_vs_0_combined.csv",
  "ckm_4_vs_0" = "logistic_results_4_vs_0_combined.csv",
  "ckm_4_vs_1" = "logistic_results_4_vs_1_combined.csv",
  "ckm_3_vs_1" = "logistic_results_3_vs_1_combined.csv",
  "ckm_2_vs_1" = "logistic_results_2_vs_1_combined.csv",
  "ckm_4_vs_2" = "logistic_results_4_vs_2_combined.csv",
  "ckm_3_vs_2" = "logistic_results_3_vs_2_combined.csv",
  "ckm_4_vs_3" = "logistic_results_4_vs_3_combined.csv"
)

# Set the effective number of tests for FDR and Bonferroni corrections
meff <- 2652

# Loop over each binary comparison
for (comparison in names(binary_comparisons)) {
  
  # Define the outcome variable for the current comparison
  outcome_var <- comparison

  # Initialize results dataframe for storing the combined results of three models
  result <- data.frame(matrix(NA, length(all_proteins), 16, dimnames = list(NULL, c(
    "PROTEIN",
    "BETA_M1", "SE_M1", "PVAL_M1", "FDR_M1", "BONF_M1",
    "BETA_M2", "SE_M2", "PVAL_M2", "FDR_M2", "BONF_M2",
    "BETA_M3", "SE_M3", "PVAL_M3", "FDR_M3", "BONF_M3"
  ))))

  # Loop over each protein
  for (i in 1:length(all_proteins)) {

    # Get the protein column name
    protein_name <- all_proteins[i]

    # Print the protein being processed
    cat("Processing protein:", protein_name, "for comparison:", comparison, "\n")

    # Check if the protein column exists in the training dataset
    if (!(protein_name %in% names(training))) {
      cat("Protein column not found for:", protein_name, "\n")
      next  # Skip this iteration if the protein column is not found
    }

    # Store the protein name
    result$PROTEIN[i] <- protein_name

    # Fit the models for the current protein
    for (m in 1:3) {
      if (m == 1) {
        form <- as.formula(paste(outcome_var, "~ `", protein_name, "` + age.y + sex.y + nonwhite+fast_hrs", sep = ""))
        column_names <- c(outcome_var, "age.y", "sex.y", "nonwhite", "fast_hrs", protein_name)
      } else if (m == 2) {
        form <- as.formula(paste(outcome_var, "~ `", protein_name, "` + age.y + sex.y + nonwhite +fast_hrs+ townsend+ PA + alcohol_cat", sep = ""))
        column_names <- c(outcome_var, "age.y", "sex.y", "nonwhite", "fast_hrs", "townsend","PA", "alcohol_cat",protein_name)
      } else if (m == 3) {
        form <- as.formula(paste(outcome_var, "~ `", protein_name, "` + age.y + sex.y + nonwhite +fast_hrs+ townsend + energy_mean + multivit + amed", sep = ""))
        column_names <- c(outcome_var, "age.y", "sex.y", "nonwhite", "fast_hrs", "townsend", "energy_mean", "multivit", "amed", protein_name)
      }

      # Subset the training data with required columns
      data_subset <- training[, ..column_names]

      # Remove missing values
      data_subset <- na.omit(data_subset)

      # Skip iteration if no rows remain after removing NAs
      if (nrow(data_subset) == 0) {
        cat("No data left after NA removal for protein:", protein_name, "in model", m, "\n")
        next
      }

      # Fit logistic regression model
      model <- glm(form, data = as.data.frame(data_subset), family = binomial)

      # Store results for the current model
      if (m == 1) {
        result$BETA_M1[i] <- coef(model)[2]
        result$SE_M1[i] <- sqrt(diag(vcov(model)))[2]
        result$PVAL_M1[i] <- summary(model)$coef[2, 4]
        result$FDR_M1[i] <- p.adjust(summary(model)$coef[2, 4], method = "fdr", n = meff)
        result$BONF_M1[i] <- p.adjust(summary(model)$coef[2, 4], method = "bonferroni", n = meff)
      } else if (m == 2) {
        result$BETA_M2[i] <- coef(model)[2]
        result$SE_M2[i] <- sqrt(diag(vcov(model)))[2]
        result$PVAL_M2[i] <- summary(model)$coef[2, 4]
        result$FDR_M2[i] <- p.adjust(summary(model)$coef[2, 4], method = "fdr", n = meff)
        result$BONF_M2[i] <- p.adjust(summary(model)$coef[2, 4], method = "bonferroni", n = meff)
      } else if (m == 3) {
        result$BETA_M3[i] <- coef(model)[2]
        result$SE_M3[i] <- sqrt(diag(vcov(model)))[2]
        result$PVAL_M3[i] <- summary(model)$coef[2, 4]
        result$FDR_M3[i] <- p.adjust(summary(model)$coef[2, 4], method = "fdr", n = meff)
        result$BONF_M3[i] <- p.adjust(summary(model)$coef[2, 4], method = "bonferroni", n = meff)
      }
    }
  }

  # Save the results to a file
  output_file <- binary_comparisons[[comparison]]
  write.table(result, file = paste0("/n/home_fasse/txia/UKB_CKM/Data/", output_file), quote = F, col.names = T, row.names = F, sep = ",")

  # Print completion message
  cat("Results saved for comparison:", comparison, "to file:", output_file, "\n")
}


