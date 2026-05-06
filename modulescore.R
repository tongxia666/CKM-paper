# Load necessary libraries
options(repos = c(CRAN = "https://cloud.r-project.org"))
myPaths <- .libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs", myPaths)  # Add new path
.libPaths(myPaths)  # Reassign them

library(glmnet, lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
library(pROC, lib.loc = "/n/home_fasse/txia/R/libs")
library(parallel, lib.loc = "/n/home_fasse/txia/R/libs")

# Step 1: Load the proteomics annotation file
print("Loading protein information...")
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")
print("Proteomics info loaded.")
print(paste("Number of rows in proteomics_info:", nrow(proteomics_info)))
print(paste("Columns in proteomics_info:", paste(colnames(proteomics_info), collapse = ", ")))

# Add ".x" to the UniProt column in proteomics_info to match logistic results
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")
print("Updated UniProt identifiers in proteomics_info.")

# Step 2: Load the dataset
print("Loading UKB datasets...")
data1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/combined_proteomics_pheno.csv")
print("Data1 loaded.")
print(paste("Number of rows in data1:", nrow(data1)))
data2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/covariate.csv")
print("Data2 loaded.")
print(paste("Number of rows in data2:", nrow(data2)))

# Merge datasets
ukb <- data1 %>% left_join(data2, by = "f.eid")
print("Datasets merged.")
print(paste("Number of rows in ukb:", nrow(ukb)))
print("Columns in ukb:")
print(colnames(ukb))

# Remove rows with missing 'ckm'
ukb <- ukb[!is.na(ukb$ckm), ]
print(paste("Number of rows in ukb after removing missing 'ckm':", nrow(ukb)))

# Step 3: Define proteomics columns (assuming they are from column 2 to column 2924)
print("Defining proteomics columns...")
proteomics_columns <- colnames(ukb)[2:2924]
print(paste("Number of proteomics columns:", length(proteomics_columns)))



# Calculate the percentage of missing values for each proteomics variable
missing_percentage <- sapply(ukb[, ..proteomics_columns], function(x) sum(is.na(x)) / length(x) * 100)
# Print the missing percentage for each proteomics variable
cat("Missing percentage for each proteomics variable:\n")
print(missing_percentage)

# Identify proteomics variables with more than 80% missing ukb
high_missing_proteins <- names(missing_percentage[missing_percentage > 80])

# Print the number of variables with more than 80% missing ukb
cat("Number of proteomics variables with more than 80% missing ukb:", length(high_missing_proteins), "\n")

# Exclude the proteomics variable with more than 80% missing ukb
proteomics_columns <- proteomics_columns[!proteomics_columns %in% "P48060.x"]



# Median imputation for missing protein values
print("Performing median imputation for missing protein values...")
ukb[, (proteomics_columns) := lapply(.SD, function(protein) {
  protein[is.na(protein)] <- median(protein, na.rm = TRUE)
  return(protein)
}), .SDcols = proteomics_columns]
print("Median imputation completed.")

# Inverse Normal Transformation (INT)
print("Performing Inverse Normal Transformation (INT)...")
ukb[, (proteomics_columns) := lapply(.SD, function(protein) {
  protein <- as.numeric(protein)
  INT_values <- qnorm(rank(protein, ties.method = "average") / (length(protein) + 1), mean = 0, sd = 1)
  return(INT_values)
}), .SDcols = proteomics_columns]
print("INT completed.")

# Remove rows with missing 'ckm' (in case any were introduced)
ukb <- ukb[!is.na(ukb$ckm), ]
print(paste("Number of rows in ukb after re-removing missing 'ckm':", nrow(ukb)))

# Check for NA values in 'ckm' column
if (any(is.na(ukb$ckm))) {
  print("The 'ckm' column has NA values.")
} else {
  print("No NA values found in the 'ckm' column.")
}

# Step 4: Filter UKB Dataset to Include Only Relevant Proteins
print("Loading weighted sums data...")
weighted_sums <- fread("/n/home_fasse/txia/UKB_CKM/Data/weighted_sums.csv")
print("Weighted sums data loaded.")
print(paste("Number of rows in weighted_sums:", nrow(weighted_sums)))
print("Columns in weighted_sums:")
print(colnames(weighted_sums))

# Match `p_Assay` in Weighted Sums with `p_Assay` in Proteomics Info and Retrieve `UniProt`
print("Matching p_Assay in weighted_sums with proteomics_info...")
matched_weights <- weighted_sums %>%
  inner_join(proteomics_info, by = "p_Assay")  # Match by `p_Assay`
print("Matching completed.")
print(paste("Number of matched proteins:", nrow(matched_weights)))

if (nrow(matched_weights) == 0) {
  stop("No matches found between weighted_sums and proteomics_info. Please check 'p_Assay' identifiers.")
}

# Get relevant protein names (UniProt) from the matched data
relevant_proteins <- matched_weights$UniProt
print(paste("Number of relevant proteins:", length(relevant_proteins)))

# Check if relevant proteins are in ukb
proteins_in_ukb <- relevant_proteins %in% colnames(ukb)
if (!all(proteins_in_ukb)) {
  missing_proteins <- relevant_proteins[!proteins_in_ukb]
  print("The following proteins are missing in ukb:")
  print(missing_proteins)
  # Optionally, remove missing proteins
  relevant_proteins <- relevant_proteins[proteins_in_ukb]
  matched_weights <- matched_weights[matched_weights$UniProt %in% relevant_proteins, ]
  print(paste("Proceeding with", length(relevant_proteins), "proteins that are present in ukb."))
}

# Filter UKB dataset to include only relevant proteins
print("Filtering UKB dataset for relevant proteins...")
ukb_filtered <- ukb %>%
  select(f.eid, all_of(relevant_proteins))  # Select only f.eid and relevant proteins
print("UKB dataset filtered.")
print(paste("Number of columns in ukb_filtered:", ncol(ukb_filtered)))
print("Columns in ukb_filtered:")
print(colnames(ukb_filtered))

# Step 5: Calculate Module Scores Using Matched Weights
print("Calculating module scores...")
module_scores <- data.table(f.eid = ukb_filtered$f.eid)

# Loop through each module column in weighted_sums
module_columns <- colnames(weighted_sums)[-1]  # Exclude 'p_Assay'
print(paste("Number of modules:", length(module_columns)))
print("Module columns:")
print(module_columns)

for (module_col in module_columns) {
  print(paste("Processing module:", module_col))
  # Get weights for the current module
  module_weights <- matched_weights[[module_col]]
  
  # Check for NA weights
  if (all(is.na(module_weights))) {
    print(paste("All weights are NA for module", module_col, ". Skipping this module."))
    next
  }
  
  # Use only proteins with non-zero weights
  non_zero_weights <- module_weights != 0
  protein_subset <- matched_weights$UniProt[non_zero_weights]
  weights_subset <- module_weights[non_zero_weights]
  
  if (length(protein_subset) == 0) {
    print(paste("No proteins with non-zero weights for module", module_col, ". Skipping this module."))
    next
  }
  
  # Subset the UKB filtered dataset to include only proteins in this module
  ukb_subset <- ukb_filtered[, ..protein_subset]
  
  # Ensure the weights and proteins align
  names(weights_subset) <- protein_subset
  
  # Calculate weighted score for the module
  module_scores[[module_col]] <- as.vector(as.matrix(ukb_subset) %*% weights_subset)
  print(paste("Module score calculated for", module_col))
}

print("Module scores calculation completed.")
print("Module scores dataframe:")
print(head(module_scores))





# Step 6: Save the Module Scores (unstandardized)
print("Saving module scores to CSV...")
fwrite(module_scores, "/n/home_fasse/txia/UKB_CKM/Data/module_scores.csv")
print("Module scores saved successfully.")

# Step 7: Standardize module scores (excluding 'f.eid')
module_names <- colnames(module_scores)[-1]  # All except f.eid
module_scores[, (module_names) := lapply(.SD, scale), .SDcols = module_names]


# Step 8: Merge standardized module scores with UKB dataset
ukb_with_scores <- ukb %>% left_join(module_scores, by = "f.eid")
print("UKB dataset merged with standardized module scores.")
print(paste("Number of rows in ukb_with_scores:", nrow(ukb_with_scores)))


#ukb_with_scores<-ukb_with_scores%>%  filter(bptreat.y == 0 & statin.y == 0)


# Step 9: Create CKM binary comparison variables
print("Preparing CKM binary comparisons...")
create_binary_comparison <- function(ckm, comparison_category1, comparison_category2) {
  return(ifelse(ckm == comparison_category1, 1, 
                ifelse(ckm == comparison_category2, 0, NA)))
}

ukb_with_scores$ckm_1_vs_0 <- create_binary_comparison(ukb_with_scores$ckm, 1, 0)
ukb_with_scores$ckm_2_vs_0 <- create_binary_comparison(ukb_with_scores$ckm, 2, 0)
ukb_with_scores$ckm_2_vs_1 <- create_binary_comparison(ukb_with_scores$ckm, 2, 1)
ukb_with_scores$ckm_3_vs_0 <- create_binary_comparison(ukb_with_scores$ckm, 3, 0)
ukb_with_scores$ckm_3_vs_1 <- create_binary_comparison(ukb_with_scores$ckm, 3, 1)
ukb_with_scores$ckm_3_vs_2 <- create_binary_comparison(ukb_with_scores$ckm, 3, 2)
ukb_with_scores$ckm_4_vs_0 <- create_binary_comparison(ukb_with_scores$ckm, 4, 0)
ukb_with_scores$ckm_4_vs_1 <- create_binary_comparison(ukb_with_scores$ckm, 4, 1)
ukb_with_scores$ckm_4_vs_2 <- create_binary_comparison(ukb_with_scores$ckm, 4, 2)
ukb_with_scores$ckm_4_vs_3 <- create_binary_comparison(ukb_with_scores$ckm, 4, 3)

# Step 10: Convert relevant variables to factors
ukb_with_scores$sex.y <- as.factor(ukb_with_scores$sex.y)
ukb_with_scores$nonwhite <- as.factor(ukb_with_scores$nonwhite)
ukb_with_scores$fast_hrs <- as.factor(ukb_with_scores$fast_hrs)

ukb_with_scores$alcohol_cat <- as.factor(ukb_with_scores$alcohol_cat)
ukb_with_scores$bptreat.y<-as.factor(ukb_with_scores$bptreat.y)
ukb_with_scores$statin.y<-as.factor(ukb_with_scores$statin.y)
# Step 11: Define comparisons
binary_comparisons <- c("ckm_1_vs_0", "ckm_2_vs_0", "ckm_2_vs_1", 
                        "ckm_3_vs_0", "ckm_3_vs_1", "ckm_3_vs_2",
                        "ckm_4_vs_0", "ckm_4_vs_1", "ckm_4_vs_2", "ckm_4_vs_3")

# Step 12: Initialize result storage
result_list <- list()  # Significant results
all_results_list <- list()  # All results

# Step 13: Logistic regression for each comparison
for (comparison in binary_comparisons) {
  print(paste("Processing comparison:", comparison))
  
  result <- data.frame(matrix(NA, nrow = length(module_names), ncol = 16,
                              dimnames = list(NULL, c("Module", 
                                                      "BETA_M1", "SE_M1", "PVAL_M1", "BONF_M1", "FDR_M1",
                                                      "BETA_M2", "SE_M2", "PVAL_M2", "BONF_M2", "FDR_M2",
                                                      "BETA_M3", "SE_M3", "PVAL_M3", "BONF_M3", "FDR_M3"))))
  
  for (i in seq_along(module_names)) {
    module <- module_names[i]
    
    # Model 1
    form1 <- as.formula(paste(comparison, "~", module, "+ age.y + sex.y + nonwhite +fast_hrs"))
    model1 <- glm(form1, data = ukb_with_scores, family = binomial)
    result$Module[i] <- module
    result$BETA_M1[i] <- coef(model1)[2]
    result$SE_M1[i] <- sqrt(diag(vcov(model1)))[2]
    result$PVAL_M1[i] <- summary(model1)$coef[2, 4]
    
    # Model 2
    form2 <- as.formula(paste(comparison, "~", module, "+ age.y + sex.y + nonwhite +fast_hrs + townsend"))
    model2 <- glm(form2, data = ukb_with_scores, family = binomial)
    result$BETA_M2[i] <- coef(model2)[2]
    result$SE_M2[i] <- sqrt(diag(vcov(model2)))[2]
    result$PVAL_M2[i] <- summary(model2)$coef[2, 4]
    
    # Model 3
    form3 <- as.formula(paste(comparison, "~", module, "+ age.y + sex.y + nonwhite +fast_hrs + townsend + PA + alcohol_cat"))
    model3 <- glm(form3, data = ukb_with_scores, family = binomial)
    result$BETA_M3[i] <- coef(model3)[2]
    result$SE_M3[i] <- sqrt(diag(vcov(model3)))[2]
    result$PVAL_M3[i] <- summary(model3)$coef[2, 4]
  }
  
  # Step 14: Apply multiple testing correction
  result$BONF_M1 <- p.adjust(result$PVAL_M1, method = "bonferroni", n = length(module_names))
  result$BONF_M2 <- p.adjust(result$PVAL_M2, method = "bonferroni", n = length(module_names))
  result$BONF_M3 <- p.adjust(result$PVAL_M3, method = "bonferroni", n = length(module_names))
  
  result$FDR_M1 <- p.adjust(result$PVAL_M1, method = "fdr", n = length(module_names))
  result$FDR_M2 <- p.adjust(result$PVAL_M2, method = "fdr", n = length(module_names))
  result$FDR_M3 <- p.adjust(result$PVAL_M3, method = "fdr", n = length(module_names))
  
  # Step 15: Save all results
  fwrite(result, paste0("/n/home_fasse/txia/UKB_CKM/Data/", comparison, "_module_associations.csv"))
  print(paste("Results saved for comparison:", comparison))
  
  # Track all results
  result$Comparison <- comparison
  all_results_list[[comparison]] <- result
  
  # Save significant results (Model 3, Bonferroni < 0.05)
  sig <- result %>% filter(BONF_M3 < 0.05)
  if (nrow(sig) > 0) {
    sig$Comparison <- comparison
    result_list[[comparison]] <- sig
  }
}

# Step 16: Save combined results
final_significant_results <- do.call(rbind, result_list)
fwrite(final_significant_results, "/n/home_fasse/txia/UKB_CKM/Data/significant_modules_model3_bonferroni.csv")

combined_all_results <- do.call(rbind, all_results_list)
fwrite(combined_all_results, "/n/home_fasse/txia/UKB_CKM/Data/combined_all_module_associations.csv")

print("Logistic regression completed. All and significant results saved.")


