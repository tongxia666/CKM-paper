# Load necessary libraries
options(repos = c(CRAN = "https://cloud.r-project.org"))
myPaths <- .libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs", myPaths)  # add new path
.libPaths(myPaths)  # reassign them

library(caret, lib.loc = "/n/home_fasse/txia/R/libs")

library(glmnet, lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
library(pROC, lib.loc = "/n/home_fasse/txia/R/libs")
library(parallel, lib.loc = "/n/home_fasse/txia/R/libs")
library(corrplot, lib.loc = "/n/home_fasse/txia/R/libs")





# Load the proteomics annotation file
print("Loading protein information...")
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")

# Add ".x" to the UniProt column in proteomics_info to match logistic results
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")

# Load logistic regression results for each comparison
print("Loading logistic regression results...")

logistic_results_1_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_1_vs_0_combined.csv")
logistic_results_2_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_2_vs_0_combined.csv")
logistic_results_3_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_3_vs_0_combined.csv")
logistic_results_4_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_4_vs_0_combined.csv")
logistic_results_4_vs_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_4_vs_1_combined.csv")
logistic_results_3_vs_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_3_vs_1_combined.csv")
logistic_results_2_vs_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_2_vs_1_combined.csv")
logistic_results_4_vs_2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_4_vs_2_combined.csv")
logistic_results_3_vs_2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_3_vs_2_combined.csv")
logistic_results_4_vs_3 <- fread("/n/home_fasse/txia/UKB_CKM/Data/logistic_results_4_vs_3_combined.csv")

# Combine each logistic regression result with protein information
print("Combining results with protein information...")

combined_1_vs_0 <- merge(logistic_results_1_vs_0, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_2_vs_0 <- merge(logistic_results_2_vs_0, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_3_vs_0 <- merge(logistic_results_3_vs_0, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_4_vs_0 <- merge(logistic_results_4_vs_0, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_4_vs_1 <- merge(logistic_results_4_vs_1, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_3_vs_1 <- merge(logistic_results_3_vs_1, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_2_vs_1 <- merge(logistic_results_2_vs_1, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_4_vs_2 <- merge(logistic_results_4_vs_2, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_3_vs_2 <- merge(logistic_results_3_vs_2, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")
combined_4_vs_3 <- merge(logistic_results_4_vs_3, proteomics_info, by.x = "PROTEIN", by.y = "UniProt")

# Filter proteins with M2 Bonferroni P-value < 0.05
filtered_bonf_1_vs_0 <- combined_1_vs_0 %>% filter(BONF_M2 < 0.05)
filtered_bonf_2_vs_0 <- combined_2_vs_0 %>% filter(BONF_M2 < 0.05)
filtered_bonf_3_vs_0 <- combined_3_vs_0 %>% filter(BONF_M2 < 0.05)
filtered_bonf_4_vs_0 <- combined_4_vs_0 %>% filter(BONF_M2 < 0.05)
filtered_bonf_4_vs_1 <- combined_4_vs_1 %>% filter(BONF_M2 < 0.05)
filtered_bonf_3_vs_1 <- combined_3_vs_1 %>% filter(BONF_M2 < 0.05)
filtered_bonf_2_vs_1 <- combined_2_vs_1 %>% filter(BONF_M2 < 0.05)
filtered_bonf_4_vs_2 <- combined_4_vs_2 %>% filter(BONF_M2 < 0.05)
filtered_bonf_3_vs_2 <- combined_3_vs_2 %>% filter(BONF_M2 < 0.05)
filtered_bonf_4_vs_3 <- combined_4_vs_3 %>% filter(BONF_M2 < 0.05)

# Export combined results for each comparison
output_folder <- "/n/home_fasse/txia/UKB_CKM/Data/"

write.csv(combined_1_vs_0, paste0(output_folder, "combined_protein_info_1_vs_0.csv"), row.names = FALSE)
write.csv(combined_2_vs_0, paste0(output_folder, "combined_protein_info_2_vs_0.csv"), row.names = FALSE)
write.csv(combined_3_vs_0, paste0(output_folder, "combined_protein_info_3_vs_0.csv"), row.names = FALSE)
write.csv(combined_4_vs_0, paste0(output_folder, "combined_protein_info_4_vs_0.csv"), row.names = FALSE)
write.csv(combined_4_vs_1, paste0(output_folder, "combined_protein_info_4_vs_1.csv"), row.names = FALSE)
write.csv(combined_3_vs_1, paste0(output_folder, "combined_protein_info_3_vs_1.csv"), row.names = FALSE)
write.csv(combined_2_vs_1, paste0(output_folder, "combined_protein_info_2_vs_1.csv"), row.names = FALSE)
write.csv(combined_4_vs_2, paste0(output_folder, "combined_protein_info_4_vs_2.csv"), row.names = FALSE)
write.csv(combined_3_vs_2, paste0(output_folder, "combined_protein_info_3_vs_2.csv"), row.names = FALSE)
write.csv(combined_4_vs_3, paste0(output_folder, "combined_protein_info_4_vs_3.csv"), row.names = FALSE)



# Export filtered results for Bonferroni P-value < 0.05
write.csv(filtered_bonf_1_vs_0, paste0(output_folder, "filtered_bonf_1_vs_0.csv"), row.names = FALSE)
write.csv(filtered_bonf_2_vs_0, paste0(output_folder, "filtered_bonf_2_vs_0.csv"), row.names = FALSE)
write.csv(filtered_bonf_3_vs_0, paste0(output_folder, "filtered_bonf_3_vs_0.csv"), row.names = FALSE)
write.csv(filtered_bonf_4_vs_0, paste0(output_folder, "filtered_bonf_4_vs_0.csv"), row.names = FALSE)
write.csv(filtered_bonf_4_vs_1, paste0(output_folder, "filtered_bonf_4_vs_1.csv"), row.names = FALSE)
write.csv(filtered_bonf_3_vs_1, paste0(output_folder, "filtered_bonf_3_vs_1.csv"), row.names = FALSE)
write.csv(filtered_bonf_2_vs_1, paste0(output_folder, "filtered_bonf_2_vs_1.csv"), row.names = FALSE)
write.csv(filtered_bonf_4_vs_2, paste0(output_folder, "filtered_bonf_4_vs_2.csv"), row.names = FALSE)
write.csv(filtered_bonf_3_vs_2, paste0(output_folder, "filtered_bonf_3_vs_2.csv"), row.names = FALSE)
write.csv(filtered_bonf_4_vs_3, paste0(output_folder, "filtered_bonf_4_vs_3.csv"), row.names = FALSE)


print("All files exported successfully.")


# Step 1: Data Preparation
# Load necessary libraries
options(repos = c(CRAN = "https://cloud.r-project.org"))
myPaths <- .libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs", myPaths)  # add new path
.libPaths(myPaths)  # reassign them

library(glmnet, lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
library(pROC, lib.loc = "/n/home_fasse/txia/R/libs")

# Load the dataset
data1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/combined_proteomics_pheno.csv")
data2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/covariate.csv")

# Merge datasets
ukb <- data1 %>% left_join(data2, by = "f.eid")

# Remove rows with missing 'ckm'
ukb <- ukb[!is.na(ukb$ckm), ]

# Define proteomics columns (assuming they are from column 2 to column 2924)
proteomics_columns <- colnames(ukb)[2:2924]

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
ukb[, (proteomics_columns) := lapply(.SD, function(protein) {
  protein[is.na(protein)] <- median(protein, na.rm = TRUE)
  return(protein)
}), .SDcols = proteomics_columns]

# Inverse Normal Transformation (INT)
ukb[, (proteomics_columns) := lapply(.SD, function(protein) {
  protein <- as.numeric(protein)
  INT_values <- qnorm(rank(protein) / (length(protein) + 1), mean = 0, sd = 1)
  return(INT_values)
}), .SDcols = proteomics_columns]


#ukb<-ukb %>% sample_n(5000)
# Remove rows with missing 'ckm'
ukb <- ukb[!is.na(ukb$ckm), ]




# Check for NA values in 'ckm' column
if (any(is.na(ukb$ckm))) {
  print("The 'ckm' column has NA values.")
} else {
  print("No NA values found in the 'ckm' column.")
}

#Step 2: Protein Selection Based on Previous Results
# Load filtered proteins for each comparison (Bonferroni-adjusted P-value < 0.05)
filtered_bonf_1_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_1_vs_0.csv")$PROTEIN
filtered_bonf_2_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_2_vs_0.csv")$PROTEIN
filtered_bonf_3_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_3_vs_0.csv")$PROTEIN
filtered_bonf_4_vs_0 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_4_vs_0.csv")$PROTEIN
filtered_bonf_4_vs_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_4_vs_1.csv")$PROTEIN
filtered_bonf_3_vs_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_3_vs_1.csv")$PROTEIN
filtered_bonf_2_vs_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_2_vs_1.csv")$PROTEIN
filtered_bonf_4_vs_2<- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_4_vs_2.csv")$PROTEIN
filtered_bonf_3_vs_2<- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_3_vs_2.csv")$PROTEIN
filtered_bonf_4_vs_3<- fread("/n/home_fasse/txia/UKB_CKM/Data/filtered_bonf_4_vs_3.csv")$PROTEIN

# Ensure the selected proteins exist in the combined_proteomics_pheno dataset
selected_proteins_1_vs_0 <- intersect(filtered_bonf_1_vs_0, colnames(ukb))
selected_proteins_2_vs_0 <- intersect(filtered_bonf_2_vs_0, colnames(ukb))
selected_proteins_3_vs_0 <- intersect(filtered_bonf_3_vs_0, colnames(ukb))
selected_proteins_4_vs_0 <- intersect(filtered_bonf_4_vs_0, colnames(ukb))
selected_proteins_4_vs_1 <- intersect(filtered_bonf_4_vs_1, colnames(ukb))
selected_proteins_3_vs_1 <- intersect(filtered_bonf_3_vs_1, colnames(ukb))
selected_proteins_2_vs_1 <- intersect(filtered_bonf_2_vs_1, colnames(ukb))
selected_proteins_4_vs_2 <- intersect(filtered_bonf_4_vs_2, colnames(ukb))
selected_proteins_3_vs_2 <- intersect(filtered_bonf_3_vs_2, colnames(ukb))
selected_proteins_4_vs_3 <- intersect(filtered_bonf_4_vs_3, colnames(ukb))





# Check for NA values in selected proteins for 3 vs 2 comparison
na_proteins_check <- sapply(selected_proteins_3_vs_2, function(col) any(is.na(ukb[[col]])))
if (any(na_proteins_check)) {
  print("The following selected proteins have NA values:")
  print(names(na_proteins_check[na_proteins_check]))
} else {
  print("No NA values found in the selected proteins for 3 vs 2 comparison.")
}

# Filter the dataset to include only rows where ckm is 2 or 0 for the `3 vs 2` comparison
comparison_column <- "ckm"  # Specify the column for ckm stages

print("Filtering dataset for 3 vs 2 comparison...")
ukb <- ukb[ckm %in% c(2,3), ]  # Only include rows with ckm 0 and 1
ukb$ckm <- as.factor(ukb$ckm)  # Convert ckm to a factor
ukb$ckm <- droplevels(ukb$ckm)





# Step 1: Split into training (50%), optimization (25%), testing (25%)
set.seed(123)
train_idx <- createDataPartition(ukb$ckm, p = 0.5, list = FALSE)
training_set <- ukb[train_idx, ]
temp_set <- ukb[-train_idx, ]

opt_idx <- createDataPartition(temp_set$ckm, p = 0.5, list = FALSE)
optimization_set <- temp_set[opt_idx, ]
testing_set <- temp_set[-opt_idx, ]

# Step 2: Bootstrap feature selection on training set
bootstrap_iterations <- 100
cf <- list()
models <- list()
jj <- lapply(1:bootstrap_iterations, function(x) sample(nrow(training_set), round(nrow(training_set) * 0.8)))

for (i in 1:bootstrap_iterations) {
  print(i)
  x_sample <- as.matrix(training_set[jj[[i]], ..selected_proteins_3_vs_2])
  y_sample <- factor(training_set[jj[[i]], ckm])

  models[[i]] <- cv.glmnet(
    x = x_sample,
    y = y_sample,
    alpha = 0.5,
    family = "binomial",
    nfolds = 10,
    type.measure = "auc"
  )
  lambda_min <- models[[i]]$lambda.min
  fitted_model <- glmnet(x_sample, y_sample, alpha = 0.5, lambda = lambda_min, family = "binomial")
  cf[[i]] <- as.matrix(coef(fitted_model, s = lambda_min))
  gc()
}

# Step 3: Summarize feature selection
prot.select <- matrix(NA, nrow = length(selected_proteins_3_vs_2) + 1, ncol = bootstrap_iterations)
for (i in 1:bootstrap_iterations) {
  temp <- ifelse(cf[[i]][, 1] != 0, 1, 0)
  prot.select[, i] <- temp
  rownames(prot.select) <- rownames(cf[[i]])
}

lasso.count <- as.data.frame(prot.select)
lasso.count$select <- rowSums(lasso.count)
lasso.count <- lasso.count[-1, ]
lasso.count$Protein <- rownames(lasso.count)
lasso.sort <- lasso.count[order(-lasso.count$select), ]
final_proteins <- filter(lasso.sort, select >= 0.95 * bootstrap_iterations)$Protein

# Step 4: Optimization model on selected features using optimization set
x_opti <- as.matrix(optimization_set[, ..final_proteins])
y_opti <- factor(optimization_set$ckm)
final_cv_model <- cv.glmnet(
  x = x_opti,
  y = y_opti,
  alpha = 0.5,
  family = "binomial",
  nfolds = 10,
  type.measure = "auc"
)
final_lambda <- final_cv_model$lambda.min

# Fit final model and coefficients using optimization set and selected lambda
final_model <- glmnet(
  x = x_opti,
  y = y_opti,
  alpha = 0.5,
  lambda = final_lambda,
  family = "binomial"
)


# Extract coefficients from the training model
coef_matrix <- coef(final_model, s = final_lambda,)
coef_vector <- as.vector(coef_matrix[-1])  # Exclude intercept


#########get protein score
x_training <- as.matrix(training_set[, ..final_proteins])
y_training <- factor(training_set$ckm)
training_set$protein_score <- predict(final_model, newx = x_training, s  = final_lambda, type = "response")
optimization_set$protein_score <- predict(final_model, newx = x_opti, s = final_lambda, type = "response")

# Step 5: Final model evaluation on test set
x_test <- as.matrix(testing_set[, ..final_proteins])
y_test <- factor(testing_set$ckm)


testing_set$protein_score <- predict(final_model, newx = x_test, s = final_lambda, type = "response")


#########get auc
roc_all <- roc(as.numeric(as.character(y_test)), testing_set$protein_score)
auc_value <- auc(roc_all)
print(paste("AUC:", round(auc_value, 4)))







# Merge with protein info
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")

coef_vector <- as.vector(coef(final_model)[-1])
important_proteins <- data.frame(
  Protein = final_proteins,
  Coefficient = coef_vector[final_proteins %in% colnames(x_opti)]
)

merged_result <- merge(important_proteins, proteomics_info, by.x = "Protein", by.y = "UniProt", all.x = TRUE)
write.csv(merged_result, "/n/home_fasse/txia/UKB_CKM/Data/3vs2_selected_proteins_with_info_bootstrap.csv", row.names = FALSE)








# Load protein annotation file
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")

# Thresholds for coefficient filtering
thresholds <- c(0.2)

for (threshold in thresholds) {
  # Filter proteins by coefficient magnitude from the final full model
  filtered_proteins_info <- important_proteins %>% filter(abs(Coefficient) >= threshold)
  filtered_proteins <- filtered_proteins_info$Protein
  
  # Refit model on optimization set using filtered proteins
  x_thresh <- as.matrix(optimization_set[, ..filtered_proteins])
  y_thresh <- factor(optimization_set$ckm)

  model_thresh_cv <- cv.glmnet(
    x = x_thresh,
    y = y_thresh,
    alpha = 0.5,
    family = "binomial",
    nfolds = 10,
    type.measure = "auc"
  )
  lambda_thresh <- model_thresh_cv$lambda.min

  model_thresh <- glmnet(
    x = x_thresh,
    y = y_thresh,
    alpha = 0.5,
    lambda = lambda_thresh,
    family = "binomial"
  )

  # Extract non-zero coefficients and protein names (excluding intercept)
  coef_matrix <- as.matrix(coef(model_thresh, s = lambda_thresh))[-1, , drop = FALSE]
  protein_names <- rownames(coef_matrix)
  coef_vals <- as.numeric(coef_matrix)


  # Create data frame
  df_coef <- data.frame(
    Protein = protein_names,
    Coefficient = as.numeric(coef_vals)
  )

  # Merge with protein info and save
  merged_result <- merge(df_coef, proteomics_info, by.x = "Protein", by.y = "UniProt", all.x = TRUE)
  out_path <- paste0("/n/home_fasse/txia/UKB_CKM/Data/3vs2_selected_proteins_coeffs_threshold_", threshold, ".csv")
  write.csv(merged_result, out_path, row.names = FALSE)

  # Predict protein scores on all three sets
  x_train_thresh <- as.matrix(training_set[, ..filtered_proteins])
  x_opti_thresh <- as.matrix(optimization_set[, ..filtered_proteins])
  x_test_thresh <- as.matrix(testing_set[, ..filtered_proteins])

  training_set[[paste0("protein_score_", threshold)]] <-
    predict(model_thresh, newx = x_train_thresh, s = lambda_thresh, type = "response")
  optimization_set[[paste0("protein_score_", threshold)]] <-
    predict(model_thresh, newx = x_opti_thresh, s = lambda_thresh, type = "response")
  testing_set[[paste0("protein_score_", threshold)]] <-
    predict(model_thresh, newx = x_test_thresh, s = lambda_thresh, type = "response")
}


#########get auc

roc_protein_score_0.2 <- roc(as.numeric(as.character(y_test)), testing_set$protein_score_0.2)
auc_protein_score_0.2 <- auc(roc_protein_score_0.2)
print(paste("AUC of auc_protein_score_0.2 :", round(auc_protein_score_0.2 , 4)))



# Ensure consistent types BEFORE modeling
optimization_set$sex.y <- as.factor(optimization_set$sex.y)
optimization_set$nonwhite <- as.factor(optimization_set$nonwhite)
optimization_set$alcohol_cat <- as.factor(optimization_set$alcohol_cat)

testing_set$sex.y <- as.factor(testing_set$sex.y)
testing_set$nonwhite <- as.factor(testing_set$nonwhite)
testing_set$alcohol_cat <- as.factor(testing_set$alcohol_cat)

# Define covariates for the models
covariates <- c("age.y", "sex.y", "nonwhite", "townsend","PA", "alcohol_cat")



#################add protein coefficinet 0.3 indivudlaly in model
# === Load protein annotation (if not already loaded) ===
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")

# === Step 1: Filter proteins with coefficient > 0.3 from final model ===
filtered_proteins_info_0.3 <- important_proteins %>% filter(abs(Coefficient) >= 0.3)
filtered_proteins_0.3 <- filtered_proteins_info_0.3$Protein

# === Step 2a: Model with proteins only ===
formula_protein_0.3_only <- as.formula(paste("ckm ~", paste(filtered_proteins_0.3, collapse = " + ")))
model_protein_0.3_only <- glm(formula_protein_0.3_only, data = optimization_set, family = binomial)
# Print the logistic regression summary for proteins with coef > 0.3
summary(model_protein_0.3_only)

predicted_prob_protein_0.3_only <- predict(model_protein_0.3_only, newdata = testing_set, type = "response")
roc_protein_0.3_only <- roc(as.numeric(as.character(testing_set$ckm)), predicted_prob_protein_0.3_only)
auc_protein_0.3 <- auc(roc_protein_0.3_only)
print(paste("AUC of auc_protein_0.3 :", round(auc_protein_0.3 , 4)))

# === Step 2b: Model with covariates + proteins ===
formula_covariates_protein_0.3 <- as.formula(paste("ckm ~", paste(c(covariates, filtered_proteins_0.3), collapse = " + ")))
model_covariates_protein_0.3 <- glm(formula_covariates_protein_0.3, data = optimization_set, family = binomial)
# Print the logistic regression summary for proteins with coef > 0.3
summary(model_covariates_protein_0.3)

predicted_prob_covariates_protein_0.3 <- predict(model_covariates_protein_0.3, newdata = testing_set, type = "response")
roc_covariates_protein_0.3 <- roc(as.numeric(as.character(testing_set$ckm)), predicted_prob_covariates_protein_0.3)
auc_covariates_protein_0.3 <- auc(roc_covariates_protein_0.3)
print(paste("AUC of auc_covariate_protein_0.3 :", round(auc_covariates_protein_0.3 , 4)))




# Save combined dataset
combined_data <- rbindlist(list(training_set, optimization_set, testing_set), use.names = TRUE, fill = TRUE)
write.csv(combined_data, "/n/home_fasse/txia/UKB_CKM/Data/3_vs_2_combined_training_optimization_testing.csv", row.names = FALSE)



##############test manuulay protein_score and direct generatelly protein_score
# Extract intercept and betas
coef_matrix <- coef(final_model, s = final_lambda)
intercept <- as.numeric(coef_matrix[1, 1])
coef_vector <- as.numeric(coef_matrix[-1, 1])

# Manual logit and probability
manual_logit_score <- intercept + x_test %*% coef_vector
manual_protein_score <- 1 / (1 + exp(-manual_logit_score))

# Compare
difference <- abs(manual_protein_score - testing_set$protein_score)
max_difference <- max(difference)
print(paste("Maximum absolute difference:", round(max_difference, 10)))




#################download coeffandintercept to validate in NHS
# Merge with protein info
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")

important_proteins <- data.frame(
  Protein = final_proteins,
  Coefficient = coef_vector[final_proteins %in% colnames(x_opti)]
)

merged_result <- merge(important_proteins, proteomics_info, by.x = "Protein", by.y = "UniProt", all.x = TRUE)

# Add intercept row manually
intercept_row <- data.frame(
  Protein = "(Intercept)",
  Coefficient = intercept,
  matrix(NA, nrow = 1, ncol = ncol(merged_result) - 2)  # Fill remaining columns with NA
)
colnames(intercept_row) <- colnames(merged_result)

# Combine intercept + merged proteins
final_output_table <- rbind(intercept_row, merged_result)

# Save
write.csv(final_output_table, "/n/home_fasse/txia/UKB_CKM/Data/3vs2_selected_proteins_coeffandintercept_with_info_bootstrap.csv", row.names = FALSE)

print("✅ Saved intercept + coefficients + protein annotations together!")


if (threshold == 0.2) {
  # Extract coefficients from the refit model
  coef_matrix_0.2 <- coef(model_thresh, s = lambda_thresh)

  # Intercept and coefficients
  intercept_0.2 <- as.numeric(coef_matrix_0.2[1, 1])
  coef_vector_0.2 <- as.numeric(coef_matrix_0.2[-1, 1])
  protein_names_0.2 <- rownames(coef_matrix_0.2)[-1]

  # Create data frame with coefficients
  important_proteins_0.2 <- data.frame(
    Protein = protein_names_0.2,
    Coefficient = coef_vector_0.2
  )

  # Merge with protein annotations
  merged_result_0.2 <- merge(important_proteins_0.2, proteomics_info, by.x = "Protein", by.y = "UniProt", all.x = TRUE)

  # Add intercept row
  intercept_row_0.2 <- data.frame(
    Protein = "(Intercept)",
    Coefficient = intercept_0.2,
    matrix(NA, nrow = 1, ncol = ncol(merged_result_0.2) - 2)
  )
  colnames(intercept_row_0.2) <- colnames(merged_result_0.2)

  # Combine and write output
  final_output_table_0.2 <- rbind(intercept_row_0.2, merged_result_0.2)

  write.csv(final_output_table_0.2, "/n/home_fasse/txia/UKB_CKM/Data/3vs2_selected_proteins_coeffandintercept_with_info_bootstrap_coef0.2_from_refit.csv", row.names = FALSE)
  print("✅ Saved intercept + coefficients from refit model (≥ 0.2 threshold)!")
}




################accuracy

# Model 1: Covariates only
formula_covariates <- as.formula(paste("ckm ~", paste(covariates, collapse = " + ")))
model_covariates <- glm(formula_covariates, data = optimization_set, family = binomial)
predicted_prob_covariates <- predict(model_covariates, newdata = testing_set, type = "response")
roc_covariates <- roc(as.numeric(as.character(testing_set$ckm)), predicted_prob_covariates)


# Model 2: Covariates + Proteomics score
formula_covariates_proteins <- as.formula(paste("ckm ~", paste(c(covariates, "protein_score"), collapse = " + ")))
model_covariates_proteins <- glm(formula_covariates_proteins, data = optimization_set, family = binomial)
predicted_prob_covariates_proteins <- predict(model_covariates_proteins, newdata = testing_set, type = "response")
roc_covariates_proteins <- roc(as.numeric(as.character(testing_set$ckm)), predicted_prob_covariates_proteins)



# Model 3: Covariates + Proteomics score_0.2
formula_covariates_proteins_0.2 <- as.formula(paste("ckm ~", paste(c(covariates, "protein_score_0.2"), collapse = " + ")))
model_covariates_proteins_0.2 <- glm(formula_covariates_proteins_0.2, data = optimization_set, family = binomial)
predicted_prob_covariates_proteins_0.2 <- predict(model_covariates_proteins_0.2, newdata = testing_set, type = "response")
roc_covariates_proteins_0.2 <- roc(as.numeric(as.character(testing_set$ckm)), predicted_prob_covariates_proteins_0.2)


##########draw auc curves
##########draw auc curves

# Combine all ROC objects into one named list
roc_curves <- list(
  "Proteomics Signature (All Proteins)" = roc_all,
  "Proteomics Signature (Proteins weight >= 0.2)" = roc_protein_score_0.2
  
)

plot_colors <- c("red", "blue")

legend_labels <- c()

# ---------- Basic AUC Plot ----------
pdf("/n/home_fasse/txia/UKB_CKM/Data/Combined_proteomicscore_covariate_ROC_Curves_3_vs_2.pdf", width = 18, height = 8)
first_curve_name <- names(roc_curves)[1]
plot(roc_curves[[first_curve_name]], col = plot_colors[1], main = "Combined ROC Curves", lwd = 2)

legend_labels <- c(legend_labels, paste0(
  first_curve_name, " (AUC: ", format(round(auc(roc_curves[[first_curve_name]]), 3), nsmall = 3), ")"
))

for (i in 2:length(roc_curves)) {
  curve_name <- names(roc_curves)[i]
  lines(roc_curves[[curve_name]], col = plot_colors[i], lwd = 2)
  legend_labels <- c(legend_labels, paste0(
    curve_name, " (AUC: ", format(round(auc(roc_curves[[curve_name]]), 3), nsmall = 3), ")"
  ))
}

legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2)
dev.off()

# ---------- AUC with 95% CI Plot ----------
pdf("/n/home_fasse/txia/UKB_CKM/Data/Combined_proteomicscore_covariate_ROC_Curves_with_CI_legend_3_vs_2.pdf", width = 18, height = 8)
legend_labels <- c()
first_roc <- roc_curves[[1]]
plot(first_roc, col = plot_colors[1], main = "Combined ROC Curves with 95% CI", lwd = 2)

ci_first <- ci.auc(first_roc, boot.n = 100)
legend_labels <- c(legend_labels, paste0(
  names(roc_curves)[1], " (AUC: ", format(round(auc(first_roc), 3), nsmall = 3),
  ", 95% CI: ", format(round(ci_first[1], 3), nsmall = 3), "-",
  format(round(ci_first[3], 3), nsmall = 3), ")"
))

for (i in 2:length(roc_curves)) {
  roc_curve <- roc_curves[[i]]
  lines(roc_curve, col = plot_colors[i], lwd = 2)

  ci_auc <- ci.auc(roc_curve, boot.n = 100)
  legend_labels <- c(legend_labels, paste0(
    names(roc_curves)[i], " (AUC: ", format(round(auc(roc_curve), 3), nsmall = 3),
    ", 95% CI: ", format(round(ci_auc[1], 3), nsmall = 3), "-",
    format(round(ci_auc[3], 3), nsmall = 3), ")"
  ))
}

legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2)
dev.off()

# ---------- Print AUC + CI ----------
cat("AUC Values with 95% CI:\n")
for (curve_name in names(roc_curves)) {
  this_roc <- roc_curves[[curve_name]]
  this_ci <- ci.auc(this_roc, boot.n = 100)
  cat(curve_name, ": AUC =", format(round(auc(this_roc), 3), nsmall = 3),
      " (95% CI:", format(round(this_ci[1], 3), nsmall = 3), "-",
      format(round(this_ci[3], 3), nsmall = 3), ")\n")
}









