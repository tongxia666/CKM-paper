options(repos = c(CRAN = "https://cloud.r-project.org"))

myPaths<-.libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs",myPaths)  # add new path
.libPaths(myPaths)  # reassign them
.libPaths()

library(brant, lib.loc = "/n/home_fasse/txia/R/libs")

#remove.packages("MASS", lib = "/n/home_fasse/txia/R/libs")
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr", lib = "/n/home_fasse/txia/R/libs", repos = "https://cloud.r-project.org")
}

library(glmnet, lib.loc = "/n/home_fasse/txia/R/libs")
library(ordinalNet, lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr, lib.loc = "/n/home_fasse/txia/R/libs")
library(matrixStats, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
library(tidyr, lib.loc = "/n/home_fasse/txia/R/libs")
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



table(ukb$ckm)
print(colnames(ukb))

# Define continuous variables
continuous_vars <- c("age.y", "townsend", "PA", "BMI", 
                     "f.48.0.0.y", "f.30750.0.0.y", "f.30740.0.0.y",
                     "f.30870.0.0.y", "f.30760.0.0.y", "f.30690.0.0.y",
                     "sbp", "dbp", "serum_creatinine_mg_dl", "eGFR")

categorical_vars <- c("sex.y", "nonwhite","fast_hrs", "alcohol_cat", "smoke", "bptreat.y", "statin.y")

# Function to compute ANOVA P-value with error handling
get_anova_p <- function(var) {
  var_data <- ukb[[var]]
  if (sum(!is.na(var_data)) < 2) return(NA)  # Ensure enough data
  aov_result <- tryCatch(aov(var_data ~ as.factor(ukb$ckm)), error = function(e) return(NA))
  if (inherits(aov_result, "try-error") || is.na(aov_result)) return(NA)  # Catch errors
  p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
  return(p_value)
}

### **1. Continuous Variables Summary (Overall & CKM-Stratified)**
# Compute summary statistics per CKM stage
continuous_summary_ckm <- ukb %>%
  group_by(ckm) %>%
  summarize(across(all_of(continuous_vars), 
                   list(mean = ~mean(. , na.rm = TRUE), 
                        sd = ~sd(. , na.rm = TRUE)), 
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(-ckm, names_to = c("Variable", ".value"), names_pattern = "(.+)_(.+)") %>%
  as.data.frame()

# Compute overall summary statistics (all participants)
continuous_summary_overall <- ukb %>%
  summarize(across(all_of(continuous_vars), 
                   list(mean = ~mean(. , na.rm = TRUE), 
                        sd = ~sd(. , na.rm = TRUE)), 
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything(), names_to = c("Variable", ".value"), names_pattern = "(.+)_(.+)") %>%
  as.data.frame()

continuous_summary_overall$ckm <- "Overall"

# Compute and add ANOVA P-values
p_values <- sapply(continuous_vars, get_anova_p)
p_value_df <- data.frame(Variable = continuous_vars, P_Value = p_values)

# Merge P-values correctly
continuous_summary_ckm <- left_join(continuous_summary_ckm, p_value_df, by = "Variable")

# Combine overall & CKM-stratified results
continuous_summary <- bind_rows(continuous_summary_overall, continuous_summary_ckm)

# Save results
write.csv(continuous_summary, "/n/home_fasse/txia/UKB_CKM/Data/continuous_summary_overall_and_by_CKM.csv", row.names = FALSE)

print(continuous_summary)

### **2. Categorical Variables Summary (Overall & CKM-Stratified)**
cat_summary <- lapply(categorical_vars, function(var) {
  table_data <- table(ukb[[var]], ukb$ckm)
  percent_data <- prop.table(table_data, margin = 2) * 100  # Convert to percentage within CKM
  
  # Use Fisher’s exact test for small sample sizes
  if (min(table_data) < 5) {
    chi_sq_test <- fisher.test(table_data)
  } else {
    chi_sq_test <- chisq.test(table_data)
  }
  
  df_ckm <- as.data.frame(table_data)
  df_ckm$Percentage <- round(as.numeric(percent_data[cbind(df_ckm$Var1, df_ckm$Var2)]), 1)
  df_ckm$P_Value <- chi_sq_test$p.value
  df_ckm$Category <- var
  colnames(df_ckm) <- c("Value", "CKM_Stage", "Count", "Percentage", "P_Value", "Variable")
  
  # Compute overall statistics (all participants)
  overall_table <- table(ukb[[var]])
  overall_percent <- prop.table(overall_table) * 100
  df_overall <- as.data.frame(overall_table)
  df_overall$Percentage <- round(as.numeric(overall_percent[df_overall$Var1]), 1)
  df_overall$P_Value <- NA  # No p-value for overall
  df_overall$CKM_Stage <- "Overall"
  df_overall$Category <- var
  colnames(df_overall) <- c("Value", "Count", "Percentage", "P_Value", "CKM_Stage", "Variable")
  
  # Combine overall & CKM-stratified results
  df_combined <- bind_rows(df_overall, df_ckm)
  
  return(df_combined)
})

# Combine categorical results into one dataframe
categorical_summary <- do.call(rbind, cat_summary)

# Save results
write.csv(categorical_summary, "/n/home_fasse/txia/UKB_CKM/Data/categorical_summary_overall_and_by_CKM.csv", row.names = FALSE)

print(categorical_summary)


