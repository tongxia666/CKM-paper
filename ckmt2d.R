

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
library(survival, lib.loc = "/n/home_fasse/txia/R/libs")

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



##################################################

# Function to select UK Biobank fields
f_select_UKB_field <- function(x, file_output, version = "July_2023") {
  cat("\nCurrent version: ", version, "\n")
  f.number = paste0("f.", x, sep = "\\.", collapse = "|")
  
  if (version == "July_2023") {
    df_head <- paste0("/n/home_fasse/txia/UKB_CKM/Data/", version, "_col.txt")
    if (!file.exists(df_head)) {
      print("Header file does not exist. Creating it...")
      df <- fread(paste0("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/July_2023_full_data/ukb674044.tab"), nrows = 1, header = FALSE)
      df <- data.table(t(df), keep.rownames = TRUE)
      df[, rn := substr(rn, 2, 100)]
      fwrite(df, df_head, sep = "\t", col.names = FALSE)
    }
    cmd <- paste0("egrep   \'", f.number, "\' ", df_head, " > ~/temp.txt")
  }
  print(cmd)
  system(cmd)
  
  temp <- fread("~/temp.txt", header = FALSE, sep = "\t")
  
  if (version == "July_2023") {
    cmd <- paste0("cut -f ", paste(c(1, temp$V1), collapse = ","), 
                  " /n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/July_2023_full_data/ukb674044.tab > ", file_output)
  }
  print(cmd)
  system(cmd)
  system("rm ~/temp.txt")
}

# Specify input field IDs and output file paths
x <- c(
  # Add your required f.xxx field IDs here
  130708, 131298, 131300, 131302, 131304, 131306, 131360, 131362, 131366, 131368, 131354, 
131364,  131296, 131058, 131056, 131386, 131350, 132032, 21842, 40000, 22189, "f.eid", 191, 42000, 42006, 42008, 42010)
file_output <- "/n/home_fasse/txia/UKB_CKM/Data/t2dvariable"
csv_output <- "/n/home_fasse/txia/UKB_CKM/Data/t2dvariable.csv"

# Apply function to select fields and save the output
f_select_UKB_field(x, file_output, version = "July_2023")

# Read the selected data into R
ukb_selected_data <- fread(file_output)

# Write the selected data to a CSV file
write.table(ukb_selected_data, file = csv_output, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")

# Load the dataset
t2d_association <- ukb
colnames(t2d_association)


ukb_data <- fread("/n/home_fasse/txia/UKB_CKM/Data/t2dvariable.csv")

# Join the datasets by participant ID (assuming 'f.eid' is the ID column in both)
t2d_association <- t2d_association %>%
  left_join(ukb_data, by = "f.eid")

# Define t2d-related outcomes (exclude CKD-related outcomes)
t2d_outcomes <- c(
  "f.130708.0.0"
)

# Convert date fields to numeric (dates are assumed to be in YYYY-MM-DD format)
for (outcome in t2d_outcomes) {
  t2d_association[[outcome]] <- as.numeric(as.Date(t2d_association[[outcome]]))
}

# Define baseline visit date (assumed to be stored in `f.21842`)
t2d_association$dt_visit0 <- as.numeric(as.Date(t2d_association$f.21842.0.0))

# Step 1: Identify first t2d diagnosis date and handle Inf values
t2d_association$t2d_dxdt <- apply(
  as.matrix(t2d_association[, ..t2d_outcomes]), 
  1, 
  function(x) {
    min_val <- min(x, na.rm = TRUE)  # Find the earliest non-NA diagnosis date
    if (is.infinite(min_val)) NA else min_val  # Replace Inf with NA
  }
)

# Replace NA with Inf where required explicitly for handling missing cases
t2d_association$t2d_dxdt[is.na(t2d_association$t2d_dxdt)] <- Inf

# Step 2: Determine baseline t2d (prevalent cases: diagnosis before baseline visit)
t2d_association$t2d_bsln <- ifelse(
  t2d_association$t2d_dxdt < t2d_association$dt_visit0,
  1, 0
)

# Step 3: Exclude participants with baseline t2d
incident_data <- t2d_association %>%
  filter(t2d_bsln == 0)

# Step 4: Define incident t2d and censoring
incident_data <- incident_data %>%
  mutate(
    incident_t2d_time = t2d_dxdt,  # First t2d diagnosis after baseline
    event = ifelse(!is.na(incident_t2d_time) & incident_t2d_time != Inf, 1, 0),  # Event indicator: 1 if incident t2d occurred
    censor_time = pmin(as.numeric(as.Date(f.40000.0.0)), as.numeric(as.Date(f.191.0.0)), as.numeric(as.Date("2022-11-30")), na.rm = TRUE),  # Censor at death or current date
    time_to_event = ifelse(event == 1, incident_t2d_time - dt_visit0, censor_time - dt_visit0)  # Time to event or censoring
  )

# Check the event distribution
print("Baseline and incident t2d have been defined successfully.")
print(table(incident_data$event))  # Check event distribution



colnames(incident_data)



# Create binary CKM comparison variables
incident_data <- incident_data %>%
  mutate(
    ckm_1_vs_0 = ifelse(ckm == 1, 1, ifelse(ckm == 0, 0, NA)),
    ckm_2_vs_1 = ifelse(ckm == 2, 1, ifelse(ckm == 1, 0, NA)),
    ckm_3_vs_2 = ifelse(ckm == 3, 1, ifelse(ckm == 2, 0, NA)),
    ckm_4_vs_3 = ifelse(ckm == 4, 1, ifelse(ckm == 3, 0, NA))
  )

ckm_comparisons <- c("ckm_1_vs_0", "ckm_2_vs_1", "ckm_3_vs_2", "ckm_4_vs_3")

incident_data$sex.y <- as.factor(incident_data$sex.y)
incident_data$nonwhite <- as.factor(incident_data$nonwhite)
incident_data$fast_hrs <- as.factor(incident_data$fast_hrs)

incident_data$alcohol_cat <- as.factor(incident_data$alcohol_cat)
incident_data$smoke <- as.factor(incident_data$smoke)

# Storage for results
ckm_t2d_results <- data.frame()

for (ckm_comp in ckm_comparisons) {
  data_subset <- incident_data %>%
    filter(!is.na(.data[[ckm_comp]]))
  
  for (i in 1:3) {
    if (i == 1) {
      model <- coxph(Surv(time_to_event, event) ~ get(ckm_comp) + age.y + sex.y, data = data_subset)
    } else if (i == 2) {
      model <- coxph(Surv(time_to_event, event) ~ get(ckm_comp) + age.y + sex.y + nonwhite +fast_hrs+ townsend, data = data_subset)
    } else if (i == 3) {
      model <- coxph(Surv(time_to_event, event) ~ get(ckm_comp) + age.y + sex.y + nonwhite +fast_hrs+ townsend  + PA  + alcohol_cat, data = data_subset)
    }
    
    coef_val <- coef(model)[1]
    se_val <- sqrt(diag(vcov(model)))[1]
    hr <- exp(coef_val)
    p_val <- summary(model)$coef[1, "Pr(>|z|)"]
    lower_ci <- exp(coef_val - 1.96 * se_val)
    upper_ci <- exp(coef_val + 1.96 * se_val)
    
    ckm_t2d_results <- rbind(ckm_t2d_results, data.frame(
      CKM_Comparison = ckm_comp,
      Model = paste0("Model ", i),
      Coefficient = coef_val,
      Hazard_Ratio = hr,
      CI_95_Lower = lower_ci,
      CI_95_Upper = upper_ci,
      SE = se_val,
      P_Value = p_val
    ))
  }
}

# Multiple testing corrections
ckm_t2d_results <- ckm_t2d_results %>%
  group_by(Model) %>%
  mutate(
    Bonferroni_P = p.adjust(P_Value, method = "bonferroni"),
    FDR = p.adjust(P_Value, method = "fdr")
  ) %>%
  ungroup()

# Save results
write.csv(ckm_t2d_results, "/n/home_fasse/txia/UKB_CKM/Data/ckmt2dassociation.csv", row.names = FALSE)
print("CKM vs. t2d associations saved.")




#######split to drug and no drug

# Keep original
incident_data_original <- incident_data

# Create drug/no-drug subsets
incident_data_nodrug <- incident_data_original %>%
  filter(bptreat.y == 0 & statin.y == 0)

incident_data_drug <- incident_data_original %>%
  filter(bptreat.y == 1 | statin.y == 1)

# Define function
run_ckm_vs_t2d <- function(data, group_label) {
  
  # Create binary CKM comparison variables
  data <- data %>%
    mutate(
      ckm_1_vs_0 = ifelse(ckm == 1, 1, ifelse(ckm == 0, 0, NA)),
      ckm_2_vs_1 = ifelse(ckm == 2, 1, ifelse(ckm == 1, 0, NA)),
      ckm_3_vs_2 = ifelse(ckm == 3, 1, ifelse(ckm == 2, 0, NA)),
      ckm_4_vs_3 = ifelse(ckm == 4, 1, ifelse(ckm == 3, 0, NA))
    )
  
  ckm_comparisons <- c("ckm_1_vs_0", "ckm_2_vs_1", "ckm_3_vs_2", "ckm_4_vs_3")
  
  # Storage
  ckm_t2d_results <- data.frame()
  
  for (ckm_comp in ckm_comparisons) {
    data_subset <- data %>%
      filter(!is.na(.data[[ckm_comp]]))
    
    for (i in 1:3) {
      if (i == 1) {
        model <- coxph(Surv(time_to_event, event) ~ get(ckm_comp) + age.y + sex.y, data = data_subset)
      } else if (i == 2) {
        model <- coxph(Surv(time_to_event, event) ~ get(ckm_comp) + age.y + sex.y + nonwhite +fast_hrs+ townsend, data = data_subset)
      } else if (i == 3) {
        model <- coxph(Surv(time_to_event, event) ~ get(ckm_comp) + age.y + sex.y + nonwhite +fast_hrs+ townsend  + PA  + alcohol_cat, data = data_subset)
      }
      
      coef_val <- coef(model)[1]
      se_val <- sqrt(diag(vcov(model)))[1]
      hr <- exp(coef_val)
      p_val <- summary(model)$coef[1, "Pr(>|z|)"]
      lower_ci <- exp(coef_val - 1.96 * se_val)
      upper_ci <- exp(coef_val + 1.96 * se_val)
      
      ckm_t2d_results <- rbind(ckm_t2d_results, data.frame(
        CKM_Comparison = ckm_comp,
        Model = paste0("Model ", i),
        Coefficient = coef_val,
        Hazard_Ratio = hr,
        CI_95_Lower = lower_ci,
        CI_95_Upper = upper_ci,
        SE = se_val,
        P_Value = p_val
      ))
    }
  }
  
  # Multiple testing corrections
  ckm_t2d_results <- ckm_t2d_results %>%
    group_by(Model) %>%
    mutate(
      Bonferroni_P = p.adjust(P_Value, method = "bonferroni"),
      FDR = p.adjust(P_Value, method = "fdr")
    ) %>%
    ungroup()
  
  # Save
  write.csv(ckm_t2d_results,
            paste0("/n/home_fasse/txia/UKB_CKM/Data/ckmt2dassociation_", group_label, ".csv"),
            row.names = FALSE)
  
  print(paste0("Finished CKM vs. t2d associations for ", group_label))
}

# Run separately
run_ckm_vs_t2d(incident_data_nodrug, "nodrug")
run_ckm_vs_t2d(incident_data_drug, "drug")


