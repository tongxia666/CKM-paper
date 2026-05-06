# Load necessary libraries
options(repos = c(CRAN = "https://cloud.r-project.org"))
myPaths <- .libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs", myPaths)  # add new path
.libPaths(myPaths)  # reassign them

library(splines, lib.loc = "/n/home_fasse/txia/R/libs")
library(glmnet, lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
library(pROC, lib.loc = "/n/home_fasse/txia/R/libs")
library(parallel, lib.loc = "/n/home_fasse/txia/R/libs")
library(survival, lib.loc = "/n/home_fasse/txia/R/libs")

# Load the proteomics annotation file
print("Loading protein information...")
proteomics_info <- fread("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv")

# Add ".x" to the UniProt column in proteomics_info to match logistic results
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")

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
t2d_association <- fread("/n/home_fasse/txia/UKB_CKM/Data/2_vs_1_combined_training_optimization_testing.csv")
ukb_data <- fread("/n/home_fasse/txia/UKB_CKM/Data/t2dvariable.csv")


t2d_association$sbp_use <- ifelse(!is.na(t2d_association$f.4080.0.0.x), t2d_association$f.4080.0.0.x, t2d_association$f.93.0.0.x)
summary(t2d_association$sbp_use)

table(t2d_association$sex.x)
table(t2d_association$sex.y)


# Remove rows with missing 'ckm'
t2d_association <- t2d_association[!is.na(t2d_association$ckm), ] 

summary(t2d_association$protein_score)
sd(t2d_association$protein_score, na.rm = TRUE)

# Standardize protein_score to mean = 0, SD = 1
#t2d_association[, protein_score := scale(protein_score)]
t2d_association[, protein_score := qnorm(rank(protein_score, na.last = "keep") / (.N + 1))]
t2d_association[, protein_score_0.2 := qnorm(rank(protein_score_0.2, na.last = "keep") / (.N + 1))]



summary(t2d_association$protein_score_0.2)
sd(t2d_association$protein_score_0.2, na.rm = TRUE)




# 1. Basic summary
summary_stats <- summary(t2d_association$protein_score)

# 2. Compute IQR and outlier bounds (boxplot method)
Q1 <- summary_stats["1st Qu."]
Q3 <- summary_stats["3rd Qu."]
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# 3. Identify outliers
outliers <- t2d_association$protein_score[
  t2d_association$protein_score < lower_bound |
  t2d_association$protein_score > upper_bound
]

# 4. Compute 1st and 99th percentiles
p1 <- quantile(t2d_association$protein_score, 0.01, na.rm = TRUE)
p99 <- quantile(t2d_association$protein_score, 0.99, na.rm = TRUE)

# 5. Display everything
cat("=== Summary ===\n")
print(summary_stats)

cat("\n=== 1st and 99th Percentiles ===\n")
cat("1% percentile:", round(p1, 4), "\n")
cat("99% percentile:", round(p99, 4), "\n")

cat("\n=== Outliers Based on IQR Rule ===\n")
print(outliers)


# Join the datasets by participant ID (assuming 'f.eid' is the ID column in both)
t2d_association <- t2d_association %>%
  left_join(ukb_data, by = "f.eid")

# Define t2d-related outcomes (exclude CKD-related outcomes)
t2d_outcomes <- c(
  "f.130708.0.0"  # diabetes
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

# Save the processed dataset for downstream analysis
#processed_output <- "/n/home_fasse/txia/UKB_CKM/Data/processed_t2d_data.csv"
#write.table(incident_data, file = processed_output, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")

#print("Processed t2d data saved successfully.")

# Load individual proteins with non-zero coefficients
important_proteins <- fread("/n/home_fasse/txia/UKB_CKM/Data/2vs1_selected_proteins_with_info_bootstrap.csv")
# Merge protein annotation with selected proteins
merged_proteins <- merge(
  important_proteins, 
  proteomics_info, 
  by.x = "Protein", 
  by.y = "UniProt", 
  all.x = TRUE
)

incident_data$sex.y <- as.factor(incident_data$sex.y)
incident_data$nonwhite <- as.factor(incident_data$nonwhite)
incident_data$fast_hrs <- as.factor(incident_data$fast_hrs)
incident_data$alcohol_cat <- as.factor(incident_data$alcohol_cat)
incident_data$smoke <- as.factor(incident_data$smoke)
  
# Thresholds for protein selection
thresholds <- c(0, 0.2)

# Prepare storage for results
results <- list()

# Perform analyses for each threshold
for (threshold in thresholds) {
  # Select proteins meeting the threshold
  proteins_threshold <- merged_proteins %>%
    filter(abs(Coefficient) >= threshold) %>%
    pull(Protein)
  
  # Check if proteins exist in the dataset
  valid_proteins <- intersect(proteins_threshold, colnames(incident_data))
  
  if (length(valid_proteins) == 0) {
    cat("No valid proteins found for threshold:", threshold, "\n")
    next
  }
  
  # Association with individual proteins
  protein_results <- data.frame()
  
  for (protein in valid_proteins) {
    # Fit Cox models
    model1 <- coxph(Surv(time_to_event, event) ~ get(protein) + age.y + sex.y, data = incident_data)
    model2 <- coxph(Surv(time_to_event, event) ~ get(protein) + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI+ alcohol_cat, data = incident_data)
    model3 <- coxph(Surv(time_to_event, event) ~ get(protein) + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI+ alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR, data = incident_data)
    
    # Extract results for all models
    models <- list(model1, model2, model3)
    for (i in seq_along(models)) {
      model <- models[[i]]
      coef_val <- coef(model)[1]
      se_val <- sqrt(diag(vcov(model)))[1]
      hr <- exp(coef_val)
      p_val <- summary(model)$coef[1, "Pr(>|z|)"]
      lower_ci <- exp(coef_val - 1.96 * se_val)
      upper_ci <- exp(coef_val + 1.96 * se_val)
      
      # Add results to protein_results
      protein_results <- rbind(protein_results, data.frame(
        Protein = protein,
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
  
  # Apply Bonferroni and FDR adjustments
  protein_results <- protein_results %>%
    group_by(Model) %>%
    mutate(
      Bonferroni_P = p.adjust(P_Value, method = "bonferroni"),
      FDR = p.adjust(P_Value, method = "fdr")
    ) %>%
    ungroup()
  
  # Merge with annotation information
  annotated_results <- merge(
    protein_results, 
    merged_proteins, 
    by.x = "Protein", 
    by.y = "Protein", 
    all.x = TRUE
  )
  
  # Save all results for this threshold
  write.csv(
    annotated_results,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/t2d_ckm2vs1_association_threshold_", threshold, ".csv"),
    row.names = FALSE
  )
  
  # Extract Model 3 significant results
  model3_results <- annotated_results %>% filter(Model == "Model 3")
  
  # Save Model 3 Bonferroni < 0.05 results
  bonferroni_results <- model3_results %>% filter(Bonferroni_P < 0.05)
  write.csv(
    bonferroni_results,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Significant_ckm2vs1_Model3_Bonferroni_Threshold_", threshold, ".csv"),
    row.names = FALSE
  )
  
  # Save Model 3 FDR < 0.05 results
  fdr_results <- model3_results %>% filter(FDR < 0.05)
  write.csv(
    fdr_results,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Significant_ckm2vs1_Model3_FDR_Threshold_", threshold, ".csv"),
    row.names = FALSE
  )
  
  cat("Analysis complete for threshold:", threshold, "\n")
}


# Step: Add Continuous, Tertile, and Spline Models for protein_score and Threshold Scores
print("Running all Cox models for continuous, tertile, and spline protein scores...")

# List of protein score variants
protein_scores <- c("protein_score",  "protein_score_0.2")

# Initialize combined result tables
all_proteomics_score_results <- data.frame()
all_tertile_results <- data.frame()

library(splines)

run_models_on_subset <- function(data, label) {
  # Ensure factor variables are properly set
  data$sex.y <- as.factor(data$sex.y)
  data$nonwhite <- as.factor(data$nonwhite)
  data$fast_hrs <- as.factor(data$fast_hrs)
  data$alcohol_cat <- as.factor(data$alcohol_cat)
  data$smoke <- as.factor(data$smoke)

  for (score_var in protein_scores) {
    print(paste("Processing:", score_var, "in group", label))

    ### 1. Continuous Cox Models
    for (i in 1:3) {
      if (i == 1) {
        model <- coxph(Surv(time_to_event, event) ~ get(score_var) + age.y + sex.y, data = data)
      } else if (i == 2) {
        model <- coxph(Surv(time_to_event, event) ~ get(score_var) + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat, data = data)
      } else if (i == 3) {
        model <- coxph(Surv(time_to_event, event) ~ get(score_var) + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR, data = data)
      }

      coef_val <- coef(model)[1]
      se_val <- sqrt(diag(vcov(model)))[1]
      hr <- exp(coef_val)
      p_val <- summary(model)$coef[1, "Pr(>|z|)"]
      lower_ci <- exp(coef_val - 1.96 * se_val)
      upper_ci <- exp(coef_val + 1.96 * se_val)

      all_proteomics_score_results <<- rbind(all_proteomics_score_results, data.frame(
        Group = label,
        Score = score_var,
        Model = paste0("Model ", i),
        Coefficient = coef_val,
        Hazard_Ratio = hr,
        CI_95_Lower = lower_ci,
        CI_95_Upper = upper_ci,
        SE = se_val,
        P_Value = p_val,
        Bonferroni_P = p.adjust(p_val, method = "bonferroni"),
        FDR = p.adjust(p_val, method = "fdr")
      ))
    }

    ### 2. Tertile Models
    tertile_var <- paste0(score_var, "_tertile")
data <- data %>%
  mutate(
    !!tertile_var := factor(
      dplyr::ntile(.data[[score_var]], 5),
      labels = c("T1","T2","T3","T4","T5"),
      ordered = FALSE
    ),
    !!tertile_var := stats::relevel(.data[[tertile_var]], ref = "T1")
  )

    for (i in 1:3) {
      if (i == 1) {
        model <- coxph(Surv(time_to_event, event) ~ get(tertile_var) + age.y + sex.y, data = data)
      } else if (i == 2) {
        model <- coxph(Surv(time_to_event, event) ~ get(tertile_var) + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat, data = data)
      } else if (i == 3) {
        model <- coxph(Surv(time_to_event, event) ~ get(tertile_var) + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR, data = data)
      }

      model_summary <- summary(model)$coef
      model_confint <- confint(model)

      for (j in 1:nrow(model_summary)) {
        all_tertile_results <<- rbind(all_tertile_results, data.frame(
          Group = label,
          Score = score_var,
          Model = paste0("Model ", i),
          Comparison = rownames(model_summary)[j],
          HR = exp(model_summary[j, "coef"]),
          CI_95_Lower = exp(model_confint[j, 1]),
          CI_95_Upper = exp(model_confint[j, 2]),
          SE = model_summary[j, "se(coef)"],
          P_Value = model_summary[j, "Pr(>|z|)"]
        ))
      }
    }

    ### 3. Spline Model
    ref_sex <- levels(data$sex.y)[1]
    ref_nonwhite <- levels(data$nonwhite)[1]
    ref_fast_hrs <- levels(data$fast_hrs)[1]
    ref_alcohol <- levels(data$alcohol_cat)[1]
    ref_smoke<- levels(data$smoke)[1]

# Fit spline model
spline_model <- coxph(Surv(time_to_event, event) ~ ns(get(score_var), df = 4) + 
  age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR, 
  data = data)

# Use full score range from full dataset
min_val <- min(data[[score_var]], na.rm = TRUE)
max_val <- max(data[[score_var]], na.rm = TRUE)

new_data <- data.frame(score = seq(min_val, max_val, length.out = 100))
colnames(new_data)[1] <- score_var

# Set covariates
new_data$age.y <- median(data$age.y, na.rm = TRUE)
new_data$sex.y <- factor(rep(levels(data$sex.y)[1], 100), levels = levels(data$sex.y))
new_data$nonwhite <- factor(rep(levels(data$nonwhite)[1], 100), levels = levels(data$nonwhite))
new_data$fast_hrs <- factor(rep(levels(data$fast_hrs)[1], 100), levels = levels(data$fast_hrs))
new_data$townsend <- median(data$townsend, na.rm = TRUE)
new_data$smoke <- factor(rep(levels(data$smoke)[1], 100), levels = levels(data$smoke))
new_data$PA <- median(data$PA, na.rm = TRUE)
new_data$BMI <- median(data$BMI, na.rm = TRUE)
new_data$alcohol_cat <- factor(rep(levels(data$alcohol_cat)[1], 100), levels = levels(data$alcohol_cat))

new_data$sbp_use <- median(as.numeric(data$sbp_use), na.rm = TRUE)

new_data$`f.30690.0.0.x` <- median(as.numeric(data$`f.30690.0.0.x`), na.rm = TRUE)

new_data$`f.30760.0.0.x` <- median(as.numeric(data$`f.30760.0.0.x`), na.rm = TRUE)

new_data$`f.30740.0.0.x` <- median(as.numeric(data$`f.30740.0.0.x`), na.rm = TRUE)

new_data$eGFR <- median(as.numeric(data$eGFR), na.rm = TRUE)

# Predict log(HR) and calculate CI correctly
spline_pred <- predict(spline_model, newdata = new_data, type = "lp", se.fit = TRUE)

new_data$HR <- exp(spline_pred$fit)
new_data$HR_lower <- exp(spline_pred$fit - 1.96 * spline_pred$se.fit)
new_data$HR_upper <- exp(spline_pred$fit + 1.96 * spline_pred$se.fit)



    fwrite(new_data, paste0("/n/home_fasse/txia/UKB_CKM/Data/spline_t2d_ckm2vs1_", score_var, "_model3_", label, ".csv"))
  }
}

# Apply the function to all data and stratified by drug use
run_models_on_subset(incident_data, "overall")
run_models_on_subset(subset(incident_data, bptreat.y == 0 & statin.y == 0), "nodrug")
run_models_on_subset(subset(incident_data, bptreat.y == 1 | statin.y == 1), "drug")

# Adjust p-values for the full result
all_tertile_results$Bonferroni_P <- p.adjust(all_tertile_results$P_Value, method = "bonferroni")
all_tertile_results$FDR <- p.adjust(all_tertile_results$P_Value, method = "fdr")

# Export results
write.csv(
  all_proteomics_score_results,
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_ckm2vs1_proteomics_score_all_models_all_thresholds_stratified.csv",
  row.names = FALSE
)

write.csv(
  all_tertile_results,
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_ckm2vs1_proteomics_score_tertile_models_all_thresholds_stratified.csv",
  row.names = FALSE
)

print("Completed all models for all protein score thresholds including spline models, stratified by drug use.")


#########EVENT AND TIME DISTIRBUTION
library(dplyr)
library(survival)

# Choose the protein score of interest (e.g., "protein_score")
score_var <- "protein_score"

# Create tertiles based on the full dataset
incident_data <- incident_data %>%
  mutate(
    tertile = cut(
      .data[[score_var]],
      breaks = quantile(.data[[score_var]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("T1", "T2", "T3"),
      include.lowest = TRUE
    )
  )

# 1. Overall event distribution summary
overall_summary <- incident_data %>%
  summarise(
    N_total = n(),
    N_events = sum(event == 1, na.rm = TRUE),
    Event_rate = mean(event == 1, na.rm = TRUE),
    Time_to_event_mean = mean(time_to_event, na.rm = TRUE),
    Time_to_event_median = median(time_to_event, na.rm = TRUE)
  )

# 2. By tertile of protein score
tertile_summary <- incident_data %>%
  group_by(tertile) %>%
  summarise(
    N = n(),
    Events = sum(event == 1, na.rm = TRUE),
    Event_rate = mean(event == 1, na.rm = TRUE),
    Time_to_event_mean = mean(time_to_event, na.rm = TRUE),
    Time_to_event_median = median(time_to_event, na.rm = TRUE)
  )

# Print results
cat("==== Overall Summary ====\n")
print(overall_summary)

cat("\n==== Summary by Protein Score Tertile ====\n")
print(tertile_summary)



#############baseline charactieristics by score
incident_data <- incident_data %>%
  mutate(
    protein_score_group = case_when(
      protein_score < -2 ~ "< -2",
      protein_score >= -2 & protein_score <= 2 ~ "-2 to 2",
      protein_score > 2 ~ "> 2"
    )
  )


summary_table <- incident_data %>%
  group_by(protein_score_group) %>%
  summarize(
    n = n(),
    mean_age = mean(age.y, na.rm = TRUE),
    sd_age = sd(age.y, na.rm = TRUE),
    mean_townsend = mean(townsend, na.rm = TRUE),
    sd_townsend = sd(townsend, na.rm = TRUE),
    mean_PA = mean(PA, na.rm = TRUE),
    sd_PA = sd(PA, na.rm = TRUE),
    mean_BMI = mean(BMI, na.rm = TRUE),
    sd_BMI = sd(BMI, na.rm = TRUE),
    
        mean_sbp_use        = mean(sbp_use, na.rm = TRUE),
sd_sbp_use          = sd(sbp_use, na.rm = TRUE),

mean_f30690_0_0_x   = mean(`f.30690.0.0.x`, na.rm = TRUE),
sd_f30690_0_0_x     = sd(`f.30690.0.0.x`, na.rm = TRUE),

mean_f30760_0_0_x   = mean(`f.30760.0.0.x`, na.rm = TRUE),
sd_f30760_0_0_x     = sd(`f.30760.0.0.x`, na.rm = TRUE),

mean_f30740_0_0_x   = mean(`f.30740.0.0.x`, na.rm = TRUE),
sd_f30740_0_0_x     = sd(`f.30740.0.0.x`, na.rm = TRUE),

mean_eGFR           = mean(eGFR, na.rm = TRUE),
sd_eGFR             = sd(eGFR, na.rm = TRUE)
  )

# Print the full table with all columns shown
print(summary_table, width = Inf)


# SEX
incident_data %>%
  group_by(protein_score_group, sex.y) %>%
  summarize(n = n()) %>%
  group_by(protein_score_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

# NONWHITE
incident_data %>%
  group_by(protein_score_group, nonwhite) %>%
  summarize(n = n()) %>%
  group_by(protein_score_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))
  
  
  
  # fast_hrs
incident_data %>%
  group_by(protein_score_group, fast_hrs) %>%
  summarize(n = n()) %>%
  group_by(protein_score_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))
  

# SMOKING
incident_data %>%
  group_by(protein_score_group, smoke) %>%
  summarize(n = n()) %>%
  group_by(protein_score_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

# ALCOHOL
incident_data %>%
  group_by(protein_score_group, alcohol_cat) %>%
  summarize(n = n()) %>%
  group_by(protein_score_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

# t2d
incident_data %>%
  group_by(protein_score_group, event) %>%
  summarize(n = n()) %>%
  group_by(protein_score_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))



########summary of protein score
# 1. Basic summary
summary_stats <- summary(incident_data$protein_score)

# 2. Compute IQR and outlier bounds (boxplot method)
Q1 <- summary_stats["1st Qu."]
Q3 <- summary_stats["3rd Qu."]
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# 3. Identify outliers
outliers <- incident_data$protein_score[
  incident_data$protein_score < lower_bound |
  incident_data$protein_score > upper_bound
]

# 4. Compute 1st and 99th percentiles
p1 <- quantile(incident_data$protein_score, 0.01, na.rm = TRUE)
p99 <- quantile(incident_data$protein_score, 0.99, na.rm = TRUE)

# 5. Display everything
cat("=== Summary ===\n")
print(summary_stats)

cat("\n=== 1st and 99th Percentiles ===\n")
cat("1% percentile:", round(p1, 4), "\n")
cat("99% percentile:", round(p99, 4), "\n")

cat("\n=== Outliers Based on IQR Rule ===\n")
print(outliers)  


########summary of protein score_0.2
# 1. Basic summary
summary_stats <- summary(incident_data$protein_score_0.2)

# 2. Compute IQR and outlier bounds (boxplot method)
Q1_0.2 <- summary_stats["1st Qu."]
Q3_0.2 <- summary_stats["3rd Qu."]
IQR_0.2 <- Q3_0.2 - Q1_0.2
lower_bound_0.2 <- Q1_0.2 - 1.5 * IQR_0.2
upper_bound_0.2 <- Q3_0.2 + 1.5 * IQR_0.2

# 3. Identify outliers
outliers_0.2 <- incident_data$protein_score_0.2[
  incident_data$protein_score_0.2 < lower_bound_0.2 |
  incident_data$protein_score_0.2 > upper_bound_0.2
]

# 4. Compute 1st and 99th percentiles
p1_0.2 <- quantile(incident_data$protein_score_0.2, 0.01, na.rm = TRUE)
p99_0.2 <- quantile(incident_data$protein_score_0.2, 0.99, na.rm = TRUE)

# 5. Display everything
cat("=== Summary ===\n")
print(summary_stats)

cat("\n=== 1st and 99th Percentiles ===\n")
cat("1% percentile:", round(p1_0.2, 4), "\n")
cat("99% percentile:", round(p99_0.2, 4), "\n")

cat("\n=== Outliers Based on IQR Rule ===\n")
print(outliers_0.2)




