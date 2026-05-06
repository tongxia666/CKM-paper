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


# Load the dataset
data1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/combined_proteomics_pheno.csv")
data2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/covariate.csv")

# Merge datasets
ukb <- data1 %>% left_join(data2, by = "f.eid")

ukb$sbp_use <- ifelse(!is.na(ukb$f.4080.0.0.x), ukb$f.4080.0.0.x, ukb$f.93.0.0.x)
summary(ukb$sbp_use)

table(ukb$sex.x)
table(ukb$sex.y)

# Remove rows with missing 'ckm'
ukb <- ukb[!is.na(ukb$ckm), ]

table(ukb$sex.x)
table(ukb$sex.y)

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

ukb$ckm <- as.factor(ukb$ckm)  # Convert ckm to a factor





# Initialize list to store scores
score_table <- data.table(f.eid = ukb$f.eid)

# === Step 2: Define comparisons ===
comparisons <- c("1vs0", "2vs1", "3vs2", "4vs3")

for (comp in comparisons) {
  # File paths
  path_full <- paste0("/n/home_fasse/txia/UKB_CKM/Data/", comp, "_selected_proteins_coeffandintercept_with_info_bootstrap.csv")
  path_0.2  <- paste0("/n/home_fasse/txia/UKB_CKM/Data/", comp, "_selected_proteins_coeffandintercept_with_info_bootstrap_coef0.2_from_refit.csv")

  # Load coefficient files
  coef_full <- fread(path_full)
  coef_0.2  <- fread(path_0.2)

  # === Full model ===
  intercept_full <- coef_full[Protein == "(Intercept)", Coefficient]
  betas_full <- coef_full[Protein != "(Intercept)"]
  proteins_full <- betas_full$Protein
  coef_vector_full <- betas_full$Coefficient
  names(coef_vector_full) <- proteins_full

  # Check missing
  missing_full <- setdiff(proteins_full, colnames(ukb))
  if (length(missing_full) > 0) stop(paste("Missing proteins for", comp, "full:", paste(missing_full, collapse = ", ")))

  # Score calculation
  X_full <- as.matrix(ukb[, ..proteins_full])
  score_full <- 1 / (1 + exp(-(intercept_full + X_full %*% coef_vector_full)))
  score_table[[paste0("manual_protein_score_", comp)]] <- score_full

  # === Threshold 0.2 model ===
  intercept_0.2 <- coef_0.2[Protein == "(Intercept)", Coefficient]
  betas_0.2 <- coef_0.2[Protein != "(Intercept)"]
  proteins_0.2 <- betas_0.2$Protein
  coef_vector_0.2 <- betas_0.2$Coefficient
  names(coef_vector_0.2) <- proteins_0.2

  missing_0.2 <- setdiff(proteins_0.2, colnames(ukb))
  if (length(missing_0.2) > 0) stop(paste("Missing proteins for", comp, "0.2 model:", paste(missing_0.2, collapse = ", ")))

  X_0.2 <- as.matrix(ukb[, ..proteins_0.2])
  score_0.2 <- 1 / (1 + exp(-(intercept_0.2 + X_0.2 %*% coef_vector_0.2)))
  score_table[[paste0("manual_protein_score_", comp, "_0.2")]] <- score_0.2
}

# === Step 3: Save scores table ===
fwrite(score_table, "/n/home_fasse/txia/UKB_CKM/Data/manual_protein_scores_all_ckm.csv")

# === Step 4: Merge scores back to UKB ===
ukb_combined <- merge(ukb, score_table, by = "f.eid", all.x = TRUE)
fwrite(ukb_combined, "/n/home_fasse/txia/UKB_CKM/Data/UKB_combined_with_all_protein_scores.csv")

cat("✅ All CKM scores calculated, saved, and merged into full UKB dataset.\n")


# Load the ckd dataset
ukb_data <- fread("/n/home_fasse/txia/UKB_CKM/Data/ckdvariable.csv")



# List of your protein score columns
score_vars <- c(
  "manual_protein_score_1vs0",
  "manual_protein_score_1vs0_0.2",
  "manual_protein_score_2vs1",
  "manual_protein_score_2vs1_0.2",
  "manual_protein_score_3vs2",
  "manual_protein_score_3vs2_0.2",
  "manual_protein_score_4vs3",
  "manual_protein_score_4vs3_0.2"
)

# Summary + SD for each
for (var in score_vars) {
  cat("-----", var, "-----\n")
  print(summary(ukb_combined[[var]]))
  cat("SD:", round(sd(ukb_combined[[var]], na.rm = TRUE), 4), "\n\n")
}



# List of score variables
score_vars <- c(
  "manual_protein_score_1vs0",
  "manual_protein_score_1vs0_0.2",
  "manual_protein_score_2vs1",
  "manual_protein_score_2vs1_0.2",
  "manual_protein_score_3vs2",
  "manual_protein_score_3vs2_0.2",
  "manual_protein_score_4vs3",
  "manual_protein_score_4vs3_0.2"
)

# Apply INT transformation directly to the same variable names
for (var in score_vars) {
  ukb_combined[, (var) := qnorm(rank(get(var), na.last = "keep") / (.N + 1))]
}




# List of your protein score columns
score_vars <- c(
  "manual_protein_score_1vs0",
  "manual_protein_score_1vs0_0.2",
  "manual_protein_score_2vs1",
  "manual_protein_score_2vs1_0.2",
  "manual_protein_score_3vs2",
  "manual_protein_score_3vs2_0.2",
  "manual_protein_score_4vs3",
  "manual_protein_score_4vs3_0.2"
)

# Summary + SD for each
for (var in score_vars) {
  cat("-----", var, "-----\n")
  print(summary(ukb_combined[[var]]))
  cat("SD:", round(sd(ukb_combined[[var]], na.rm = TRUE), 4), "\n\n")
}



ckd_association<-ukb_combined



# Join the datasets by participant ID (assuming 'f.eid' is the ID column in both)
ckd_association <- ckd_association %>%
  left_join(ukb_data, by = "f.eid")

# Define ckd-related outcomes (exclude CKD-related outcomes)
ckd_outcomes <- c(
  "f.132032.0.0", "f.42026.0.0"
  
)

# Convert date fields to numeric (dates are assumed to be in YYYY-MM-DD format)
for (outcome in ckd_outcomes) {
  ckd_association[[outcome]] <- as.numeric(as.Date(ckd_association[[outcome]]))
}

# Define baseline visit date (assumed to be stored in `f.21842`)
ckd_association$dt_visit0 <- as.numeric(as.Date(ckd_association$f.21842.0.0))

# Step 1: Identify first ckd diagnosis date and handle Inf values
ckd_association$ckd_dxdt <- apply(
  as.matrix(ckd_association[, ..ckd_outcomes]), 
  1, 
  function(x) {
    min_val <- min(x, na.rm = TRUE)  # Find the earliest non-NA diagnosis date
    if (is.infinite(min_val)) NA else min_val  # Replace Inf with NA
  }
)

# Replace NA with Inf where required explicitly for handling missing cases
ckd_association$ckd_dxdt[is.na(ckd_association$ckd_dxdt)] <- Inf

# Step 2: Determine baseline ckd (prevalent cases: diagnosis before baseline visit)
ckd_association$ckd_bsln <- ifelse(
  ckd_association$ckd_dxdt < ckd_association$dt_visit0,
  1, 0
)

# Step 3: Exclude participants with baseline ckd
incident_data <- ckd_association %>%
  filter(ckd_bsln == 0)

#############delete ckm 4
        #incident_data <-incident_data[ckm != 4]

# Step 4: Define incident ckd and censoring
incident_data <- incident_data %>%
  mutate(
    incident_ckd_time = ckd_dxdt,  # First ckd diagnosis after baseline
    event = ifelse(!is.na(incident_ckd_time) & incident_ckd_time != Inf, 1, 0),  # Event indicator: 1 if incident ckd occurred
    censor_time = pmin(as.numeric(as.Date(f.40000.0.0)), as.numeric(as.Date(f.191.0.0)), as.numeric(as.Date("2022-11-30")), na.rm = TRUE),  # Censor at death or current date
    time_to_event = ifelse(event == 1, incident_ckd_time - dt_visit0, censor_time - dt_visit0)  # Time to event or censoring
  )

# Check the event distribution
print("Baseline and incident ckd have been defined successfully.")
print(table(incident_data$event))  # Check event distribution

# Save the processed dataset for downstream analysis
#processed_output <- "/n/home_fasse/txia/UKB_CKM/Data/processed_ckd_data.csv"
#write.table(incident_data, file = processed_output, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")

#print("Processed ckd data saved successfully.")






library(data.table)
library(dplyr)
library(survival)
library(splines)

# List of protein score variables
protein_scores <- c(
  "manual_protein_score_1vs0",
  "manual_protein_score_1vs0_0.2",
  "manual_protein_score_2vs1",
  "manual_protein_score_2vs1_0.2",
  "manual_protein_score_3vs2",
  "manual_protein_score_3vs2_0.2",
  "manual_protein_score_4vs3",
  "manual_protein_score_4vs3_0.2"
)

# Initialize result tables
all_proteomics_score_results <- data.frame()
all_tertile_results <- data.frame()

run_models_on_subset <- function(data, label) {
  data$sex.y <- as.factor(data$sex.y)
  data$nonwhite <- as.factor(data$nonwhite)
  data$fast_hrs <- as.factor(data$fast_hrs)
  
  data$alcohol_cat <- as.factor(data$alcohol_cat)
  data$smoke <- as.factor(data$smoke)

  for (score_var in protein_scores) {
    cat("Processing:", score_var, "in group", label, "\n")

    # --- 1. Continuous Cox Models ---
    for (i in 1:3) {
      formula_i <- switch(i,
        as.formula(paste("Surv(time_to_event, event) ~", score_var, "+ age.y + sex.y")),
        as.formula(paste("Surv(time_to_event, event) ~", score_var, "+ age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat")),
        as.formula(paste("Surv(time_to_event, event) ~", score_var, "+ age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR"))
      )
      model <- coxph(formula_i, data = data)

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

    # --- 2. Tertile Models ---
    tertile_var <- paste0(score_var, "_tertile")
    data[[tertile_var]] <- cut(
      data[[score_var]],
      breaks = quantile(data[[score_var]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("T1", "T2", "T3"),
      include.lowest = TRUE
    )
    data[[tertile_var]] <- relevel(as.factor(data[[tertile_var]]), ref = "T1")

    for (i in 1:3) {
      formula_i <- switch(i,
        as.formula(paste("Surv(time_to_event, event) ~", tertile_var, "+ age.y + sex.y")),
        as.formula(paste("Surv(time_to_event, event) ~", tertile_var, "+ age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat")),
        as.formula(paste("Surv(time_to_event, event) ~", tertile_var, "+ age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR"))
      )
      model <- coxph(formula_i, data = data)
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

    # --- 3. Spline Model (Model 3 only) ---
    spline_model <- coxph(Surv(time_to_event, event) ~ ns(get(score_var), df = 4) +
      age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR, data = data)

    min_val <- min(data[[score_var]], na.rm = TRUE)
    max_val <- max(data[[score_var]], na.rm = TRUE)
    new_data <- data.frame(score = seq(min_val, max_val, length.out = 100))
    colnames(new_data)[1] <- score_var

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

    spline_pred <- predict(spline_model, newdata = new_data, type = "lp", se.fit = TRUE)
    new_data$HR <- exp(spline_pred$fit)
    new_data$HR_lower <- exp(spline_pred$fit - 1.96 * spline_pred$se.fit)
    new_data$HR_upper <- exp(spline_pred$fit + 1.96 * spline_pred$se.fit)

    fwrite(new_data, paste0("/n/home_fasse/txia/UKB_CKM/Data/spline_ckd_alldata_", score_var, "_model3_", label, ".csv"))
  }
}

# === Apply models ===
run_models_on_subset(incident_data, "overall")
run_models_on_subset(subset(incident_data, bptreat.y == 0 & statin.y == 0), "nodrug")
run_models_on_subset(subset(incident_data, bptreat.y == 1 | statin.y == 1), "drug")

# === Adjust and export results ===
all_tertile_results$Bonferroni_P <- p.adjust(all_tertile_results$P_Value, method = "bonferroni")
all_tertile_results$FDR <- p.adjust(all_tertile_results$P_Value, method = "fdr")

fwrite(all_proteomics_score_results, "/n/home_fasse/txia/UKB_CKM/Data/ckd_alldata_all_scores_all_models.csv")
fwrite(all_tertile_results, "/n/home_fasse/txia/UKB_CKM/Data/ckd_alldata_all_scores_tertile_models.csv")

cat("✅ All models completed for all protein scores (continuous, tertile, spline).\n")



###########proteins and ckd

# ------------------------------------------------------------------
# Required inputs present in your session:
# - proteomics_info : annotation table with a 'UniProt' column (plus names you want to carry along)
# - incident_data   : analysis dataset with protein columns + covariates:
#                     time_to_event, event, age.y, sex.y, nonwhite, townsend,
#                     smoke, PA, BMI, alcohol_cat
# ------------------------------------------------------------------

# 1) Files for all contrasts
contrast_files <- c(
  "1vs0" = "/n/home_fasse/txia/UKB_CKM/Data/1vs0_selected_proteins_with_info_bootstrap.csv",
  "2vs1" = "/n/home_fasse/txia/UKB_CKM/Data/2vs1_selected_proteins_with_info_bootstrap.csv",
  "3vs2" = "/n/home_fasse/txia/UKB_CKM/Data/3vs2_selected_proteins_with_info_bootstrap.csv",
  "4vs3" = "/n/home_fasse/txia/UKB_CKM/Data/4vs3_selected_proteins_with_info_bootstrap.csv"
)

# 2) Read and merge annotation for each contrast
merged_list <- lapply(names(contrast_files), function(ct){
  imp <- fread(contrast_files[[ct]])  # must contain at least: Protein, Coefficient
  # Merge on UniProt to bring names/labels; keep all selected proteins
  out <- merge(imp, proteomics_info, by.x = "Protein", by.y = "UniProt", all.x = TRUE)
  out$Contrast <- ct
  out
})
names(merged_list) <- names(contrast_files)

# 3) Thresholds for selecting proteins by absolute coefficient
thresholds <- c(0, 0.2)

# 4) Helper: fit 3 Cox models for ONE protein with PER-SD scaling
fit_three_models_perSD <- function(protein, dat){
  x  <- dat[[protein]]
  mu <- mean(x, na.rm = TRUE)
  sg <- sd(x,   na.rm = TRUE)
  if (!is.finite(sg) || sg == 0) return(NULL)

  # Standardize protein -> per 1 SD
  dat$._prot_z <- as.numeric((x - mu) / sg)
  term <- "._prot_z"

  f1 <- as.formula(sprintf("Surv(time_to_event, event) ~ %s + age.y + sex.y", term))
  f2 <- as.formula(sprintf("Surv(time_to_event, event) ~ %s + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat", term))
  f3 <- as.formula(sprintf("Surv(time_to_event, event) ~ %s + age.y + sex.y + nonwhite +fast_hrs+ townsend + smoke + PA + BMI + alcohol_cat +sbp_use +  f.30690.0.0.x + f.30760.0.0.x +  f.30740.0.0.x + eGFR", term))

  fits <- list(
    "Model 1" = try(coxph(f1, data = dat), silent = TRUE),
    "Model 2" = try(coxph(f2, data = dat), silent = TRUE),
    "Model 3" = try(coxph(f3, data = dat), silent = TRUE)
  )

  pull <- function(fit, mn){
    if (inherits(fit, "try-error")) return(NULL)
    s <- summary(fit)
    est <- s$coef[term, "coef"]; se <- s$coef[term, "se(coef)"]; p <- s$coef[term, "Pr(>|z|)"]
    hr <- exp(est); lo <- exp(est - 1.96*se); hi <- exp(est + 1.96*se)
    data.frame(
      Protein = protein,
      Model = mn,
      Coefficient = est,
      SE = se,
      Hazard_Ratio = hr,
      CI_95_Lower = lo,
      CI_95_Upper = hi,
      P_Value = p,
      Protein_Mean = mu,
      Protein_SD = sg,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, Map(pull, fits, names(fits)))
}

# 5) Main loop over thresholds and contrasts
for (threshold in thresholds) {
  all_contrasts_results <- list()

  for (ct in names(merged_list)) {
    merged_proteins <- merged_list[[ct]]  # has Protein, Coefficient, annotation cols, Contrast

    # Select proteins meeting the threshold |Coefficient| >= threshold
    proteins_threshold <- merged_proteins %>%
      filter(abs(Coefficient) >= threshold) %>%
      pull(Protein) %>% unique()

    # Keep only proteins present in incident_data
    valid_proteins <- intersect(proteins_threshold, colnames(incident_data))
    if (length(valid_proteins) == 0) {
      cat("No valid proteins for", ct, "at threshold", threshold, "\n")
      next
    }

    # Fit per-SD models protein-by-protein
    protein_results <- do.call(rbind, lapply(valid_proteins, fit_three_models_perSD, dat = incident_data))
    if (is.null(protein_results) || nrow(protein_results) == 0) {
      cat("No model results for", ct, "at threshold", threshold, "\n")
      next
    }

    # Merge annotation and contrast label
    annotated_results <- merge(
      protein_results,
      merged_proteins,
      by = "Protein",
      all.x = TRUE
    )
    annotated_results$Contrast <- ct  # ensure present

    # P-value adjustments WITHIN each Contrast × Model
    annotated_results <- annotated_results %>%
      group_by(Contrast, Model) %>%
      mutate(
        Bonferroni_P = p.adjust(P_Value, method = "bonferroni"),
        FDR          = p.adjust(P_Value, method = "fdr")
      ) %>%
      ungroup()

    # Save per-contrast full table
    write.csv(
      annotated_results,
      sprintf("/n/home_fasse/txia/UKB_CKM/Data/wholedata_ckd_ckm%s_association_threshold_%s.csv", ct, threshold),
      row.names = FALSE
    )

    # Model 3 significant exports (per contrast)
    model3_res <- annotated_results %>% filter(Model == "Model 3")
    write.csv(
      model3_res %>% filter(Bonferroni_P < 0.05),
      sprintf("/n/home_fasse/txia/UKB_CKM/Data/wholedata_ckd_Significant_ckm%s_Model3_Bonferroni_Threshold_%s.csv", ct, threshold),
      row.names = FALSE
    )
    write.csv(
      model3_res %>% filter(FDR < 0.05),
      sprintf("/n/home_fasse/txia/UKB_CKM/Data/wholedata_ckd_Significant_ckm%s_Model3_FDR_Threshold_%s.csv", ct, threshold),
      row.names = FALSE
    )

    all_contrasts_results[[ct]] <- annotated_results
    cat("Done:", ct, "at threshold", threshold, "\n")
  }

  # Optional: one combined file across all contrasts for this threshold
  if (length(all_contrasts_results)) {
    combined <- dplyr::bind_rows(all_contrasts_results)
    write.csv(
      combined,
      sprintf("/n/home_fasse/txia/UKB_CKM/Data/wholedata_ckd_ckm_allcontrasts_association_threshold_%s.csv", threshold),
      row.names = FALSE
    )
  } else {
    cat("No contrasts produced results at threshold", threshold, "\n")
  }
}


####event and time distribution

# List of protein score variables
protein_scores <- c(
  "manual_protein_score_1vs0",
  "manual_protein_score_1vs0_0.2",
  "manual_protein_score_2vs1",
  "manual_protein_score_2vs1_0.2",
  "manual_protein_score_3vs2",
  "manual_protein_score_3vs2_0.2",
  "manual_protein_score_4vs3",
  "manual_protein_score_4vs3_0.2"
)

# Initialize containers
all_overall_summaries <- list()
all_tertile_summaries <- list()

for (score_var in protein_scores) {
  cat("\n============================\n")
  cat("Processing:", score_var, "\n")

  # Skip if not present or all NA
  if (!score_var %in% colnames(incident_data) || all(is.na(incident_data[[score_var]]))) {
    cat("Skipping", score_var, "- missing or all NA.\n")
    next
  }

  # Create tertile grouping
  incident_data <- incident_data %>%
    mutate(
      tertile_group = cut(
        .data[[score_var]],
        breaks = quantile(.data[[score_var]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
        labels = c("T1", "T2", "T3"),
        include.lowest = TRUE
      )
    )

  # Overall event/time summary
  overall_summary <- incident_data %>%
    summarise(
      N_total = n(),
      N_events = sum(event == 1, na.rm = TRUE),
      Event_rate = mean(event == 1, na.rm = TRUE),
      Time_to_event_mean = mean(time_to_event, na.rm = TRUE),
      Time_to_event_median = median(time_to_event, na.rm = TRUE)
    ) %>%
    mutate(Protein_Score = score_var)

  # Summary by tertile
  tertile_summary <- incident_data %>%
    group_by(tertile_group) %>%
    summarise(
      N = n(),
      Events = sum(event == 1, na.rm = TRUE),
      Event_rate = mean(event == 1, na.rm = TRUE),
      Time_to_event_mean = mean(time_to_event, na.rm = TRUE),
      Time_to_event_median = median(time_to_event, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Protein_Score = score_var)

  # Save
  all_overall_summaries[[score_var]] <- overall_summary
  all_tertile_summaries[[score_var]] <- tertile_summary

  # Print to console
  cat("---- Overall Summary ----\n")
  print(overall_summary)

  cat("---- Tertile Summary ----\n")
  print(tertile_summary)
}



