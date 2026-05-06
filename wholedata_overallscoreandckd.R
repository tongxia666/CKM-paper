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
cat("Missing percentage for each proteomics variable:\n")
print(missing_percentage)

# Identify proteomics variables with more than 80% missing ukb
high_missing_proteins <- names(missing_percentage[missing_percentage > 80])
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

# Remove rows with missing 'ckm' (again, same as your script)
ukb <- ukb[!is.na(ukb$ckm), ]

# Check for NA values in 'ckm' column
if (any(is.na(ukb$ckm))) {
  print("The 'ckm' column has NA values.")
} else {
  print("No NA values found in the 'ckm' column.")
}
ukb$ckm <- as.factor(ukb$ckm)  # Convert ckm to a factor

# ============================================================
# SCORE PART (CHANGED ONLY HERE)
# - Replace stage comparisons with OVERALL scores
# - Score = sum(beta * protein) (NO intercept, NO exp/log)
# ============================================================

# Initialize list to store scores
score_table <- data.table(f.eid = ukb$f.eid)

# Overall score files
overall_score_files <- list(
  "polr_full"  = "/n/home_fasse/txia/UKB_CKM/Data/polr_full_selected_proteins_with_info.csv",
  "polr_top07" = "/n/home_fasse/txia/UKB_CKM/Data/polr_top_selected_proteins_with_info_coef0.7.csv"
)

for (nm in names(overall_score_files)) {
  coef_path <- overall_score_files[[nm]]
  coef_dt <- fread(coef_path)

  # keep only needed cols, drop intercept if present, drop zero coeff
  coef_dt <- coef_dt[, .(Protein, Coefficient)]
  coef_dt <- coef_dt[Protein != "(Intercept)"]
  coef_dt <- coef_dt[!is.na(Coefficient) & Coefficient != 0]

  proteins <- coef_dt$Protein
  coef_vec <- coef_dt$Coefficient
  names(coef_vec) <- proteins

  # Check missing proteins in ukb
  missing_p <- setdiff(proteins, colnames(ukb))
  if (length(missing_p) > 0) stop(paste("Missing proteins for", nm, ":", paste(missing_p, collapse = ", ")))

  # Score calculation: sum(beta * protein) ONLY
  X <- as.matrix(ukb[, ..proteins])
  score <- as.numeric(X %*% coef_vec)

  score_table[[paste0("manual_protein_score_", nm)]] <- score
}

# === Step 3: Save scores table ===
fwrite(score_table, "/n/home_fasse/txia/UKB_CKM/Data/manual_protein_scores_overall.csv")

# === Step 4: Merge scores back to UKB ===
ukb_combined <- merge(ukb, score_table, by = "f.eid", all.x = TRUE)
fwrite(ukb_combined, "/n/home_fasse/txia/UKB_CKM/Data/UKB_combined_with_overall_protein_scores.csv")

cat("✅ Overall protein scores calculated (sum beta*protein), saved, and merged into full UKB dataset.\n")

# Load the ckd dataset
ukb_data <- fread("/n/home_fasse/txia/UKB_CKM/Data/ckdvariable.csv")

# List of your protein score columns (UPDATED)
score_vars <- c(
  "manual_protein_score_polr_full",
  "manual_protein_score_polr_top07"
)

# Summary + SD for each
for (var in score_vars) {
  cat("-----", var, "-----\n")
  print(summary(ukb_combined[[var]]))
  cat("SD:", round(sd(ukb_combined[[var]], na.rm = TRUE), 4), "\n\n")
}

# Apply INT transformation directly to the same variable names (KEEP same behavior as your code)
for (var in score_vars) {
  ukb_combined[, (var) := qnorm(rank(get(var), na.last = "keep") / (.N + 1))]
}

# Summary + SD again (same as your script)
for (var in score_vars) {
  cat("-----", var, "-----\n")
  print(summary(ukb_combined[[var]]))
  cat("SD:", round(sd(ukb_combined[[var]], na.rm = TRUE), 4), "\n\n")
}

ckd_association <- ukb_combined

# Join the datasets by participant ID
ckd_association <- ckd_association %>%
  left_join(ukb_data, by = "f.eid")

# Define ckd-related outcomes
ckd_outcomes <- c("f.132032.0.0", "f.42026.0.0")

# Convert date fields to numeric
for (outcome in ckd_outcomes) {
  ckd_association[[outcome]] <- as.numeric(as.Date(ckd_association[[outcome]]))
}

# Define baseline visit date
ckd_association$dt_visit0 <- as.numeric(as.Date(ckd_association$f.21842.0.0))

# Identify first ckd diagnosis date and handle Inf values
ckd_association$ckd_dxdt <- apply(
  as.matrix(ckd_association[, ..ckd_outcomes]),
  1,
  function(x) {
    min_val <- min(x, na.rm = TRUE)
    if (is.infinite(min_val)) NA else min_val
  }
)

ckd_association$ckd_dxdt[is.na(ckd_association$ckd_dxdt)] <- Inf

# Determine baseline ckd
ckd_association$ckd_bsln <- ifelse(
  ckd_association$ckd_dxdt < ckd_association$dt_visit0,
  1, 0
)

# Exclude participants with baseline ckd
incident_data <- ckd_association %>% filter(ckd_bsln == 0)

# Define incident ckd and censoring
incident_data <- incident_data %>%
  mutate(
    incident_ckd_time = ckd_dxdt,
    event = ifelse(!is.na(incident_ckd_time) & incident_ckd_time != Inf, 1, 0),
    censor_time = pmin(as.numeric(as.Date(f.40000.0.0)),
                       as.numeric(as.Date(f.191.0.0)),
                       as.numeric(as.Date("2022-11-30")),
                       na.rm = TRUE),
    time_to_event = ifelse(event == 1, incident_ckd_time - dt_visit0, censor_time - dt_visit0)
  )

print("Baseline and incident ckd have been defined successfully.")
print(table(incident_data$event))

library(data.table)
library(dplyr)
library(survival)
library(splines)

# List of protein score variables (UPDATED)
protein_scores <- c(
  "manual_protein_score_polr_full",
  "manual_protein_score_polr_top07"
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

    fwrite(new_data, paste0("/n/home_fasse/txia/UKB_CKM/Data/spline_ckd_overallscore_alldata_", score_var, "_model3_", label, ".csv"))
  }
}

# === Apply models ===
run_models_on_subset(incident_data, "overall")
run_models_on_subset(subset(incident_data, bptreat.y == 0 & statin.y == 0), "nodrug")
run_models_on_subset(subset(incident_data, bptreat.y == 1 | statin.y == 1), "drug")

# === Adjust and export results ===
all_tertile_results$Bonferroni_P <- p.adjust(all_tertile_results$P_Value, method = "bonferroni")
all_tertile_results$FDR <- p.adjust(all_tertile_results$P_Value, method = "fdr")

fwrite(all_proteomics_score_results, "/n/home_fasse/txia/UKB_CKM/Data/ckd_alldata_all_overallscores_all_models.csv")
fwrite(all_tertile_results, "/n/home_fasse/txia/UKB_CKM/Data/ckd_alldata_all_overallscores_tertile_models.csv")

cat("✅ All models completed for all protein scores (continuous, tertile, spline).\n")


