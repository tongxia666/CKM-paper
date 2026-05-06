# ================================================================
# Proteomic score with ridge-adjusted weights (λ chosen to MAXIMIZE
# mean macro one-vs-rest AUC across CKM stages)
# - Stratified 5-fold CV by CKM stage
# - Macro one-vs-rest AUC (pROC) as CV metric
# - Robust to p=1; ordered-factor fixes
# - Compares λ choices on TEST; FINALIZES with λ_maxAUC
# - Builds parsimonious overall CKM score_top (|w_full| >= 0.7)
# Files tagged *_auc_max
# ================================================================

# -------------------------------
# 0) Libraries & lib paths
# -------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
myPaths <- .libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs", myPaths)
.libPaths(myPaths)

suppressPackageStartupMessages({
  library(dplyr,      lib.loc = "/n/home_fasse/txia/R/libs")
  library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
  library(MASS)       # polr()
  library(Hmisc,      lib.loc = "/n/home_fasse/txia/R/libs")  # rcorr.cens (optional)
  library(pROC,       lib.loc = "/n/home_fasse/txia/R/libs")  # AUC
})

# -------------------------------
# Helpers
# -------------------------------
as_ordered_like <- function(x, ref_levels) {
  factor(as.character(x), levels = ref_levels, ordered = TRUE)
}

ensure_matrix <- function(S, p_hint = NULL) {
  if (is.null(dim(S))) {
    p <- if (is.null(p_hint)) 1L else as.integer(p_hint)
    matrix(as.numeric(S), nrow = p, ncol = p)
  } else S
}

safe_diag <- function(val, p) {
  if (length(p) == 0L || p < 1L) stop("safe_diag: p must be >=1")
  vv <- if (length(val) == 0L || !is.finite(val)) 1e-6 else val
  diag(vv, p, p)
}

weights_for_lambda <- function(lambda, b, Sigma, jitterI = 1e-6) {
  p <- length(b)
  if (p < 1L) stop("weights_for_lambda: no proteins (p=0).")
  Sigma <- ensure_matrix(Sigma, p_hint = p)
  lam  <- if (length(lambda) == 0L || !is.finite(lambda)) 1e-6 else lambda
  w_ridge <- as.numeric(solve(Sigma + safe_diag(lam + jitterI, p), b))
  names(w_ridge) <- names(b)
  denom <- sum(w_ridge * w_ridge)
  scale_factor <- if (denom <= 0) 1 else sum(b * w_ridge) / denom
  w_ridge * scale_factor
}

# Stratified K-folds by y (ordered factor)
make_stratified_folds <- function(y, K = 5L, seed = 123L) {
  set.seed(seed)
  y_chr <- as.character(y)
  folds <- rep(NA_integer_, length(y_chr))
  for (lev in unique(y_chr)) {
    idx <- which(y_chr == lev)
    if (length(idx) == 0L) next
    assign <- rep(1:K, length.out = length(idx))
    folds[idx] <- sample(assign)  # shuffle within level
  }
  if (anyNA(folds)) {
    rem <- which(is.na(folds))
    folds[rem] <- sample(rep(1:K, length.out = length(rem)))
  }
  folds
}

# Pick λ: maximize mean metric; also compute 5% band (optional)
choose_lambdas_metric <- function(m, lam_grid) {
  ok <- which(is.finite(m))
  if (!length(ok)) {
    return(list(lam_max = 1e-3, lam_5pct = 1e-3, max_val = NA_real_))
  }
  idx_max <- ok[which.max(m[ok])]
  lam_max <- lam_grid[idx_max]
  M_max   <- m[idx_max]
  cand_5pct <- which(is.finite(m) & m >= 0.95 * M_max)
  lam_5pct  <- if (length(cand_5pct)) max(lam_grid[cand_5pct]) else lam_max
  list(lam_max = lam_max, lam_5pct = lam_5pct, max_val = M_max)
}

# Macro one-vs-rest AUC on (score, y_factor)
macro_auc_one_vs_rest <- function(score, y_fac) {
  levs <- levels(y_fac)
  aucs <- c()
  for (lv in levs) {
    yb <- as.numeric(y_fac == lv)
    if (length(unique(yb)) < 2L) next
    roc_obj <- try(pROC::roc(yb, as.numeric(score), quiet = TRUE), silent = TRUE)
    if (!inherits(roc_obj, "try-error")) {
      aucs <- c(aucs, as.numeric(pROC::auc(roc_obj)))
    }
  }
  if (!length(aucs)) return(NA_real_)
  mean(aucs, na.rm = TRUE)
}

# -------------------------------
# 1) Data Preparation
# -------------------------------
cat(">>> Loading data...\n")
data1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/combined_proteomics_pheno.csv")
data2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/covariate.csv")
ukb <- data1 %>% left_join(data2, by = "f.eid")
ukb <- ukb[!is.na(ckm)]

proteomics_columns <- colnames(ukb)[2:2924]
proteomics_columns <- proteomics_columns[proteomics_columns %in% names(ukb)]

missing_percentage <- sapply(ukb[, ..proteomics_columns],
                             function(x) sum(is.na(x)) / length(x) * 100)
cat("Missing percentage computed for", length(missing_percentage), "proteins.\n")
cat("Proteins with >80% missing:", sum(missing_percentage > 80), "\n")
proteomics_columns <- setdiff(proteomics_columns, "P48060.x")

cat(">>> Median imputation...\n")
ukb[, (proteomics_columns) := lapply(.SD, function(protein) {
  protein[is.na(protein)] <- median(protein, na.rm = TRUE)
  protein
}), .SDcols = proteomics_columns]

cat(">>> Rank-based INT transform...\n")
ukb[, (proteomics_columns) := lapply(.SD, function(protein) {
  qnorm(rank(as.numeric(protein)) / (length(protein) + 1), mean = 0, sd = 1)
}), .SDcols = proteomics_columns]

# Drop stage 4 (established CVD)
ukb <- ukb[ckm != 4 & !is.na(ckm)]

# -------------------------------
# 2) Selected protein lists (from stage-specific analyses)
# -------------------------------
cat(">>> Loading selected protein lists...\n")
df_1 <- fread("/n/home_fasse/txia/UKB_CKM/Data/1vs0_selected_proteins_with_info_bootstrap.csv") %>%
  mutate(comparison = "1vs0")
df_2 <- fread("/n/home_fasse/txia/UKB_CKM/Data/2vs1_selected_proteins_with_info_bootstrap.csv") %>%
  mutate(comparison = "2vs1")
df_3 <- fread("/n/home_fasse/txia/UKB_CKM/Data/3vs2_selected_proteins_with_info_bootstrap.csv") %>%
  mutate(comparison = "3vs2")
combined_all <- bind_rows(df_1, df_2, df_3)

prot_pool <- intersect(unique(combined_all$Protein), colnames(ukb))
stopifnot(length(prot_pool) > 0)

# -------------------------------
# 3) Train/Test split
# -------------------------------
cat(">>> Train/Test split (75/25)...\n")
ukb_filtered <- data.table::copy(ukb)
set.seed(123)
rows <- sample(seq_len(nrow(ukb_filtered)))
ukb_filtered <- ukb_filtered[rows, ]
train_size <- round(0.75 * nrow(ukb_filtered))
training <- ukb_filtered[1:train_size, ]
testing  <- ukb_filtered[(train_size + 1):nrow(ukb_filtered), ]

training$ckm <- factor(training$ckm, ordered = TRUE)
testing$ckm  <- factor(testing$ckm,  levels = levels(training$ckm), ordered = TRUE)

covariate_names <- c("age.y", "sex.y", "nonwhite", "fast_hrs","townsend", "PA", "alcohol_cat")
covariate_names <- covariate_names[covariate_names %in% names(training)]
for (nm in intersect(c("nonwhite", "fast_hrs","alcohol_cat"), names(training))) {
  training[[nm]] <- as.factor(training[[nm]])
  testing[[nm]]  <- factor(testing[[nm]], levels = levels(training[[nm]]))
}

prot_keep <- prot_pool[prot_pool %in% names(training)]
stopifnot(length(prot_keep) > 0)

Xtr <- as.matrix(training[, ..prot_keep])
Xte <- as.matrix(testing [, ..prot_keep])

# -------------------------------
# 4) Per-protein marginal betas (covariate-adjusted)
# -------------------------------
cat(">>> Estimating per-protein adjusted betas via polr()...\n")
get_beta_raw <- function(pn){
  df <- cbind.data.frame(
    y    = training$ckm,
    prot = Xtr[, pn],
    training[, ..covariate_names]
  )
  fit <- try(MASS::polr(y ~ prot + ., data = df, method = "logistic", Hess = TRUE),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  b <- try(coef(fit)[["prot"]], silent = TRUE)
  if (inherits(b, "try-error") || is.na(b)) NA_real_ else as.numeric(b)
}
b <- vapply(colnames(Xtr), get_beta_raw, numeric(1))
names(b) <- colnames(Xtr)
b[is.na(b)] <- 0
cat("Proteins with non-zero marginal beta:", sum(b != 0), "of", length(b), "\n")

# ================================================================
# 5) Cross-validation over λ, MAXIMIZING macro one-vs-rest AUC
# ================================================================
cat(">>> Cross-validation over lambda (metric: macro AUC)...\n")
p <- ncol(Xtr); stopifnot(p >= 1L)
Sigma <- ensure_matrix(cov(Xtr), p_hint = p)

K <- 5
fold_id <- make_stratified_folds(training$ckm, K = K, seed = 123)

lam_grid <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5,
              1, 2, 5, 10, 50, 100, 500, 1e3, 1e4)
nlam <- length(lam_grid)
jitterI <- 1e-6

cv_auc <- matrix(NA_real_, nrow = K, ncol = nlam,
                 dimnames = list(paste0("fold", 1:K), paste0("lam_", lam_grid)))

for (k in seq_len(K)) {
  trk <- which(fold_id != k)
  vak <- which(fold_id == k)

  X_trk <- Xtr[trk, , drop = FALSE]
  X_vak <- Xtr[vak, , drop = FALSE]
  y_trk <- training$ckm[trk]
  y_vak <- training$ckm[vak]
  C_trk <- training[trk, ..covariate_names]

  if (length(unique(y_vak)) < 2L) next
  p_trk <- ncol(X_trk); if (p_trk < 1L) next

  # recompute b and Sigma on inner-train
  get_beta_trk <- function(j){
    df <- cbind.data.frame(y = y_trk, prot = X_trk[, j], C_trk)
    fit <- try(MASS::polr(y ~ prot + ., data = df, method = "logistic", Hess = TRUE),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(0)
    bb <- try(coef(fit)[["prot"]], silent = TRUE)
    if (inherits(bb, "try-error") || is.na(bb)) 0 else as.numeric(bb)
  }
  b_trk   <- vapply(colnames(X_trk), get_beta_trk, numeric(1))
  Sigma_k <- ensure_matrix(cov(X_trk), p_hint = p_trk)

  for (i in seq_len(nlam)) {
    lam <- lam_grid[i]
    w_trk <- try(
      as.numeric(solve(Sigma_k + safe_diag(lam + jitterI, length(b_trk)), b_trk)),
      silent = TRUE
    )
    if (inherits(w_trk, "try-error")) next
    s_vak <- as.numeric(X_vak %*% w_trk)

    auc_val <- macro_auc_one_vs_rest(
      s_vak,
      factor(y_vak, levels = levels(training$ckm), ordered = TRUE)
    )
    if (is.finite(auc_val)) cv_auc[k, i] <- auc_val
  }
}

mean_auc <- colMeans(cv_auc, na.rm = TRUE)

choice <- choose_lambdas_metric(mean_auc, lam_grid)
lam_maxAUC <- choice$lam_max
lam_5pct   <- choice$lam_5pct

cat(sprintf("λ_maxAUC = %g (CV mean macro AUC ≈ %.6f)\n", lam_maxAUC, choice$max_val))
cat(sprintf("λ_5%%     = %g\n", lam_5pct))

cv_curve_df <- data.frame(lambda = lam_grid, cv_mean_macro_auc = mean_auc)
fwrite(cv_curve_df,
       "/n/home_fasse/txia/UKB_CKM/Data/cv_curve_by_lambda_auc_max.csv")

# -------------------------------
# 6) Compare λ choices on TEST
# -------------------------------
eval_at_lambda <- function(lambda_choice, label = "") {
  w <- weights_for_lambda(lambda_choice, b, Sigma, jitterI = 1e-6)
  s_tr <- as.numeric(Xtr %*% w)
  s_te <- as.numeric(Xte %*% w)

  # Cutpoints from TRAIN (score-only)
  cut_model_tmp <- MASS::polr(
    ckm ~ scale(s_tr),
    data = data.frame(ckm = training$ckm, s_tr = s_tr),
    method = "logistic", Hess = TRUE
  )
  pred_class <- predict(cut_model_tmp,
                        newdata = data.frame(s_tr = s_te),
                        type = "class")
  pred_class <- as_ordered_like(pred_class, levels(testing$ckm))
  acc <- mean(as.character(pred_class) == as.character(testing$ckm))

  # TEST macro AUC & ordinal C-index
  test_auc <- macro_auc_one_vs_rest(s_te, testing$ckm)
  cin <- try(Hmisc::rcorr.cens(s_te, as.numeric(testing$ckm)), silent = TRUE)
  cidx <- if (!inherits(cin, "try-error")) as.numeric(cin["C Index"]) else NA_real_

  i <- which(lam_grid == lambda_choice)
  data.frame(
    lambda            = lambda_choice,
    which             = label,
    cv_mean_macro_auc = if (length(i)) mean_auc[i] else NA_real_,
    test_macro_auc    = test_auc,
    test_c_index      = cidx,
    test_accuracy     = acc,
    stringsAsFactors  = FALSE
  )
}

cat(">>> Comparing lambda choices on TEST...\n")
lambda_comparison <- rbind(
  eval_at_lambda(lam_maxAUC, "lambda_maxAUC"),
  eval_at_lambda(lam_5pct,   "lambda_5pct")
)
print(lambda_comparison)
fwrite(lambda_comparison,
       "/n/home_fasse/txia/UKB_CKM/Data/lambda_comparison_summary_auc_max.csv")

# -------------------------------
# 7) FINALIZE with λ_maxAUC (FULL MODEL)
# -------------------------------
best_lambda <- lam_maxAUC
cat(sprintf(">>> Using λ = %g (MAX mean macro AUC) for final model.\n",
            best_lambda))
w_raw <- weights_for_lambda(best_lambda, b, Sigma, jitterI = 1e-6)

cat("Corr(w_raw, b) =", round(cor(w_raw, b), 3), "\n")
cat("||w_raw|| / ||b|| =",
    round(sqrt(sum(w_raw^2)) / max(1e-12, sqrt(sum(b^2))), 3), "\n")

training$manual_protein_score <- as.numeric(Xtr %*% w_raw)
testing$manual_protein_score  <- as.numeric(Xte %*% w_raw)

# Final TEST metrics (FULL)
ckm_num_test <- as.numeric(testing$ckm)
c_index_full <- Hmisc::rcorr.cens(testing$manual_protein_score, ckm_num_test)
cat("✅ (FULL) C-index (TEST):",
    round(c_index_full["C Index"], 4), "\n")
final_macro_auc_full <- macro_auc_one_vs_rest(
  testing$manual_protein_score, testing$ckm
)
cat("✅ (FULL) Macro AUC (TEST):",
    round(final_macro_auc_full, 4), "\n")

# Train cutpoints; classify TEST (FULL)
cut_model_full <- MASS::polr(
  ckm ~ scale(manual_protein_score),
  data = data.frame(ckm = training$ckm,
                    manual_protein_score = training$manual_protein_score),
  method = "logistic", Hess = TRUE
)
pred_class_full <- predict(
  cut_model_full,
  newdata = data.frame(
    manual_protein_score = testing$manual_protein_score
  ),
  type = "class"
)
testing$predicted_ckm_full <- as_ordered_like(pred_class_full,
                                              levels(testing$ckm))
acc_full <- mean(as.character(testing$predicted_ckm_full) ==
                 as.character(testing$ckm))
cat("✅ (FULL) Accuracy (TEST):", round(acc_full, 4), "\n")
cat("✅ Confusion Matrix (FULL):\n")
print(table(Predicted = testing$predicted_ckm_full,
            Observed  = testing$ckm))

# ================================================================
# 8) Parsimonious overall CKM score_top (|w_full| >= 0.6)
#    - Repeat same procedure as full model but restricted to top proteins
# ================================================================

cat("\n>>> Building parsimonious overall CKM score_top (|w| >= 0.6)...\n")

# 8.1 Select top proteins by absolute ridge-adjusted weight from full model
thr_top <- 0.6
top_prots <- names(w_raw)[abs(w_raw) >= thr_top]
cat("Number of top proteins (|w| >=", thr_top, "):", length(top_prots), "\n")
stopifnot(length(top_prots) > 0)

# Subset design matrices to top proteins
Xtr_top <- as.matrix(training[, ..top_prots])
Xte_top <- as.matrix(testing [, ..top_prots])

# 8.2 Re-estimate marginal betas (β_top) with covariate-adjusted ordinal models
cat(">>> Estimating per-protein adjusted betas for top proteins via polr()...\n")
get_beta_top <- function(pn) {
  df <- cbind.data.frame(
    y    = training$ckm,
    prot = Xtr_top[, pn],
    training[, ..covariate_names]
  )
  fit <- try(MASS::polr(y ~ prot + ., data = df, method = "logistic", Hess = TRUE),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(0)
  b <- try(coef(fit)[["prot"]], silent = TRUE)
  if (inherits(b, "try-error") || is.na(b)) 0 else as.numeric(b)
}
b_top <- vapply(colnames(Xtr_top), get_beta_top, numeric(1))
names(b_top) <- colnames(Xtr_top)
cat("Top proteins with non-zero marginal beta:", sum(b_top != 0), "of",
    length(b_top), "\n")

# 8.3 Covariance matrix for top proteins
p_top <- ncol(Xtr_top); stopifnot(p_top >= 1L)
Sigma_top <- ensure_matrix(cov(Xtr_top), p_hint = p_top)

# 8.4 Stratified 5-fold CV over λ (same grid, macro one-vs-rest AUC)
cat(">>> Cross-validation over lambda for score_top (metric: macro AUC)...\n")
K_top <- 5L
fold_id_top <- make_stratified_folds(training$ckm, K = K_top, seed = 456)

lam_grid_top <- lam_grid        # reuse the same grid as full model
nlam_top <- length(lam_grid_top)

cv_auc_top <- matrix(NA_real_, nrow = K_top, ncol = nlam_top,
                     dimnames = list(paste0("fold", 1:K_top),
                                     paste0("lam_", lam_grid_top)))

for (k in seq_len(K_top)) {
  trk <- which(fold_id_top != k)
  vak <- which(fold_id_top == k)

  X_trk <- Xtr_top[trk, , drop = FALSE]
  X_vak <- Xtr_top[vak, , drop = FALSE]
  y_trk <- training$ckm[trk]
  y_vak <- training$ckm[vak]
  C_trk <- training[trk, ..covariate_names]

  if (length(unique(y_vak)) < 2L) next
  p_trk <- ncol(X_trk); if (p_trk < 1L) next

  # recompute betas and Sigma within inner-train
  get_beta_trk_top <- function(j) {
    df <- cbind.data.frame(y = y_trk, prot = X_trk[, j], C_trk)
    fit <- try(MASS::polr(y ~ prot + ., data = df, method = "logistic", Hess = TRUE),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(0)
    bb <- try(coef(fit)[["prot"]], silent = TRUE)
    if (inherits(bb, "try-error") || is.na(bb)) 0 else as.numeric(bb)
  }
  b_trk_top <- vapply(colnames(X_trk), get_beta_trk_top, numeric(1))
  Sigma_k_top <- ensure_matrix(cov(X_trk), p_hint = p_trk)

  for (i in seq_len(nlam_top)) {
    lam <- lam_grid_top[i]
    w_trk_top <- try(
      as.numeric(solve(Sigma_k_top + safe_diag(lam + jitterI, length(b_trk_top)),
                       b_trk_top)),
      silent = TRUE
    )
    if (inherits(w_trk_top, "try-error")) next
    s_vak <- as.numeric(X_vak %*% w_trk_top)

    auc_val <- macro_auc_one_vs_rest(
      s_vak,
      factor(y_vak, levels = levels(training$ckm), ordered = TRUE)
    )
    if (is.finite(auc_val)) cv_auc_top[k, i] <- auc_val
  }
}

mean_auc_top <- colMeans(cv_auc_top, na.rm = TRUE)
choice_top   <- choose_lambdas_metric(mean_auc_top, lam_grid_top)
lam_maxAUC_top <- choice_top$lam_max
cat(sprintf("λ_maxAUC (top) = %g (CV mean macro AUC ≈ %.6f)\n",
            lam_maxAUC_top, choice_top$max_val))

# 8.5 Final ridge-adjusted weights and score_top
w_top <- weights_for_lambda(lam_maxAUC_top, b_top, Sigma_top, jitterI = 1e-6)
names(w_top) <- names(b_top)

training$manual_protein_score_top <- as.numeric(Xtr_top %*% w_top)
testing$manual_protein_score_top  <- as.numeric(Xte_top %*% w_top)

# 8.6 Final TEST metrics for score_top
ckm_num_test_top <- as.numeric(testing$ckm)
c_index_top <- Hmisc::rcorr.cens(testing$manual_protein_score_top,
                                 ckm_num_test_top)
cat("✅ (TOP) C-index (TEST):",
    round(c_index_top["C Index"], 4), "\n")
final_macro_auc_top <- macro_auc_one_vs_rest(
  testing$manual_protein_score_top, testing$ckm
)
cat("✅ (TOP) Macro AUC (TEST):",
    round(final_macro_auc_top, 4), "\n")

# Train cutpoints; classify TEST (TOP)
cut_model_top <- MASS::polr(
  ckm ~ scale(manual_protein_score_top),
  data = data.frame(ckm = training$ckm,
                    manual_protein_score_top = training$manual_protein_score_top),
  method = "logistic", Hess = TRUE
)
pred_class_top <- predict(
  cut_model_top,
  newdata = data.frame(
    manual_protein_score_top = testing$manual_protein_score_top
  ),
  type = "class"
)
testing$predicted_ckm_top <- as_ordered_like(pred_class_top,
                                             levels(testing$ckm))
acc_top <- mean(as.character(testing$predicted_ckm_top) ==
                as.character(testing$ckm))
cat("✅ (TOP) Accuracy (TEST):", round(acc_top, 4), "\n")
cat("✅ Confusion Matrix (TOP):\n")
print(table(Predicted = testing$predicted_ckm_top,
            Observed  = testing$ckm))

cat(">>> Done building parsimonious overall CKM score_top.\n")

# -------------------------------
# 9) Save outputs (FULL + TOP)
# -------------------------------
cat(">>> Saving outputs (FULL + TOP)...\n")

## 9a) Full-model coefficients + annotations
coef_full_df <- data.frame(
  Protein     = names(w_raw),
  Coefficient = as.numeric(w_raw),
  stringsAsFactors = FALSE
)
proteomics_info <- fread(
  "/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/Olink_info_11.10.2023.csv"
)
proteomics_info$UniProt <- paste0(proteomics_info$UniProt, ".x")
merged_full <- merge(coef_full_df, proteomics_info,
                     by.x = "Protein", by.y = "UniProt", all.x = TRUE)
write.csv(
  merged_full,
  "/n/home_fasse/txia/UKB_CKM/Data/polr_full_selected_proteins_with_info.csv",
  row.names = FALSE
)

## 9b) Full-model cutpoints
intercepts_full_df <- data.frame(
  Threshold = names(cut_model_full$zeta),
  Value     = as.numeric(cut_model_full$zeta)
)
write.csv(
  intercepts_full_df,
  "/n/home_fasse/txia/UKB_CKM/Data/polr_full_model_intercepts.csv",
  row.names = FALSE
)

## 9c) Top-model coefficients + annotations
coef_top_df <- data.frame(
  Protein     = names(w_top),
  Coefficient = as.numeric(w_top),
  stringsAsFactors = FALSE
)
merged_top <- merge(coef_top_df, proteomics_info,
                    by.x = "Protein", by.y = "UniProt", all.x = TRUE)
write.csv(
  merged_top,
  "/n/home_fasse/txia/UKB_CKM/Data/polr_top_selected_proteins_with_info_coef0.7.csv",
  row.names = FALSE
)

## 9d) Top-model cutpoints
intercepts_top_df <- data.frame(
  Threshold = names(cut_model_top$zeta),
  Value     = as.numeric(cut_model_top$zeta)
)
write.csv(
  intercepts_top_df,
  "/n/home_fasse/txia/UKB_CKM/Data/polr_top_model_intercepts_coef0.7.csv",
  row.names = FALSE
)

## 9e) Add FULL and TOP scores for all participants in ukb_filtered
Xall_full <- as.matrix(ukb_filtered[, ..prot_keep])
ukb_filtered$manual_protein_score     <- as.numeric(Xall_full %*% w_raw)

Xall_top <- as.matrix(ukb_filtered[, ..top_prots])
ukb_filtered$manual_protein_score_top <- as.numeric(Xall_top %*% w_top)

write.csv(
  ukb_filtered,
  "/n/home_fasse/txia/UKB_CKM/Data/ukb_filtered_with_polr_manual_protein_score.csv",
  row.names = FALSE
)

cat("\n>>> Done (FULL + TOP overall CKM scores saved).\n")


# ================================================================
# 9f) Quintiles & deciles of overall scores, by CKM stage
# ================================================================

cat(">>> Creating quintiles/deciles for FULL and TOP scores and tabulating by CKM...\n")

# Make sure CKM is a factor (ordered)
ukb_filtered$ckm <- factor(ukb_filtered$ckm, ordered = TRUE)

# ---- FULL score: quintiles & deciles ----
ukb_filtered$manual_protein_score_q5  <- Hmisc::cut2(ukb_filtered$manual_protein_score, g = 5)
ukb_filtered$manual_protein_score_d10 <- Hmisc::cut2(ukb_filtered$manual_protein_score, g = 10)

# Quintile table: counts and row % within CKM stage
tab_full_q5 <- ukb_filtered %>%
  group_by(ckm, manual_protein_score_q5) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(ckm) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

# Decile table: counts and row % within CKM stage
tab_full_d10 <- ukb_filtered %>%
  group_by(ckm, manual_protein_score_d10) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(ckm) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

# ---- TOP score: quintiles & deciles ----
ukb_filtered$manual_protein_score_top_q5  <- Hmisc::cut2(ukb_filtered$manual_protein_score_top, g = 5)
ukb_filtered$manual_protein_score_top_d10 <- Hmisc::cut2(ukb_filtered$manual_protein_score_top, g = 10)

tab_top_q5 <- ukb_filtered %>%
  group_by(ckm, manual_protein_score_top_q5) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(ckm) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

tab_top_d10 <- ukb_filtered %>%
  group_by(ckm, manual_protein_score_top_d10) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(ckm) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

# ---- Save summary tables ----
write.csv(tab_full_q5,
          "/n/home_fasse/txia/UKB_CKM/Data/overall_full_score_quintile_by_ckm.csv",
          row.names = FALSE)

write.csv(tab_full_d10,
          "/n/home_fasse/txia/UKB_CKM/Data/overall_full_score_decile_by_ckm.csv",
          row.names = FALSE)

write.csv(tab_top_q5,
          "/n/home_fasse/txia/UKB_CKM/Data/overall_top_score_quintile_by_ckm.csv",
          row.names = FALSE)

write.csv(tab_top_d10,
          "/n/home_fasse/txia/UKB_CKM/Data/overall_top_score_decile_by_ckm.csv",
          row.names = FALSE)

cat(">>> Saved full/top score quintile & decile distributions by CKM stage.\n")


