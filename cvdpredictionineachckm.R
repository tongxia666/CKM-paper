############################################################
## FULL STANDALONE SCRIPT
## CKM + FULL PROTEOMICS SCORE
## - Follow-up summary
## - Observed risk heatmaps for ALL yearly horizons
## - Overall observed risk heatmap (not by year)
## - Q5-Q1 risk contrast plots for ALL yearly horizons
## - Threshold heatmaps
## - Corrected risk-gradient figure
## - Within-stage 5-year reclassification
## - CKM on x-axis, quintile on y-axis
## - Tile labels show risk %, 95% CI, and Events/N
## - Wider rectangular cells for readable labels
## - Nature-like lighter heatmap palette
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
myPaths <- .libPaths()
myPaths <- c("/n/home_fasse/txia/R/libs", myPaths)
.libPaths(myPaths)

library(data.table,   lib.loc = "/n/home_fasse/txia/R/libs")
library(dplyr,        lib.loc = "/n/home_fasse/txia/R/libs")
library(survival,     lib.loc = "/n/home_fasse/txia/R/libs")
library(ggplot2,      lib.loc = "/n/home_fasse/txia/R/libs")
library(RColorBrewer, lib.loc = "/n/home_fasse/txia/R/libs")

############################################################
## 1. LOAD DATA
############################################################

cvd_association <- fread("/n/home_fasse/txia/UKB_CKM/Data/ukb_filtered_with_polr_manual_protein_score.csv")
ukb_data        <- fread("/n/home_fasse/txia/UKB_CKM/Data/cvdvariable.csv")

############################################################
## 2. BASIC PREP
############################################################

cvd_association <- cvd_association[!is.na(ckm), ]

cvd_association$sbp_use <- ifelse(
  !is.na(cvd_association$f.4080.0.0.x),
  cvd_association$f.4080.0.0.x,
  cvd_association$f.93.0.0.x
)

cvd_association$protein_score <- cvd_association$manual_protein_score
cvd_association[, protein_score := qnorm(rank(protein_score, na.last = "keep") / (.N + 1))]

cvd_association <- cvd_association %>%
  left_join(ukb_data, by = "f.eid")

############################################################
## 3. DEFINE INCIDENT CVD
############################################################

CVD_outcomes <- c(
  "f.131298.0.0","f.131300.0.0","f.131302.0.0","f.131304.0.0","f.131306.0.0",
  "f.131360.0.0","f.131362.0.0","f.131366.0.0","f.131368.0.0","f.131354.0.0",
  "f.131296.0.0","f.131364.0.0","f.131058.0.0","f.131056.0.0","f.131386.0.0",
  "f.131350.0.0","f.42000.0.0","f.42006.0.0","f.42008.0.0","f.42010.0.0"
)

for (outcome in CVD_outcomes) {
  cvd_association[[outcome]] <- as.numeric(as.Date(cvd_association[[outcome]]))
}

cvd_association$dt_visit0 <- as.numeric(as.Date(cvd_association$f.21842.0.0))

cvd_association$cvd_dxdt <- apply(
  as.matrix(cvd_association[, ..CVD_outcomes]),
  1,
  function(x) {
    min_val <- min(x, na.rm = TRUE)
    if (is.infinite(min_val)) NA else min_val
  }
)

cvd_association$cvd_dxdt[is.na(cvd_association$cvd_dxdt)] <- Inf

cvd_association$cvd_bsln <- ifelse(
  cvd_association$cvd_dxdt < cvd_association$dt_visit0,
  1, 0
)

incident_data <- cvd_association %>%
  filter(cvd_bsln == 0)

incident_data <- incident_data[ckm != 4]

incident_data <- incident_data %>%
  mutate(
    incident_cvd_time = cvd_dxdt,
    event = ifelse(!is.na(incident_cvd_time) & incident_cvd_time != Inf, 1, 0),
    censor_time = pmin(
      as.numeric(as.Date(f.40000.0.0)),
      as.numeric(as.Date(f.191.0.0)),
      as.numeric(as.Date("2022-11-30")),
      na.rm = TRUE
    ),
    time_to_event = ifelse(event == 1, incident_cvd_time - dt_visit0, censor_time - dt_visit0)
  )

cat("Incident dataset created.\n")
print(table(incident_data$event, useNA = "ifany"))

############################################################
## 4. ANALYSIS DATASET FOR OBSERVED-RISK ANALYSES
############################################################

dat <- as.data.table(incident_data)
dat$ckm <- factor(dat$ckm)

vars_needed <- c("time_to_event", "event", "ckm", "protein_score")
dat <- dat[, ..vars_needed]
dat <- dat[complete.cases(dat)]
dat <- as.data.frame(dat)

############################################################
## 5. FOLLOW-UP SUMMARY
############################################################

followup_years <- dat$time_to_event / 365.25

followup_summary <- data.frame(
  N = length(followup_years),
  Min_Years = min(followup_years, na.rm = TRUE),
  Median_Years = median(followup_years, na.rm = TRUE),
  P90_Years = quantile(followup_years, 0.90, na.rm = TRUE),
  Max_Years = max(followup_years, na.rm = TRUE)
)

write.csv(
  followup_summary,
  "/n/home_fasse/txia/UKB_CKM/Data/cvd_followup_summary_years.csv",
  row.names = FALSE
)

p_followup <- ggplot(
  data.frame(followup_years = followup_years),
  aes(x = followup_years)
) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  labs(
    title = "Distribution of follow-up time",
    x = "Follow-up time (years)",
    y = "Number of participants"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_followup_distribution_years.png",
  p_followup, width = 8, height = 5, dpi = 300, bg = "white"
)

############################################################
## 6. DEFINE ALL YEARLY HORIZONS
############################################################

max_followup_year <- floor(max(dat$time_to_event, na.rm = TRUE) / 365.25)
years_to_run <- seq(1, max_followup_year, by = 1)

cat("Running observed-risk analyses for yearly horizons:\n")
print(years_to_run)

############################################################
## 7. HELPER: HEATMAP STYLE
## Large tile numbers + CI moved to supplement table
############################################################

make_heatmap_plot <- function(df, fill_title, plot_title) {

  nature_cols <- colorRampPalette(
    brewer.pal(9, "YlGnBu")[2:7]
  )(100)

  ## stage-level N for x-axis labels: 0 (n=XX)
  stage_n <- df %>%
    dplyr::group_by(CKM_Stage) %>%
    dplyr::summarise(Stage_N = sum(N, na.rm = TRUE), .groups = "drop")

  stage_labels <- setNames(
    paste0(stage_n$CKM_Stage, " (n=", format(stage_n$Stage_N, big.mark = ","), ")"),
    stage_n$CKM_Stage
  )

  ## Keep tile text short and large
  df$tile_label <- paste0(
    sprintf("%.1f%%", 100 * df$Observed_Risk),
    "\n",
    df$Events, "/", df$N
  )

  ggplot(
    df,
    aes(x = factor(CKM_Stage), y = Quintile, fill = Observed_Risk)
  ) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(
      aes(label = tile_label),
      size = 11.0,
      fontface = "bold",
      color = "black",
      lineheight = 0.90
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_gradientn(
      colours = nature_cols,
      name = fill_title,
      labels = function(x) paste0(sprintf("%.1f", 100 * x), "%")
    ) +
    labs(
      title = plot_title,
      x = "CKM stage",
      y = "Proteomics score quintile"
    ) +
    theme_bw(base_size = 26) +
    theme(
      plot.title = element_text(face = "bold", size = 30, hjust = 0),
      axis.title = element_text(size = 28, face = "bold"),
      axis.text.x = element_text(size = 22, color = "black", face = "bold"),
      axis.text.y = element_text(size = 24, color = "black", face = "bold"),
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 20),
      panel.grid = element_blank(),
      plot.margin = margin(25, 90, 25, 25)
    )
}
############################################################
## 8. OBSERVED RISK BY STAGE AND QUINTILE
## Horizons: ALL YEARS
############################################################

all_risk_results <- data.frame()

for (yr in years_to_run) {
  horizon_days <- yr * 365.25
  risk_results_this_year <- data.frame()

  for (s in levels(dat$ckm)) {
    ds <- dat[dat$ckm == s, ]
    if (nrow(ds) < 50 || sum(ds$event == 1, na.rm = TRUE) < 10) next

    qcuts <- quantile(ds$protein_score, probs = seq(0, 1, 0.20), na.rm = TRUE)
    qcuts <- unique(qcuts)
    if (length(qcuts) < 6) next

    ds$Quintile <- cut(
      ds$protein_score,
      breaks = qcuts,
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    )

    for (q in levels(ds$Quintile)) {
      sub <- ds[ds$Quintile == q, ]
      if (nrow(sub) == 0) next

      km <- survfit(Surv(time_to_event, event) ~ 1, data = sub)
      ssum <- summary(km, times = horizon_days, extend = TRUE)

      surv  <- ssum$surv[1]
      lower <- ssum$lower[1]
      upper <- ssum$upper[1]

      obs_risk <- 1 - surv
      ci_lower <- 1 - upper
      ci_upper <- 1 - lower

      risk_results_this_year <- rbind(
        risk_results_this_year,
        data.frame(
          Horizon_Years = yr,
          CKM_Stage = s,
          Quintile = q,
          N = nrow(sub),
          Events = sum(sub$event == 1),
          Observed_Risk = obs_risk,
          CI_lower = ci_lower,
          CI_upper = ci_upper
        )
      )
    }
  }

  all_risk_results <- rbind(all_risk_results, risk_results_this_year)

  write.csv(
    risk_results_this_year,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/cvd_observed_", yr, "y_risk_by_stage_and_quintile.csv"),
    row.names = FALSE
  )
}

write.csv(
  all_risk_results,
  "/n/home_fasse/txia/UKB_CKM/Data/cvd_observed_risk_all_horizons_by_stage_and_quintile.csv",
  row.names = FALSE
)

############################################################
## 9. HEATMAPS FOR ALL YEARLY HORIZONS
############################################################

for (yr in years_to_run) {
  temp <- all_risk_results %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  temp$CKM_Stage <- factor(temp$CKM_Stage, levels = sort(unique(temp$CKM_Stage)))
  temp$Quintile  <- factor(temp$Quintile, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

  p_heatmap <- make_heatmap_plot(
    df = temp,
    fill_title = paste0(yr, "-year risk"),
    plot_title = paste0("Observed ", yr, "-year CVD risk by CKM stage and proteomics quintile")
  )

  ggsave(
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_heatmap_", yr, "y_risk_stage_by_quintile_fullscore.png"),
    p_heatmap,
    width = 18,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}


############################################################
## 10. OVERALL HEATMAP + SUPPLEMENTAL CI TABLE
############################################################

overall_horizon_days <- max_followup_year * 365.25
overall_risk_results <- data.frame()

for (s in levels(dat$ckm)) {
  ds <- dat[dat$ckm == s, ]
  if (nrow(ds) < 50 || sum(ds$event == 1, na.rm = TRUE) < 10) next

  qcuts <- quantile(ds$protein_score, probs = seq(0, 1, 0.20), na.rm = TRUE)
  qcuts <- unique(qcuts)
  if (length(qcuts) < 6) next

  ds$Quintile <- cut(
    ds$protein_score,
    breaks = qcuts,
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
  )

  for (q in levels(ds$Quintile)) {
    sub <- ds[ds$Quintile == q, ]
    if (nrow(sub) == 0) next

    km <- survfit(Surv(time_to_event, event) ~ 1, data = sub)
    ssum <- summary(km, times = overall_horizon_days, extend = TRUE)

    surv  <- ssum$surv[1]
    lower <- ssum$lower[1]
    upper <- ssum$upper[1]

    obs_risk <- 1 - surv
    ci_lower <- 1 - upper
    ci_upper <- 1 - lower

    overall_risk_results <- rbind(
      overall_risk_results,
      data.frame(
        CKM_Stage = s,
        Quintile = q,
        N = nrow(sub),
        Events = sum(sub$event == 1),
        Observed_Risk = obs_risk,
        CI_lower = ci_lower,
        CI_upper = ci_upper
      )
    )
  }
}

overall_risk_results$CKM_Stage <- factor(
  overall_risk_results$CKM_Stage,
  levels = sort(unique(overall_risk_results$CKM_Stage))
)

overall_risk_results$Quintile <- factor(
  overall_risk_results$Quintile,
  levels = c("Q1", "Q2", "Q3", "Q4", "Q5")
)

write.csv(
  overall_risk_results,
  "/n/home_fasse/txia/UKB_CKM/Data/cvd_observed_overall_risk_by_stage_and_quintile.csv",
  row.names = FALSE
)

############################################################
## Supplemental CI table
############################################################

overall_ci_table <- overall_risk_results %>%
  dplyr::mutate(
    CKM_Stage_Label = paste0(
      CKM_Stage,
      " (n=",
      format(ave(N, CKM_Stage, FUN = sum), big.mark = ","),
      ")"
    ),
    Risk_Percent = sprintf("%.1f%%", 100 * Observed_Risk),
    CI_95 = paste0(
      sprintf("%.1f", 100 * CI_lower),
      "–",
      sprintf("%.1f", 100 * CI_upper),
      "%"
    ),
    Cases_N = paste0(Events, "/", N)
  ) %>%
  dplyr::select(
    CKM_Stage_Label,
    Quintile,
    Risk_Percent,
    CI_95,
    Cases_N
  )

write.csv(
  overall_ci_table,
  "/n/home_fasse/txia/UKB_CKM/Data/Supplemental_Table_overall_heatmap_95CI.csv",
  row.names = FALSE
)

############################################################
## Overall heatmap with large numbers
############################################################

p_heatmap_overall <- make_heatmap_plot(
  df = overall_risk_results,
  fill_title = "Overall risk",
  plot_title = paste0(
    "Overall cumulative CVD risk by CKM stage and proteomics quintile\n",
    "(through ", max_followup_year, "-year follow-up)"
  )
)

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_heatmap_overall_risk_stage_by_quintile_fullscore.png",
  p_heatmap_overall,
  width = 15,
  height = 10,
  dpi = 300,
  bg = "white"
)

############################################################
## 11. OPTIONAL: PERSON-YEARS TABLE
## Supplementary table, not the heatmap label
############################################################

person_year_summary <- all_risk_results
person_year_summary$Person_Years <- NA_real_
person_year_summary$Incidence_Rate_per_1000_PY <- NA_real_

for (i in seq_len(nrow(person_year_summary))) {
  yr <- person_year_summary$Horizon_Years[i]
  s  <- as.character(person_year_summary$CKM_Stage[i])
  q  <- as.character(person_year_summary$Quintile[i])

  ds <- dat[dat$ckm == s, ]
  qcuts <- quantile(ds$protein_score, probs = seq(0, 1, 0.20), na.rm = TRUE)
  qcuts <- unique(qcuts)
  if (length(qcuts) < 6) next

  ds$Quintile <- cut(
    ds$protein_score,
    breaks = qcuts,
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
  )

  sub <- ds[ds$Quintile == q, ]
  if (nrow(sub) == 0) next

  horizon_days <- yr * 365.25
  followup_used <- pmin(sub$time_to_event, horizon_days)
  py <- sum(followup_used, na.rm = TRUE) / 365.25

  person_year_summary$Person_Years[i] <- py
  person_year_summary$Incidence_Rate_per_1000_PY[i] <-
    ifelse(py > 0, 1000 * person_year_summary$Events[i] / py, NA)
}

write.csv(
  person_year_summary,
  "/n/home_fasse/txia/UKB_CKM/Data/cvd_observed_risk_with_person_years_by_stage_and_quintile.csv",
  row.names = FALSE
)

############################################################
## 12. Q5 - Q1 RISK CONTRASTS FOR ALL HORIZONS
############################################################

all_risk_contrasts <- data.frame()

for (yr in years_to_run) {
  temp <- all_risk_results %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  risk_wide <- reshape(
    temp[, c("CKM_Stage", "Quintile", "Observed_Risk")],
    idvar = "CKM_Stage",
    timevar = "Quintile",
    direction = "wide"
  )

  colnames(risk_wide) <- gsub("Observed_Risk.", "", colnames(risk_wide), fixed = TRUE)

  if (!all(c("Q1", "Q5") %in% colnames(risk_wide))) next

  risk_wide$Horizon_Years <- yr
  risk_wide$Risk_Difference_Q5_minus_Q1 <- risk_wide$Q5 - risk_wide$Q1
  risk_wide$Risk_Ratio_Q5_div_Q1 <- risk_wide$Q5 / risk_wide$Q1

  all_risk_contrasts <- rbind(all_risk_contrasts, risk_wide)

  write.csv(
    risk_wide,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/cvd_stage_specific_", yr, "y_risk_contrast_Q5_vs_Q1.csv"),
    row.names = FALSE
  )
}

write.csv(
  all_risk_contrasts,
  "/n/home_fasse/txia/UKB_CKM/Data/cvd_stage_specific_risk_contrast_all_horizons.csv",
  row.names = FALSE
)

############################################################
## 13. BAR PLOTS OF Q5 - Q1 DIFFERENCE FOR ALL HORIZONS
############################################################

for (yr in years_to_run) {
  temp <- all_risk_contrasts %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  temp$CKM_Stage <- factor(temp$CKM_Stage, levels = temp$CKM_Stage)

  p_bar <- ggplot(temp, aes(x = CKM_Stage, y = Risk_Difference_Q5_minus_Q1)) +
    geom_col(fill = "#5AAE61") +
    labs(
      title = paste0(yr, "-year CVD risk difference between highest and lowest proteomics quintiles"),
      x = "CKM stage",
      y = "Observed risk difference (Q5 - Q1)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )

  ggsave(
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_stage_specific_", yr, "y_risk_difference_fullscore.png"),
    p_bar, width = 8, height = 5, dpi = 300, bg = "white"
  )
}

############################################################
## 14. THRESHOLD ANALYSIS
############################################################

threshold_years <- c(5, 10)
threshold_years <- threshold_years[threshold_years %in% years_to_run]

threshold_results <- data.frame()

categorize_risk <- function(x) {
  ifelse(x < 0.05, "<5%",
         ifelse(x < 0.10, "5-<10%", ">=10%"))
}

for (yr in threshold_years) {
  temp <- all_risk_results %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  temp$Risk_Category <- categorize_risk(temp$Observed_Risk)
  threshold_results <- rbind(threshold_results, temp)
}

write.csv(
  threshold_results,
  "/n/home_fasse/txia/UKB_CKM/Data/cvd_threshold_categories_5y_10y.csv",
  row.names = FALSE
)

############################################################
## 15. THRESHOLD HEATMAPS
############################################################

for (yr in threshold_years) {
  temp <- threshold_results %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  temp$tile_label <- paste0(
    sprintf("%.1f%%", 100 * temp$Observed_Risk),
    "\n95% CI ",
    sprintf("%.1f", 100 * temp$CI_lower), "–",
    sprintf("%.1f", 100 * temp$CI_upper), "%",
    "\n",
    temp$Events, "/", temp$N
  )

  p_thr <- ggplot(
    temp,
    aes(x = factor(CKM_Stage), y = Quintile, fill = Risk_Category)
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = tile_label),
      size = 3.8,
      color = "black",
      lineheight = 0.92
    ) +
    labs(
      title = paste0("Clinically actionable ", yr, "-year CVD risk categories"),
      x = "CKM stage",
      y = "Proteomics score quintile",
      fill = "Risk category"
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 22),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      panel.grid = element_blank(),
      plot.margin = margin(20, 50, 20, 20)
    )

  ggsave(
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_threshold_heatmap_", yr, "y_fullscore.png"),
    p_thr, width = 18, height = 8, dpi = 300, bg = "white"
  )
}

############################################################
## 16. STACKED BAR PLOTS OF RISK CATEGORIES
############################################################

summary_threshold <- threshold_results %>%
  group_by(Horizon_Years, CKM_Stage, Risk_Category) %>%
  summarise(n_groups = n(), .groups = "drop")

for (yr in threshold_years) {
  temp <- summary_threshold %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  p_stack <- ggplot(temp, aes(x = factor(CKM_Stage), y = n_groups, fill = Risk_Category)) +
    geom_col(position = "stack") +
    labs(
      title = paste0("Risk-category distribution across CKM stages (", yr, "-year horizon)"),
      x = "CKM stage",
      y = "Number of proteomics quintile groups",
      fill = "Risk category"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )

  ggsave(
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_threshold_stacked_", yr, "y_fullscore.png"),
    p_stack, width = 8, height = 5, dpi = 300, bg = "white"
  )
}

############################################################
## 17. COMBINED LINE PLOT ACROSS ALL HORIZONS
############################################################

plot_contrast <- all_risk_contrasts
plot_contrast$CKM_Stage <- factor(
  plot_contrast$CKM_Stage,
  levels = sort(unique(plot_contrast$CKM_Stage))
)

p_line <- ggplot(
  plot_contrast,
  aes(
    x = Horizon_Years,
    y = Risk_Difference_Q5_minus_Q1,
    color = CKM_Stage,
    group = CKM_Stage
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = years_to_run) +
  labs(
    title = "Risk stratification by proteomics quintiles across follow-up horizons",
    x = "Prediction horizon (years)",
    y = "Observed risk difference (Q5 - Q1)",
    color = "CKM stage"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_risk_difference_across_horizons_fullscore.png",
  p_line, width = 10, height = 6, dpi = 300, bg = "white"
)

############################################################
## 18. CORRECTED RISK-GRADIENT FIGURE
############################################################

gradient_year <- if (5 %in% years_to_run) 5 else min(years_to_run)

risk_data <- read.csv(
  paste0("/n/home_fasse/txia/UKB_CKM/Data/cvd_observed_", gradient_year, "y_risk_by_stage_and_quintile.csv")
)

risk_data$CKM_Stage <- factor(risk_data$CKM_Stage)
risk_data$Quintile <- factor(risk_data$Quintile, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

p_gradient <- ggplot(
  risk_data,
  aes(
    x = CKM_Stage,
    y = Observed_Risk,
    color = Quintile,
    group = Quintile
  )
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.2) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(
    title = paste0("Observed ", gradient_year, "-year CVD risk across CKM stage and proteomics quintile"),
    x = "CKM stage",
    y = paste0("Observed ", gradient_year, "-year CVD risk"),
    color = "Proteomics score quintile"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_risk_gradient_CKM_proteomics_", gradient_year, "y.png"),
  p_gradient,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

############################################################
## 19. WITHIN-STAGE 5-YEAR RECLASSIFICATION
############################################################

reclass_dat <- as.data.table(incident_data)

reclass_dat$sex.y       <- as.factor(reclass_dat$sex.y)
reclass_dat$nonwhite    <- as.factor(reclass_dat$nonwhite)
reclass_dat$fast_hrs    <- as.factor(reclass_dat$fast_hrs)
reclass_dat$alcohol_cat <- as.factor(reclass_dat$alcohol_cat)
reclass_dat$ckm         <- factor(reclass_dat$ckm)

vars_reclass <- c(
  "time_to_event", "event", "ckm", "protein_score",
  "age.y", "sex.y", "nonwhite", "fast_hrs", "townsend", "PA", "alcohol_cat"
)

reclass_dat <- reclass_dat[, ..vars_reclass]
reclass_dat <- reclass_dat[complete.cases(reclass_dat)]
reclass_dat <- as.data.frame(reclass_dat)
reclass_dat$protein_score_z <- as.numeric(scale(reclass_dat$protein_score))

predict_cox_risk <- function(model, newdata, horizon_days) {
  bh <- basehaz(model, centered = FALSE)
  H0_t <- approx(
    x = bh$time,
    y = bh$hazard,
    xout = horizon_days,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y

  lp <- predict(model, newdata = newdata, type = "lp")
  risk <- 1 - exp(-H0_t * exp(lp))
  as.numeric(risk)
}

cut_risk_cat <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 0.05, 0.10, Inf),
    labels = c("<5%", "5-<10%", ">=10%"),
    right = FALSE
  )
}

reclass_individual <- data.frame()
reclass_summary <- data.frame()
reclass_counts_plot <- data.frame()

horizon_days <- 5 * 365.25

for (s in levels(reclass_dat$ckm)) {

  ds <- subset(reclass_dat, ckm == s)

  if (nrow(ds) < 100 || sum(ds$event == 1, na.rm = TRUE) < 20) next

  fit_base <- coxph(
    Surv(time_to_event, event) ~ age.y + sex.y + nonwhite + fast_hrs + townsend + PA + alcohol_cat,
    data = ds,
    x = TRUE, y = TRUE
  )

  fit_prot <- coxph(
    Surv(time_to_event, event) ~ age.y + sex.y + nonwhite + fast_hrs + townsend + PA + alcohol_cat + protein_score_z,
    data = ds,
    x = TRUE, y = TRUE
  )

  ds$risk_base_5y <- predict_cox_risk(fit_base, ds, horizon_days)
  ds$risk_prot_5y <- predict_cox_risk(fit_prot, ds, horizon_days)

  ds$cat_base_5y <- cut_risk_cat(ds$risk_base_5y)
  ds$cat_prot_5y <- cut_risk_cat(ds$risk_prot_5y)

  ds$movement <- ifelse(
    as.numeric(ds$cat_prot_5y) > as.numeric(ds$cat_base_5y), "Up",
    ifelse(as.numeric(ds$cat_prot_5y) < as.numeric(ds$cat_base_5y), "Down", "Same")
  )

  ds$CKM_Stage <- s

  reclass_individual <- rbind(
    reclass_individual,
    ds[, c("CKM_Stage", "time_to_event", "event",
           "risk_base_5y", "risk_prot_5y",
           "cat_base_5y", "cat_prot_5y", "movement")]
  )

  tab <- as.data.frame.matrix(table(ds$cat_base_5y, ds$cat_prot_5y))
  tab$Base_Category <- rownames(tab)
  tab$CKM_Stage <- s
  rownames(tab) <- NULL

  write.csv(
    tab,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/reclassification_table_stage_", s, "_5y.csv"),
    row.names = FALSE
  )

  n_total <- nrow(ds)
  n_up <- sum(ds$movement == "Up", na.rm = TRUE)
  n_down <- sum(ds$movement == "Down", na.rm = TRUE)
  n_same <- sum(ds$movement == "Same", na.rm = TRUE)

  n_up_event <- sum(ds$movement == "Up" & ds$event == 1, na.rm = TRUE)
  n_down_event <- sum(ds$movement == "Down" & ds$event == 1, na.rm = TRUE)
  n_up_nonevent <- sum(ds$movement == "Up" & ds$event == 0, na.rm = TRUE)
  n_down_nonevent <- sum(ds$movement == "Down" & ds$event == 0, na.rm = TRUE)

  n_event <- sum(ds$event == 1, na.rm = TRUE)
  n_nonevent <- sum(ds$event == 0, na.rm = TRUE)

  nri_events <- ifelse(n_event > 0, (n_up_event - n_down_event) / n_event, NA)
  nri_nonevents <- ifelse(n_nonevent > 0, (n_down_nonevent - n_up_nonevent) / n_nonevent, NA)
  nri_total <- nri_events + nri_nonevents

  reclass_summary <- rbind(
    reclass_summary,
    data.frame(
      CKM_Stage = s,
      N = n_total,
      Events = n_event,
      NonEvents = n_nonevent,
      Up = n_up,
      Down = n_down,
      Same = n_same,
      Pct_Up = 100 * n_up / n_total,
      Pct_Down = 100 * n_down / n_total,
      Pct_Same = 100 * n_same / n_total,
      NRI_Events = nri_events,
      NRI_NonEvents = nri_nonevents,
      NRI_Total = nri_total
    )
  )

  reclass_counts_plot <- rbind(
    reclass_counts_plot,
    data.frame(CKM_Stage = s, Movement = "Up", Count = n_up),
    data.frame(CKM_Stage = s, Movement = "Same", Count = n_same),
    data.frame(CKM_Stage = s, Movement = "Down", Count = n_down)
  )
}

write.csv(
  reclass_individual,
  "/n/home_fasse/txia/UKB_CKM/Data/reclassification_individual_5y_within_stage.csv",
  row.names = FALSE
)

write.csv(
  reclass_summary,
  "/n/home_fasse/txia/UKB_CKM/Data/reclassification_summary_5y_within_stage.csv",
  row.names = FALSE
)

############################################################
## 20. RECLASSIFICATION FIGURES
############################################################

reclass_counts_plot$CKM_Stage <- factor(reclass_counts_plot$CKM_Stage)

p_reclass_counts <- ggplot(
  reclass_counts_plot,
  aes(x = CKM_Stage, y = Count, fill = Movement)
) +
  geom_col(position = "stack") +
  labs(
    title = "Within-stage 5-year risk reclassification after adding proteomics",
    x = "CKM stage",
    y = "Number of individuals",
    fill = "Risk category movement"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_reclassification_counts_5y_within_stage.png",
  p_reclass_counts,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

reclass_pct_plot <- reclass_summary[, c("CKM_Stage", "Pct_Up", "Pct_Same", "Pct_Down")]
reclass_pct_plot <- reshape(
  reclass_pct_plot,
  varying = c("Pct_Up", "Pct_Same", "Pct_Down"),
  v.names = "Percent",
  timevar = "Movement",
  times = c("Up", "Same", "Down"),
  direction = "long"
)

reclass_pct_plot$CKM_Stage <- factor(reclass_pct_plot$CKM_Stage)
reclass_pct_plot$Movement <- factor(reclass_pct_plot$Movement, levels = c("Down", "Same", "Up"))

p_reclass_pct <- ggplot(
  reclass_pct_plot,
  aes(x = CKM_Stage, y = Percent, fill = Movement)
) +
  geom_col(position = "stack") +
  labs(
    title = "Within-stage 5-year risk reclassification after adding proteomics",
    subtitle = "Percent moving to lower, same, or higher risk category",
    x = "CKM stage",
    y = "Percent of individuals",
    fill = "Risk category movement"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_reclassification_percent_5y_within_stage.png",
  p_reclass_pct,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

############################################################
## 21. OPTIONAL BUBBLE PLOT FOR CKM STAGE 3
############################################################

if ("3" %in% reclass_summary$CKM_Stage) {
  ds3 <- subset(reclass_individual, CKM_Stage == "3")
  tab3 <- as.data.frame(table(ds3$cat_base_5y, ds3$cat_prot_5y))
  colnames(tab3) <- c("Base", "With_Proteomics", "Count")

  p_bubble3 <- ggplot(tab3, aes(x = Base, y = With_Proteomics, size = Count)) +
    geom_point(alpha = 0.7) +
    labs(
      title = "Reclassification within CKM stage 3",
      subtitle = "5-year risk categories before vs after adding proteomics",
      x = "Base clinical model",
      y = "Clinical + proteomics model"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )

  ggsave(
    "/n/home_fasse/txia/UKB_CKM/Data/Figure_reclassification_bubble_stage3_5y.png",
    p_bubble3,
    width = 7,
    height = 5.5,
    dpi = 300,
    bg = "white"
  )
}

############################################################
## 22. OUTPUT SUMMARY
############################################################

cat("\nDone.\n")
cat("Saved key files:\n")
cat(" - cvd_followup_summary_years.csv\n")
cat(" - Figure_followup_distribution_years.png\n")
cat(" - cvd_observed_[1...max]y_risk_by_stage_and_quintile.csv\n")
cat(" - cvd_observed_risk_all_horizons_by_stage_and_quintile.csv\n")
cat(" - Figure_heatmap_[1...max]y_risk_stage_by_quintile_fullscore.png\n")
cat(" - cvd_observed_overall_risk_by_stage_and_quintile.csv\n")
cat(" - Figure_heatmap_overall_risk_stage_by_quintile_fullscore.png\n")
cat(" - cvd_observed_risk_with_person_years_by_stage_and_quintile.csv\n")
cat(" - cvd_stage_specific_[1...max]y_risk_contrast_Q5_vs_Q1.csv\n")
cat(" - cvd_stage_specific_risk_contrast_all_horizons.csv\n")
cat(" - Figure_stage_specific_[1...max]y_risk_difference_fullscore.png\n")
cat(" - cvd_threshold_categories_5y_10y.csv\n")
cat(" - Figure_threshold_heatmap_5y_fullscore.png\n")
cat(" - Figure_threshold_heatmap_10y_fullscore.png\n")
cat(" - Figure_threshold_stacked_5y_fullscore.png\n")
cat(" - Figure_threshold_stacked_10y_fullscore.png\n")
cat(" - Figure_risk_difference_across_horizons_fullscore.png\n")
cat(" - Figure_risk_gradient_CKM_proteomics_[chosen year].png\n")
cat(" - reclassification_individual_5y_within_stage.csv\n")
cat(" - reclassification_summary_5y_within_stage.csv\n")
cat(" - reclassification_table_stage_0_5y.csv, etc.\n")
cat(" - Figure_reclassification_counts_5y_within_stage.png\n")
cat(" - Figure_reclassification_percent_5y_within_stage.png\n")
cat(" - Figure_reclassification_bubble_stage3_5y.png (if stage 3 available)\n")

