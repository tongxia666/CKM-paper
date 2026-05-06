############################################################
## FULL STANDALONE SCRIPT: T2D RISK HEATMAP WORKFLOW
## - Follow-up summary
## - Observed risk heatmaps for ALL yearly horizons
## - Overall observed risk heatmap (not by year)
## - Q5-Q1 risk contrast plots for ALL yearly horizons
## - Threshold heatmaps
## - CKM on x-axis, quintile on y-axis
## - Tile labels show risk %, 95% CI, and Events/N
## - Wider rectangular cells for readable labels
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

t2d_association <- fread("/n/home_fasse/txia/UKB_CKM/Data/ukb_filtered_with_polr_manual_protein_score.csv")
ukb_data        <- fread("/n/home_fasse/txia/UKB_CKM/Data/t2dvariable.csv")

############################################################
## 2. BASIC PREP
############################################################

t2d_association <- t2d_association[!is.na(ckm), ]

t2d_association$sbp_use <- ifelse(
  !is.na(t2d_association$f.4080.0.0.x),
  t2d_association$f.4080.0.0.x,
  t2d_association$f.93.0.0.x
)

t2d_association$protein_score <- t2d_association$manual_protein_score
t2d_association[, protein_score := qnorm(rank(protein_score, na.last = "keep") / (.N + 1))]

t2d_association <- t2d_association %>%
  left_join(ukb_data, by = "f.eid")

############################################################
## 3. DEFINE INCIDENT T2D
############################################################

t2d_outcomes <- c("f.130708.0.0")

for (outcome in t2d_outcomes) {
  t2d_association[[outcome]] <- as.numeric(as.Date(t2d_association[[outcome]]))
}

t2d_association$dt_visit0 <- as.numeric(as.Date(t2d_association$f.21842.0.0))

t2d_association$t2d_dxdt <- apply(
  as.matrix(t2d_association[, ..t2d_outcomes]),
  1,
  function(x) {
    min_val <- min(x, na.rm = TRUE)
    if (is.infinite(min_val)) NA else min_val
  }
)

t2d_association$t2d_dxdt[is.na(t2d_association$t2d_dxdt)] <- Inf

t2d_association$t2d_bsln <- ifelse(
  t2d_association$t2d_dxdt < t2d_association$dt_visit0,
  1, 0
)

incident_data <- t2d_association %>%
  filter(t2d_bsln == 0)

incident_data <- incident_data[ckm != 4]

incident_data <- incident_data %>%
  mutate(
    incident_t2d_time = t2d_dxdt,
    event = ifelse(!is.na(incident_t2d_time) & incident_t2d_time != Inf, 1, 0),
    censor_time = pmin(
      as.numeric(as.Date(f.40000.0.0)),
      as.numeric(as.Date(f.191.0.0)),
      as.numeric(as.Date("2022-11-30")),
      na.rm = TRUE
    ),
    time_to_event = ifelse(event == 1, incident_t2d_time - dt_visit0, censor_time - dt_visit0)
  )

cat("Incident T2D dataset created.\n")
print(table(incident_data$event, useNA = "ifany"))

############################################################
## 4. ANALYSIS DATASET
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
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_followup_summary_years.csv",
  row.names = FALSE
)

p_followup <- ggplot(data.frame(followup_years = followup_years),
                     aes(x = followup_years)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  labs(
    title = "Distribution of follow-up time for T2D",
    x = "Follow-up time (years)",
    y = "Number of participants"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_t2d_followup_distribution_years.png",
  p_followup, width = 8, height = 5, dpi = 300, bg = "white"
)

############################################################
## 6. DEFINE ALL YEARLY HORIZONS
############################################################

max_followup_year <- floor(max(dat$time_to_event, na.rm = TRUE) / 365.25)
years_to_run <- seq(1, max_followup_year, by = 1)

############################################################
## 7. HELPER: HEATMAP STYLE
## Large tile numbers + CI moved to supplement table
############################################################

make_heatmap_plot <- function(df, fill_title, plot_title) {

  nature_cols <- colorRampPalette(
    brewer.pal(9, "YlGnBu")[2:7]
  )(100)

  df$CKM_Stage <- factor(df$CKM_Stage, levels = sort(unique(df$CKM_Stage)))
  df$Quintile  <- factor(df$Quintile, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

  ## Large tile label: risk % and cases/N only
  df$tile_label <- paste0(
    sprintf("%.1f%%", 100 * df$Observed_Risk),
    "\n",
    df$Events, "/", df$N
  )

  ggplot(
    df,
    aes(x = CKM_Stage, y = Quintile, fill = Observed_Risk)
  ) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(
      aes(label = tile_label),
      size = 11.0,
      fontface = "bold",
      color = "black",
      lineheight = 0.90
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
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
    paste0("/n/home_fasse/txia/UKB_CKM/Data/t2d_observed_", yr, "y_risk_by_stage_and_quintile.csv"),
    row.names = FALSE
  )
}

write.csv(
  all_risk_results,
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_observed_risk_all_horizons_by_stage_and_quintile.csv",
  row.names = FALSE
)

############################################################
## 9. YEARLY HEATMAPS + YEARLY CI TABLES
############################################################

for (yr in years_to_run) {
  temp <- all_risk_results %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  temp$CKM_Stage <- factor(temp$CKM_Stage, levels = sort(unique(temp$CKM_Stage)))
  temp$Quintile  <- factor(temp$Quintile, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

  ## Supplemental 95% CI table for each yearly heatmap
  temp_ci_table <- temp %>%
    dplyr::mutate(
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
      Horizon_Years,
      CKM_Stage,
      Quintile,
      Risk_Percent,
      CI_95,
      Cases_N
    )

  write.csv(
    temp_ci_table,
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Supplemental_Table_t2d_heatmap_", yr, "y_95CI.csv"),
    row.names = FALSE
  )

  p_heatmap <- make_heatmap_plot(
    df = temp,
    fill_title = paste0(yr, "-year risk"),
    plot_title = paste0("Observed ", yr, "-year T2D risk by CKM stage and proteomics quintile")
  )

  ggsave(
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_t2d_heatmap_", yr, "y_risk_stage_by_quintile_fullscore.png"),
    p_heatmap,
    width = 15,
    height = 10,
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
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_observed_overall_risk_by_stage_and_quintile.csv",
  row.names = FALSE
)

overall_ci_table <- overall_risk_results %>%
  dplyr::mutate(
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
    CKM_Stage,
    Quintile,
    Risk_Percent,
    CI_95,
    Cases_N
  )

write.csv(
  overall_ci_table,
  "/n/home_fasse/txia/UKB_CKM/Data/Supplemental_Table_t2d_overall_heatmap_95CI.csv",
  row.names = FALSE
)

p_heatmap_overall <- make_heatmap_plot(
  df = overall_risk_results,
  fill_title = "Overall risk",
  plot_title = paste0(
    "Overall cumulative T2D risk by CKM stage and proteomics quintile\n",
    "(through ", max_followup_year, "-year follow-up)"
  )
)

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_t2d_heatmap_overall_risk_stage_by_quintile_fullscore.png",
  p_heatmap_overall,
  width = 15,
  height = 10,
  dpi = 300,
  bg = "white"
)
############################################################
## 11. PERSON-YEARS TABLE
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

  ds$Quintile <- cut(ds$protein_score,
                     breaks = qcuts,
                     include.lowest = TRUE,
                     labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

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
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_observed_risk_with_person_years_by_stage_and_quintile.csv",
  row.names = FALSE
)

############################################################
## 12. Q5 - Q1 CONTRASTS
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
    paste0("/n/home_fasse/txia/UKB_CKM/Data/t2d_stage_specific_", yr, "y_risk_contrast_Q5_vs_Q1.csv"),
    row.names = FALSE
  )
}

write.csv(
  all_risk_contrasts,
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_stage_specific_risk_contrast_all_horizons.csv",
  row.names = FALSE
)

############################################################
## 13. THRESHOLD HEATMAPS
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
  "/n/home_fasse/txia/UKB_CKM/Data/t2d_threshold_categories_5y_10y.csv",
  row.names = FALSE
)

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

  p_thr <- ggplot(temp, aes(x = factor(CKM_Stage), y = Quintile, fill = Risk_Category)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = tile_label), size = 3.8, color = "black", lineheight = 0.92) +
    labs(
      title = paste0("Clinically actionable ", yr, "-year T2D risk categories"),
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
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_t2d_threshold_heatmap_", yr, "y_fullscore.png"),
    p_thr, width = 18, height = 8, dpi = 300, bg = "white"
  )
}

############################################################
## 14. RISK GRADIENT FIGURE
############################################################

gradient_year <- if (5 %in% years_to_run) 5 else min(years_to_run)

risk_data <- read.csv(
  paste0("/n/home_fasse/txia/UKB_CKM/Data/t2d_observed_", gradient_year, "y_risk_by_stage_and_quintile.csv")
)

risk_data$CKM_Stage <- factor(risk_data$CKM_Stage)
risk_data$Quintile <- factor(risk_data$Quintile, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

p_gradient <- ggplot(
  risk_data,
  aes(x = CKM_Stage, y = Observed_Risk, color = Quintile, group = Quintile)
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.2) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(
    title = paste0("Observed ", gradient_year, "-year T2D risk across CKM stage and proteomics quintile"),
    x = "CKM stage",
    y = paste0("Observed ", gradient_year, "-year T2D risk"),
    color = "Proteomics score quintile"
  ) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

ggsave(
  paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_t2d_risk_gradient_CKM_proteomics_", gradient_year, "y.png"),
  p_gradient, width = 8, height = 5.5, dpi = 300, bg = "white"
)

