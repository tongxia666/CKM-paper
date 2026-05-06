############################################################
## CKD HEATMAP-ONLY CODE
## - Incident CKD definition
## - Observed risk by CKM stage and proteomics quintile
## - All yearly horizons
## - Overall heatmap
## - CKM on x-axis, quintile on y-axis
## - Tile labels: risk %, 95% CI, Events/N
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

ckd_association <- fread("/n/home_fasse/txia/UKB_CKM/Data/ukb_filtered_with_polr_manual_protein_score.csv")
ukb_data        <- fread("/n/home_fasse/txia/UKB_CKM/Data/ckdvariable.csv")

############################################################
## 2. BASIC PREP
############################################################

ckd_association <- ckd_association[!is.na(ckm), ]

ckd_association$sbp_use <- ifelse(
  !is.na(ckd_association$f.4080.0.0.x),
  ckd_association$f.4080.0.0.x,
  ckd_association$f.93.0.0.x
)

ckd_association$protein_score <- ckd_association$manual_protein_score
ckd_association[, protein_score := qnorm(rank(protein_score, na.last = "keep") / (.N + 1))]

ckd_association <- ckd_association %>%
  left_join(ukb_data, by = "f.eid")

############################################################
## 3. DEFINE INCIDENT CKD
############################################################

ckd_outcomes <- c(
  "f.132032.0.0", "f.42026.0.0"
)

for (outcome in ckd_outcomes) {
  ckd_association[[outcome]] <- as.numeric(as.Date(ckd_association[[outcome]]))
}

ckd_association$dt_visit0 <- as.numeric(as.Date(ckd_association$f.21842.0.0))

ckd_association$ckd_dxdt <- apply(
  as.matrix(ckd_association[, ..ckd_outcomes]),
  1,
  function(x) {
    min_val <- min(x, na.rm = TRUE)
    if (is.infinite(min_val)) NA else min_val
  }
)

ckd_association$ckd_dxdt[is.na(ckd_association$ckd_dxdt)] <- Inf

ckd_association$ckd_bsln <- ifelse(
  ckd_association$ckd_dxdt < ckd_association$dt_visit0,
  1, 0
)

incident_data <- ckd_association %>%
  filter(ckd_bsln == 0)

incident_data <- incident_data[ckm != 4]

incident_data <- incident_data %>%
  mutate(
    incident_ckd_time = ckd_dxdt,
    event = ifelse(!is.na(incident_ckd_time) & incident_ckd_time != Inf, 1, 0),
    censor_time = pmin(
      as.numeric(as.Date(f.40000.0.0)),
      as.numeric(as.Date(f.191.0.0)),
      as.numeric(as.Date("2022-11-30")),
      na.rm = TRUE
    ),
    time_to_event = ifelse(event == 1, incident_ckd_time - dt_visit0, censor_time - dt_visit0)
  )

cat("Baseline and incident CKD have been defined successfully.\n")
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
  "/n/home_fasse/txia/UKB_CKM/Data/ckd_followup_summary_years.csv",
  row.names = FALSE
)

############################################################
## 6. DEFINE ALL YEARLY HORIZONS
############################################################

max_followup_year <- floor(max(dat$time_to_event, na.rm = TRUE) / 365.25)
years_to_run <- seq(1, max_followup_year, by = 1)

cat("Running CKD observed-risk heatmaps for yearly horizons:\n")
print(years_to_run)

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

  df$tile_label <- paste0(
    sprintf("%.1f%%", 100 * df$Observed_Risk),
    "\n",
    df$Events, "/", df$N
  )

  ggplot(df, aes(x = CKM_Stage, y = Quintile, fill = Observed_Risk)) +
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
    paste0("/n/home_fasse/txia/UKB_CKM/Data/ckd_observed_", yr, "y_risk_by_stage_and_quintile.csv"),
    row.names = FALSE
  )
}

write.csv(
  all_risk_results,
  "/n/home_fasse/txia/UKB_CKM/Data/ckd_observed_risk_all_horizons_by_stage_and_quintile.csv",
  row.names = FALSE
)

############################################################
## 9. YEARLY HEATMAPS + YEARLY 95% CI TABLES
############################################################

for (yr in years_to_run) {
  temp <- all_risk_results %>% filter(Horizon_Years == yr)
  if (nrow(temp) == 0) next

  temp$CKM_Stage <- factor(temp$CKM_Stage, levels = sort(unique(temp$CKM_Stage)))
  temp$Quintile  <- factor(temp$Quintile, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

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
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Supplemental_Table_ckd_heatmap_", yr, "y_95CI.csv"),
    row.names = FALSE
  )

  p_heatmap <- make_heatmap_plot(
    df = temp,
    fill_title = paste0(yr, "-year risk"),
    plot_title = paste0("Observed ", yr, "-year CKD risk by CKM stage and proteomics quintile")
  )

  ggsave(
    paste0("/n/home_fasse/txia/UKB_CKM/Data/Figure_heatmap_ckd_", yr, "y_risk_stage_by_quintile_fullscore.png"),
    p_heatmap,
    width = 15,
    height = 10,
    dpi = 300,
    bg = "white"
  )
}

############################################################
## 10. OVERALL HEATMAP + SUPPLEMENTAL 95% CI TABLE
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
  "/n/home_fasse/txia/UKB_CKM/Data/ckd_observed_overall_risk_by_stage_and_quintile.csv",
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
  "/n/home_fasse/txia/UKB_CKM/Data/Supplemental_Table_ckd_overall_heatmap_95CI.csv",
  row.names = FALSE
)

p_heatmap_overall <- make_heatmap_plot(
  df = overall_risk_results,
  fill_title = "Overall risk",
  plot_title = paste0(
    "Overall cumulative CKD risk by CKM stage and proteomics quintile\n",
    "(through ", max_followup_year, "-year follow-up)"
  )
)

ggsave(
  "/n/home_fasse/txia/UKB_CKM/Data/Figure_heatmap_overall_ckd_risk_stage_by_quintile_fullscore.png",
  p_heatmap_overall,
  width = 15,
  height = 10,
  dpi = 300,
  bg = "white"
)

