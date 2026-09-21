# ============================================================
# single molecule ~ time
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(caret)

load("allid.RData")
load("sparse_med.RData")

sparse_med <- as.matrix(sparse_med)

time_data  <- fread("bloodtime_participant.csv", data.table = FALSE)
time_data  <- time_data[, c(1, 2, 6)]
colnames(time_data) <- c("id", "datetime", "fasting")

chem_data  <- fread("bloodchem_participant.csv",  data.table = FALSE)
cell_data  <- fread("bloodcell_participant.csv",  data.table = FALSE)
covar_data <- fread("democovar.csv",             data.table = FALSE)
lifestyle_data <- fread("lifestylecovar.csv",    data.table = FALSE)
metab_data <- fread("metabolome.csv",            data.table = FALSE)

clean_names <- function(df) {
  colnames(df) <- gsub(" \\| Instance 0", "", colnames(df))
  colnames(df) <- gsub(" ", "_", colnames(df))
  df
}
cell_data  <- clean_names(cell_data)
chem_data  <- clean_names(chem_data)

parsed_times <- ymd_hms(time_data$datetime)

time_data$year    <- year(parsed_times)
time_data$hour    <- hour(parsed_times) +
                     minute(parsed_times) / 60 +
                     second(parsed_times) / 3600

day_of_year       <- yday(parsed_times)
theta             <- day_of_year * 2 * pi / 365
time_data$sinday  <- sin(theta)
time_data$cosday  <- cos(theta)

match_by_id <- function(df, ids) {
  df[match(ids, df[[1]]), , drop = FALSE]
}

covar_aligned     <- match_by_id(covar_data,     time_data$id)
lifestyle_aligned <- match_by_id(lifestyle_data, time_data$id)
sparse_med_aligned <- sparse_med[match(time_data$id, rownames(sparse_med)), , drop = FALSE]

stopifnot(!any(is.na(covar_aligned[[1]])))
stopifnot(!any(is.na(lifestyle_aligned[[1]])))
stopifnot(!any(is.na(rownames(sparse_med_aligned))))

analyze_omics <- function(omics_data, time_data, covar_aligned,
                          lifestyle_aligned, sparse_med_aligned,
                          period = 24) {

  feature_names <- colnames(omics_data)[-1]
  results_list  <- vector("list", length(feature_names))

  sin_time <- sin(2 * pi * time_data$hour / period)
  cos_time <- cos(2 * pi * time_data$hour / period)

  for (i in seq_along(feature_names)) {
    feature <- feature_names[i]

    feature_value <- omics_data[match(time_data$id, omics_data[[1]]), feature]

    dat_demo <- data.frame(
      feature_value = feature_value,
      fasting       = time_data$fasting,
      sinday        = time_data$sinday,
      cosday        = time_data$cosday,
      covar_aligned
    )
    dat_demo <- dat_demo[dat_demo$ethnicity == "British", , drop = FALSE]
    dat_demo$ethnicity <- NULL
    fit_demo <- lm(feature_value ~ ., data = dat_demo)
    r2_demo  <- summary(fit_demo)$r.squared

    dat_life <- data.frame(
      feature_value = feature_value,
      fasting       = time_data$fasting,
      sinday        = time_data$sinday,
      cosday        = time_data$cosday,
      lifestyle_aligned,
      covar_aligned
    )
    dat_life <- dat_life[dat_life$ethnicity == "British", , drop = FALSE]
    dat_life$ethnicity <- NULL
    fit_life <- lm(feature_value ~ ., data = dat_life)
    r2_life  <- summary(fit_life)$r.squared

    dat_med <- data.frame(
      feature_value = feature_value,
      fasting       = time_data$fasting,
      sinday        = time_data$sinday,
      cosday        = time_data$cosday,
      sparse_med_aligned,
      lifestyle_aligned,
      covar_aligned
    )
    dat_med <- dat_med[dat_med$ethnicity == "British", , drop = FALSE]
    dat_med$ethnicity <- NULL
    fit_med <- lm(feature_value ~ ., data = dat_med)
    r2_med  <- summary(fit_med)$r.squared

    dat_full <- data.frame(
      sin_time      = sin_time,
      cos_time      = cos_time,
      feature_value = feature_value,
      fasting       = time_data$fasting,
      sinday        = time_data$sinday,
      cosday        = time_data$cosday,
      sparse_med_aligned,
      lifestyle_aligned,
      covar_aligned
    )
    dat_full <- dat_full[dat_full$ethnicity == "British", , drop = FALSE]
    dat_full$ethnicity <- NULL
    fit_full <- lm(feature_value ~ ., data = dat_full)
    r2_all   <- summary(fit_full)$r.squared

    coef_tab <- summary(fit_full)$coefficients
    ry_res   <- coef_tab["sin_time", ]
    ry_type  <- "sin"
    if (coef_tab["sin_time", 4] > coef_tab["cos_time", 4]) {
      ry_res  <- coef_tab["cos_time", ]
      ry_type <- "cos"
    }

    results_list[[i]] <- data.frame(
      feature      = feature,
      R2_demo      = r2_demo,
      R2_lifestyle = r2_life - r2_demo,
      R2_med       = r2_med - r2_life,
      R2_time      = r2_all - r2_med,
      rhythm_type  = ry_type,
      beta         = ry_res[1],
      se           = ry_res[2],
      t            = ry_res[3],
      p            = ry_res[4]
    )

    gc()
  }

  do.call(rbind, results_list)
}

# ----------分别分析三类组学 ----------
bloodcyc_cell  <- analyze_omics(cell_data,  time_data, covar_aligned,
                                lifestyle_aligned, sparse_med_aligned)
bloodcyc_chem  <- analyze_omics(chem_data,  time_data, covar_aligned,
                                lifestyle_aligned, sparse_med_aligned)
bloodcyc_metab <- analyze_omics(metab_data, time_data, covar_aligned,
                                lifestyle_aligned, sparse_med_aligned)

# ---------- 保存结果 ----------
fwrite(bloodcyc_cell,  "bloodcyc_med.csv")
fwrite(bloodcyc_chem,  "bloodchem_med.csv")
fwrite(bloodcyc_metab, "bloodmetab_med.csv")
