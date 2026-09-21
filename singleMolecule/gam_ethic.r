# ============================================================
# 种族差异 GAM 分析脚本
# 目的：对 cell / chem / metab 的每个特征，比较 British 与
#       AFR / Caribbean / EAS / NBW / SAS 的昼夜曲线差异。
#       用 s(time, by = grp_ind) 的 p-value 和 edf 衡量差异。
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(biglm)
library(lightgbm)
library(stochtree)

# ---------- 1. 读取数据 ----------
load("sparse_med.RData")
sparse_med <- as.matrix(sparse_med)

time_data      <- fread("bloodtime_participant.csv", data.table = FALSE)
chem_data      <- fread("bloodchem_participant.csv", data.table = FALSE)
cell_data      <- fread("bloodcell_participant.csv", data.table = FALSE)
covar_data     <- fread("democovar.csv",             data.table = FALSE)
lifestyle_data <- fread("lifestylecovar.csv",        data.table = FALSE)
metab_data     <- fread("../ukbmetab.csv",           data.table = FALSE)

# ---------- 2. 清洗列名 ----------
clean_names <- function(df) {
  colnames(df) <- gsub(" \\| Instance 0", "", colnames(df))
  colnames(df) <- gsub(" ", "_", colnames(df))
  df
}
cell_data <- clean_names(cell_data)
chem_data <- clean_names(chem_data)

# ---------- 3. 时间变量派生 ----------
time_data <- time_data[, c(1, 2, 6)]
colnames(time_data) <- c("id", "datetime", "fasting")

parsed_times <- ymd_hms(time_data$datetime)

time_data$year    <- year(parsed_times)
time_data$hour    <- hour(parsed_times) +
                     minute(parsed_times) / 60 +
                     second(parsed_times) / 3600
time_data$day     <- yday(parsed_times)

theta             <- time_data$day * 2 * pi / 365
time_data$sinday  <- sin(theta)
time_data$cosday  <- cos(theta)

# ---------- 4. 按 id 对齐协变量 ----------
match_by_id <- function(df, ids) {
  df[match(ids, df[[1]]), , drop = FALSE]
}

lifestyle_aligned  <- match_by_id(lifestyle_data, time_data$id)
covar_aligned      <- match_by_id(covar_data,     time_data$id)
sparse_med_aligned <- sparse_med[match(time_data$id, rownames(sparse_med)), , drop = FALSE]

stopifnot(!any(is.na(covar_aligned[[1]])))
stopifnot(!any(is.na(lifestyle_aligned[[1]])))
stopifnot(!any(is.na(rownames(sparse_med_aligned))))

# ---------- 5. 公共参数 ----------
comparison_groups <- c("AFR", "Caribbean", "EAS", "NBW", "SAS")
ethnicity_col     <- "ethnicity"  
k_spline          <- 25
n_threads         <- 15

# ---------- 6. 定义种族差异分析函数 ----------
# 返回：data.frame(group, metric, p, edf)
analyze_ethnicity_diff <- function(omics_data, time_data,
                                   sparse_med_aligned,
                                   lifestyle_aligned, covar_aligned,
                                   comparison_groups,
                                   ethnicity_col = "ethic",
                                   include_season = TRUE,
                                   k_spline = 25, n_threads = 15) {

  feature_names <- colnames(omics_data)[-1]
  sample_ids    <- time_data$id
  results_list  <- list()
  idx           <- 0

  # ---- 预先构造协变量部分（不含特征值），只做一次 ----
  base_covars <- data.frame(
    fasting = time_data$fasting,
    sparse_med_aligned,
    lifestyle_aligned,
    covar_aligned
  )
  if (include_season) {
    base_covars$sinday <- time_data$sinday
    base_covars$cosday <- time_data$cosday
  }
  rownames(base_covars) <- sample_ids

  # 协变量名（排除种族列，因为要单独用 grp_ind 处理）
  cov_names <- setdiff(names(base_covars), ethnicity_col)

  for (feature in feature_names) {
    feature_value <- omics_data[match(sample_ids, omics_data[[1]]), feature]
    names(feature_value) <- sample_ids

    dat_all <- data.frame(
      time   = time_data$hour,
      metric = feature_value,
      base_covars
    )
    rownames(dat_all) <- sample_ids

    for (grp in comparison_groups) {
      dat <- dat_all[dat_all[[ethnicity_col]] %in% c("British", grp), , drop = FALSE]
      dat$grp_ind <- as.numeric(dat[[ethnicity_col]] == grp)

      formula_str <- paste(
        "metric ~ grp_ind +",
        "s(time, bs = 'cc', k =", k_spline, ") +",
        "s(time, by = grp_ind, bs = 'cc', k =", k_spline, ") +",
        paste(cov_names, collapse = " + ")
      )

      gam_model <- bam(as.formula(formula_str),
                       data = dat,
                       family = gaussian(),
                       method = "fREML",
                       discrete = TRUE,
                       nthreads = n_threads)

      s_table  <- summary(gam_model)$s.table
      diff_row <- grep("grp_ind", rownames(s_table))[1]

      idx <- idx + 1
      results_list[[idx]] <- data.frame(
        group  = grp,
        metric = feature,
        p      = s_table[diff_row, "p-value"],
        edf    = s_table[diff_row, "edf"],
        stringsAsFactors = FALSE
      )

      rm(dat, gam_model, s_table)
      gc()
    }

    rm(dat_all)
    gc()
    message("Done: ", feature)
  }

  do.call(rbind, results_list)
}

# ---------- 7. 分别分析三类组学 ----------
# 7.1 cell
res_cell <- analyze_ethnicity_diff(
  cell_data, time_data,
  sparse_med_aligned, lifestyle_aligned, covar_aligned,
  comparison_groups,
  ethnicity_col = ethnicity_col,
  include_season = TRUE,
  k_spline = k_spline, n_threads = n_threads
)
fwrite(res_cell, "cell_ethic_diff.csv")

# 7.2 chem
res_chem <- analyze_ethnicity_diff(
  chem_data, time_data,
  sparse_med_aligned, lifestyle_aligned, covar_aligned,
  comparison_groups,
  ethnicity_col = ethnicity_col,
  include_season = TRUE,
  k_spline = k_spline, n_threads = n_threads
)
fwrite(res_chem, "chem_ethic_diff.csv")

# 7.3 metab
res_metab <- analyze_ethnicity_diff(
  metab_data, time_data,
  sparse_med_aligned, lifestyle_aligned, covar_aligned,
  comparison_groups,
  ethnicity_col = ethnicity_col,
  include_season = TRUE,
  k_spline = k_spline, n_threads = n_threads
)
fwrite(res_metab, "metab_ethic_diff.csv")

