# ============================================================
# sensitivity analysis
#         - nomed：x medication
#         - main：all
#         - time：9–21
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(caret)

# ---------- 1. 读取数据 ----------
load("allid.RData")
load("sparse_med.RData")
load("medsvd.RData")

sparse_med <- as.matrix(sparse_med)

time_data      <- fread("bloodtime_participant.csv", data.table = FALSE)
chem_data      <- fread("bloodchem_participant.csv", data.table = FALSE)
cell_data      <- fread("bloodcell_participant.csv", data.table = FALSE)
covar_data     <- fread("democovar.csv",             data.table = FALSE)
lifestyle_data <- fread("lifestylecovar.csv",        data.table = FALSE)

metab_data <- fread("metabolome.csv", data.table = FALSE)
metab_id   <- fread("~/metabid",      data.table = FALSE)
metab_data <- data.frame(id = metab_id[, 1], metab_data)

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
period        <- 24
ethnicity_col <- "ethnicity"  

# 提前算好时间项
sin_time <- sin(2 * pi * time_data$hour / period)

# ---------- 6. 定义敏感性分析函数 ----------
# 返回：data.frame(cell, model, beta, se, t, p)
run_sensitivity <- function(omics_data, time_data,
                            sparse_med_aligned,
                            lifestyle_aligned, covar_aligned,
                            sin_time,
                            ethnicity_col = "ethic",
                            period = 24,
                            time_window = c(9, 21)) {

  feature_names <- colnames(omics_data)[-1]
  sample_ids    <- time_data$id
  results_list  <- list()
  idx           <- 0

  for (feature in feature_names) {
    feature_value <- omics_data[match(sample_ids, omics_data[[1]]), feature]
    names(feature_value) <- sample_ids

    # ---- 基础 data.frame----
    dat_base <- data.frame(
      sintime = sin_time,
      cell    = feature_value,
      fast    = time_data$fasting,
      sinday  = time_data$sinday,
      cosday  = time_data$cosday,
      lifestyle_aligned,
      covar_aligned
    )
    rownames(dat_base) <- sample_ids

    # ---- 按种族筛选（British），并删除种族列 ----
    dat_british <- dat_base[dat_base[[ethnicity_col]] == "British", , drop = FALSE]
    dat_british[[ethnicity_col]] <- NULL

    # ---- 模型 1：nomed----
    fit_nomed <- lm(cell ~ ., data = dat_british)
    coef_nomed <- summary(fit_nomed)$coefficients["sintime", ]

    idx <- idx + 1
    results_list[[idx]] <- data.frame(
      cell  = feature,
      model = "nomed",
      beta  = coef_nomed[1],
      se    = coef_nomed[2],
      t     = coef_nomed[3],
      p     = coef_nomed[4],
      stringsAsFactors = FALSE
    )

    # ---- 模型 2：main----
    dat_med <- data.frame(
      sintime = sin_time,
      cell    = feature_value,
      fast    = time_data$fasting,
      sinday  = time_data$sinday,
      cosday  = time_data$cosday,
      sparse_med_aligned,
      lifestyle_aligned,
      covar_aligned
    )
    rownames(dat_med) <- sample_ids

    dat_med_british <- dat_med[dat_med[[ethnicity_col]] == "British", , drop = FALSE]
    dat_med_british[[ethnicity_col]] <- NULL

    fit_main <- lm(cell ~ ., data = dat_med_british)
    coef_main <- summary(fit_main)$coefficients["sintime", ]

    idx <- idx + 1
    results_list[[idx]] <- data.frame(
      cell  = feature,
      model = "main",
      beta  = coef_main[1],
      se    = coef_main[2],
      t     = coef_main[3],
      p     = coef_main[4],
      stringsAsFactors = FALSE
    )

    # ---- 模型 3：time（限制采血时间 9–21 点）----
    # 注意：原脚本用全局 time$hour 筛选，行数可能不匹配；
    #       这里用 dat_med_british 自带的 time 信息筛选
    keep_time <- time_data$hour > time_window[1] &
                 time_data$hour < time_window[2]
    names(keep_time) <- sample_ids
    keep_time <- keep_time[rownames(dat_med_british)]

    dat_time <- dat_med_british[keep_time, , drop = FALSE]

    fit_time <- lm(cell ~ ., data = dat_time)
    coef_time <- summary(fit_time)$coefficients["sintime", ]

    idx <- idx + 1
    results_list[[idx]] <- data.frame(
      cell  = feature,
      model = "time",
      beta  = coef_time[1],
      se    = coef_time[2],
      t     = coef_time[3],
      p     = coef_time[4],
      stringsAsFactors = FALSE
    )

    rm(dat_base, dat_british, dat_med, dat_med_british, dat_time,
       fit_nomed, fit_main, fit_time)
    gc()
  }

  do.call(rbind, results_list)
}

# ---------- 7. 分别分析三类组学 ----------
res_cell <- run_sensitivity(
  cell_data, time_data,
  sparse_med_aligned, lifestyle_aligned, covar_aligned,
  sin_time,
  ethnicity_col = ethnicity_col,
  period = period
)

res_chem <- run_sensitivity(
  chem_data, time_data,
  sparse_med_aligned, lifestyle_aligned, covar_aligned,
  sin_time,
  ethnicity_col = ethnicity_col,
  period = period
)

res_metab <- run_sensitivity(
  metab_data, time_data,
  sparse_med_aligned, lifestyle_aligned, covar_aligned,
  sin_time,
  ethnicity_col = ethnicity_col,
  period = period
)


fwrite(res_all[res_all$model == "nomed", ], "sen_nomed.csv")
fwrite(res_all[res_all$model == "main",  ], "sen_main.csv")
fwrite(res_all[res_all$model == "time",  ], "sen_time.csv")

