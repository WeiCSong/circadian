# ============================================================
# gam to wide data frame
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

lifestyle_aligned <- match_by_id(lifestyle_data, time_data$id)
covar_aligned     <- match_by_id(covar_data,     time_data$id)
sparse_med_aligned <- sparse_med[match(time_data$id, rownames(sparse_med)), , drop = FALSE]

stopifnot(!any(is.na(covar_aligned[[1]])))
stopifnot(!any(is.na(lifestyle_aligned[[1]])))
stopifnot(!any(is.na(rownames(sparse_med_aligned))))

# ---------- 5. 预测网格 ----------
pred_times <- seq(from = 9, to = 21, length.out = 200)

# ---------- 6. 定义 GAM 拟合 + 预测函数 ----------
# 返回：矩阵，行 = 预测时间点，列 = 特征名
fit_gam_curves <- function(omics_data, time_data,
                           sparse_med_aligned,
                           lifestyle_aligned, covar_aligned,
                           pred_times,
                           ethnicity_col = "ethnicity",
                           k_spline = 25, n_threads = 19) {

  feature_names <- colnames(omics_data)[-1]
  sample_ids    <- time_data$id

  # ---- 预先构造协变量部分（不含特征值），只做一次 ----
  base_covars <- data.frame(
    fasting = time_data$fasting,
    sinday  = time_data$sinday,
    cosday  = time_data$cosday,
    sparse_med_aligned,
    lifestyle_aligned,
    covar_aligned
  )
  rownames(base_covars) <- sample_ids

  # 按种族筛选（一次完成）
  keep <- base_covars[[ethnicity_col]] == "British"
  base_covars <- base_covars[keep, , drop = FALSE]
  base_covars[[ethnicity_col]] <- NULL
  keep_ids <- rownames(base_covars)

  # 识别因子列及其水平（用于 newdata 中保留水平，避免 contrasts 报错）
  factor_levels <- lapply(base_covars, function(x) {
    if (is.character(x) || is.factor(x)) unique(as.character(x)) else NULL
  })
  factor_levels <- factor_levels[!vapply(factor_levels, is.null, logical(1))]

  # 构造“典型个体”协变量（数值取均值，因子取众数）
  typical_values <- base_covars %>%
    summarise(
      across(where(is.numeric),   ~ mean(., na.rm = TRUE)),
      across(where(is.character), ~ {
        tab <- table(.)
        names(tab)[which.max(tab)]
      }),
      across(where(is.factor),    ~ {
        tab <- table(.)
        factor(names(tab)[which.max(tab)], levels = levels(.))
      })
    )

  # 构造预测用 newdata（只做一次，所有特征共用）
  newdata <- typical_values[rep(1, length(pred_times)), , drop = FALSE]
  newdata$time <- pred_times
  # 保留因子水平
  for (col in names(factor_levels)) {
    newdata[[col]] <- factor(newdata[[col]], levels = factor_levels[[col]])
  }
  newdata <- newdata[, c("time", setdiff(names(newdata), "time")), drop = FALSE]

  # 结果矩阵：行 = 预测时间点，列 = 特征
  curve_mat <- matrix(NA_real_, nrow = length(pred_times), ncol = length(feature_names),
                      dimnames = list(NULL, feature_names))

  cov_names <- setdiff(names(base_covars), character(0))
  rhs_terms <- paste(c("s(time, bs = 'cc', k =", paste0(k_spline, ")")),
                     paste(cov_names, collapse = " + "))

  for (i in seq_along(feature_names)) {
    feature <- feature_names[i]

    feature_value <- omics_data[match(sample_ids, omics_data[[1]]), feature]
    names(feature_value) <- sample_ids
    feature_value <- feature_value[keep_ids]

    dat <- data.frame(time = time_data$hour[keep],
                      metric = feature_value,
                      base_covars)
    rownames(dat) <- keep_ids

    formula_str <- paste("metric ~", rhs_terms)
    gam_model <- bam(as.formula(formula_str),
                     data = dat,
                     family = gaussian(),
                     method = "fREML",
                     discrete = TRUE,
                     nthreads = n_threads)

    pred <- predict(gam_model, newdata = newdata,
                    type = "response", se.fit = TRUE)
    curve_mat[, i] <- pred$fit

    rm(dat, gam_model, pred)
    gc()
    message("Done: ", feature)
  }

  curve_mat
}

# ---------- 7. 分别拟合三类组学 ----------
curve_cell  <- fit_gam_curves(cell_data,  time_data,
                              sparse_med_aligned,
                              lifestyle_aligned, covar_aligned,
                              pred_times)

curve_chem  <- fit_gam_curves(chem_data,  time_data,
                              sparse_med_aligned,
                              lifestyle_aligned, covar_aligned,
                              pred_times)

curve_metab <- fit_gam_curves(metab_data, time_data,
                              sparse_med_aligned,
                              lifestyle_aligned, covar_aligned,
                              pred_times)

# ---------- 8. 写出结果 ----------
write_curves <- function(curve_mat, pred_times, file) {
  out <- data.frame(time = pred_times, curve_mat, check.names = FALSE)
  fwrite(out, file)
}

write_curves(curve_cell,  pred_times, "gam_cell.csv")
write_curves(curve_chem,  pred_times, "gam_chem.csv")
write_curves(curve_metab, pred_times, "gam_metab.csv")
