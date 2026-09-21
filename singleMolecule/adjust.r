# ============================================================
# 血液/代谢组学残差提取脚本（完善版）
# 目的：对每个特征，用协变量（生活方式 + 人口学）做线性回归，
#       取残差并输出为按样本 id 排列的矩阵（每个特征一行）。
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(biglm)

# ---------- 1. 读取数据 ----------
load("sparse_med.RData")
sparse_med <- as.matrix(sparse_med)

time_data      <- fread("bloodtime_participant.csv", data.table = FALSE)
chem_data      <- fread("bloodchem_participant.csv", data.table = FALSE)
cell_data      <- fread("bloodcell_participant.csv", data.table = FALSE)
covar_data     <- fread("democovar.csv",             data.table = FALSE)
lifestyle_data <- fread("lifestylecovar.csv",        data.table = FALSE)

metab_data <- fread("../metabolome.csv", data.table = FALSE)
metab_id   <- fread("~/metabid",         data.table = FALSE)
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

theta             <- time_data$hour * 2 * pi / 24
time_data$sinhour <- sin(theta)
time_data$coshour <- cos(theta)

# ---------- 4. 缺失值填补 ----------
# 原脚本 f() 用于把 NA 替换为列中位数
impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

time_data$fasting <- impute_median(time_data$fasting)

# bcovar / covar 是 data.frame，apply 后返回矩阵，需还原为 data.frame
lifestyle_data <- as.data.frame(apply(lifestyle_data, 2, impute_median))
covar_data     <- as.data.frame(covar_data)
covar_cols     <- c(2, 4:14, 16)
covar_data[, covar_cols] <- apply(covar_data[, covar_cols], 2, impute_median)

# ---------- 5. 按 id 对齐样本 ----------
match_by_id <- function(df, ids) {
  df[match(ids, df[[1]]), , drop = FALSE]
}

lifestyle_aligned <- match_by_id(lifestyle_data, time_data$id)
covar_aligned     <- match_by_id(covar_data,     time_data$id)
med_aligned     <- match_by_id(sparse_med,     time_data$id)
stopifnot(!any(is.na(covar_aligned[[1]])))
stopifnot(!any(is.na(lifestyle_aligned[[1]])))

# ---------- 6. 定义残差提取函数 ----------
# 返回：矩阵，行 = 特征名，列 = 样本 id
extract_residuals <- function(omics_data, time_data,
                              lifestyle_aligned, covar_aligned,
                              ethnicity_col = "ethnicity") {

  feature_names <- colnames(omics_data)[-1]
  sample_ids    <- time_data$id
  n_samples     <- length(sample_ids)

  # 预先构造协变量部分（不含特征值），只做一次匹配和筛选
  base_covars <- data.frame(
    fasting = time_data$fasting,
    sinhour = time_data$sinhour,
    coshour = time_data$coshour,
    lifestyle_aligned,
    covar_aligned,
    med_aligned
  )
  rownames(base_covars) <- sample_ids

  # 按种族筛选（一次完成，后续特征共用）
  keep <- base_covars[[ethnicity_col]] == "British"
  base_covars <- base_covars[keep, , drop = FALSE]
  base_covars[[ethnicity_col]] <- NULL
  keep_ids <- rownames(base_covars)
  n_keep   <- length(keep_ids)

  # 结果矩阵：每个特征一行
  residual_mat <- matrix(NA_real_, nrow = length(feature_names), ncol = n_keep,
                         dimnames = list(feature_names, keep_ids))

  for (i in seq_along(feature_names)) {
    feature <- feature_names[i]

    feature_value <- omics_data[match(sample_ids, omics_data[[1]]), feature]
    names(feature_value) <- sample_ids
    feature_value <- feature_value[keep_ids]

    dat <- data.frame(feature_value = feature_value, base_covars)
    fit <- lm(feature_value ~ ., data = dat)

    residual_mat[i, ] <- round(resid(fit), 3)

    rm(dat, fit)
    gc()
  }

  residual_mat
}

# ---------- 7. 分别提取三类组学的残差 ----------
resid_cell  <- extract_residuals(cell_data,  time_data,
                                 lifestyle_aligned, covar_aligned)
resid_chem  <- extract_residuals(chem_data,  time_data,
                                 lifestyle_aligned, covar_aligned)
resid_metab <- extract_residuals(metab_data, time_data,
                                 lifestyle_aligned, covar_aligned)

# ---------- 8. 写出结果 ----------
# 每个文件：第一行为样本 id（列名），第一列为特征名
write_residuals <- function(resid_mat, file) {
  out <- data.frame(feature = rownames(resid_mat), resid_mat,
                    check.names = FALSE)
  fwrite(out, file, sep = "\t", quote = FALSE, na = "NA")
}

write_residuals(resid_cell,  "adjcell.txt.gz")
write_residuals(resid_chem,  "adjchem.txt.gz")
write_residuals(resid_metab, "adjmetab.txt.gz")
