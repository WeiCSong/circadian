# ============================================================
# PWAS：protein -> shift/disturbance
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(caret)

# ---------- 1. 读取数据 ----------
protein_data   <- fread("ukbprotein.csv",       data.table = FALSE)
time_data      <- fread("bloodtime_participant.csv", data.table = FALSE)
covar_data     <- fread("democovar.csv",         data.table = FALSE)
lifestyle_data <- fread("lifestylecovar.csv",    data.table = FALSE)

time_data <- time_data[, c(1, 2, 6)]
colnames(time_data) <- c("id", "datetime", "fasting")

# 只保留有蛋白数据的样本
time_data <- time_data[time_data$id %in% protein_data[[1]], ]

# ---------- 2. 读取 BART 预测 ----------
# BART 输出列：ID, IID, pred, true, sd, interval
read_bart_shift_disturb <- function(file) {
  d <- fread(file, data.table = FALSE)
  data.frame(
    ID          = d[[1]],
    shift       = d[[4]] - d[[3]],   # true - pred（保留原脚本定义）
    disturbance = d[[6]],             # interval
    stringsAsFactors = FALSE
  )
}

pred_all <- read_bart_shift_disturb("bart/allpred.csv")
pred_g1  <- read_bart_shift_disturb("bart/g1pred.csv")
pred_g2  <- read_bart_shift_disturb("bart/g2pred.csv")

# ---------- 3. 按 id 对齐协变量 ----------
match_by_id <- function(df, ids) {
  df[match(ids, df[[1]]), , drop = FALSE]
}

covar_sub     <- match_by_id(covar_data,     time_data$id)[, -1, drop = FALSE]
lifestyle_sub <- match_by_id(lifestyle_data, time_data$id)[, -1, drop = FALSE]

# 种族列名（请确认）
ethnicity_col <- "ethic"

# 对齐三个预测
pred_all_m <- pred_all[match(time_data$id, pred_all$ID), ]
pred_g1_m  <- pred_g1[match(time_data$id,  pred_g1$ID), ]
pred_g2_m  <- pred_g2[match(time_data$id,  pred_g2$ID), ]

# ---------- 4. 定义 PWAS 函数 ----------
# 对每个蛋白，对指定的 pred_df 的 shift / disturbance 各拟合一次 lm
run_pwas <- function(protein_data, time_data, pred_df,
                     covar_sub, lifestyle_sub, ethnicity_col,
                     group_tag, results_list, idx) {

  sample_ids <- time_data$id

  for (prot in colnames(protein_data)[-1]) {
    protein_value <- protein_data[match(sample_ids, protein_data[[1]]), prot]
    names(protein_value) <- sample_ids

    # ---- 公共协变量部分 ----
    base_covars <- data.frame(
      fasting = time_data$fasting,
      lifestyle_sub,
      covar_sub
    )
    rownames(base_covars) <- sample_ids

    # 按种族筛选 British，并删除种族列
    keep <- base_covars[[ethnicity_col]] == "British"
    base_covars_british <- base_covars[keep, , drop = FALSE]
    base_covars_british[[ethnicity_col]] <- NULL
    keep_ids <- rownames(base_covars_british)

    # ---- 模型 1：shift ----
    dat_shift <- data.frame(
      protein = protein_value[keep_ids],
      shift   = pred_df$shift[match(keep_ids, pred_df$ID)],
      base_covars_british
    )
    fit_shift <- lm(protein ~ ., data = dat_shift)
    coef_shift <- summary(fit_shift)$coefficients["shift", ]

    idx <- idx + 1
    results_list[[idx]] <- data.frame(
      protein = prot,
      group   = group_tag,
      metric  = "shift",
      beta    = coef_shift[1],
      se      = coef_shift[2],
      t       = coef_shift[3],
      p       = coef_shift[4],
      stringsAsFactors = FALSE
    )

    # ---- 模型 2：disturbance ----
    dat_dist <- data.frame(
      protein     = protein_value[keep_ids],
      disturbance = pred_df$disturbance[match(keep_ids, pred_df$ID)],
      base_covars_british
    )
    fit_dist <- lm(protein ~ ., data = dat_dist)
    coef_dist <- summary(fit_dist)$coefficients["disturbance", ]

    idx <- idx + 1
    results_list[[idx]] <- data.frame(
      protein = prot,
      group   = group_tag,
      metric  = "disturbance",
      beta    = coef_dist[1],
      se      = coef_dist[2],
      t       = coef_dist[3],
      p       = coef_dist[4],
      stringsAsFactors = FALSE
    )
  }

  list(results_list = results_list, idx = idx)
}

# ---------- 5. 三次 PWAS ----------
results_list <- list()
idx <- 0

out <- run_pwas(protein_data, time_data, pred_all_m,
                covar_sub, lifestyle_sub, ethnicity_col,
                group_tag = "all", results_list, idx)
results_list <- out$results_list
idx          <- out$idx

out <- run_pwas(protein_data, time_data, pred_g1_m,
                covar_sub, lifestyle_sub, ethnicity_col,
                group_tag = "g1", results_list, idx)
results_list <- out$results_list
idx          <- out$idx

out <- run_pwas(protein_data, time_data, pred_g2_m,
                covar_sub, lifestyle_sub, ethnicity_col,
                group_tag = "g2", results_list, idx)
results_list <- out$results_list
idx          <- out$idx

# ---------- 6. 写出结果 ----------
pwas_table <- do.call(rbind, results_list)
fwrite(pwas_table, "pwas.csv")
