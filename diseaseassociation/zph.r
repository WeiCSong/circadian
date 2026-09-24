# ============================================================
# Cox zph test
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(CMAverse)

# ---------- 1. 读取基础数据 ----------
time_data      <- fread("bloodtime_participant.csv", data.table = FALSE)
covar_data     <- fread("democovar.csv",             data.table = FALSE)
lifestyle_data <- fread("lifestylecovar.csv",        data.table = FALSE)
chro_data      <- fread("common.regenie.pheno",      data.table = FALSE)

icd_files        <- system("ls ICD/*.csv", intern = TRUE)
end_date         <- ymd("2025-01-01")
blood_collection <- ymd_hms(time_data[, 2])

chro_data <- chro_data[match(time_data[, 1], chro_data[, 1]), c(1, 7)]

# ---------- 2. 辅助函数 ----------
rank_qnorm <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# 从 BART 输出读取（列：ID, IID, pred, true, sd, interval）
read_pred_bart <- function(file) {
  p <- fread(file, data.table = FALSE)
  data.frame(
    ID          = p[[1]],
    shift       = rank_qnorm(p[[3]] - p[[4]]),   # pred - true
    disturbance = rank_qnorm(p[[6]]),             # interval
    stringsAsFactors = FALSE
  )
}

# 安全地从 cox.zph 结果中取 p 值
safe_zph_p <- function(zph_table, var) {
  if (var %in% rownames(zph_table)) zph_table[var, "p"] else NA_real_
}

# ---------- 3. 主分析函数 ----------
analyze_zph <- function(pred, chro_file, output_file,
                        time_data, covar_data, lifestyle_data, chro_data,
                        icd_files, blood_collection, end_date,
                        bonferroni_n = 3618) {

  # ---- 3.1 从 *_chro.csv 读出显著疾病 ----
  res <- fread(chro_file, data.table = FALSE)
  # 期望列：dis, metric, coef, exp.coef, se, z, p
  # 按列名取，避免位置索引
  sig_dis <- unique(res$dis[res$metric == "disturbance" & res$p < 0.05 / bonferroni_n])

  if (length(sig_dis) == 0) {
    message("No significant diseases in ", chro_file)
    return(invisible(NULL))
  }

  message("Analyzing ", length(sig_dis), " significant diseases for ", output_file)

  # ---- 3.2 对齐预测与协变量 ----
  pred_matched <- pred[match(time_data[, 1], pred$ID), ]
  covar_sub     <- covar_data[match(time_data[, 1], covar_data[, 1]), 2:16]
  lifestyle_sub <- lifestyle_data[match(time_data[, 1], lifestyle_data[, 1]), -1]
  chro_sub      <- chro_data[match(time_data[, 1], chro_data[, 1]), 2]

  results_list <- list()
  idx <- 0

  for (path in icd_files) {
    diag <- fread(path, data.table = FALSE)

    for (dis in colnames(diag)[-1]) {
      if (!dis %in% sig_dis) next

      diagnosis_date <- ymd(diag[match(time_data[, 1], diag[, 1]), dis])

      status <- ifelse(!is.na(diagnosis_date) & diagnosis_date > blood_collection, 1, 0)
      dtime  <- ifelse(status == 1,
                       as.numeric(difftime(diagnosis_date, blood_collection, units = "days")),
                       as.numeric(difftime(end_date, blood_collection, units = "days")))

      former_diag <- which(!is.na(diagnosis_date) & diagnosis_date <= blood_collection)
      dtime[former_diag]  <- NA
      status[former_diag] <- NA

      df <- data.frame(
        status      = status,
        time        = dtime,
        covar_sub,
        lifestyle_sub,
        shift       = pred_matched$shift,
        disturbance = pred_matched$disturbance,
        chronotype  = chro_sub
      )
      df <- df[!is.na(df$status) & !is.na(df$time) &
               !is.na(df$shift)  & !is.na(df$chronotype), ]
      df[is.na(df)] <- 0

      cox_model <- coxph(Surv(time, status) ~ ., data = df)
      zph       <- cox.zph(cox_model)
      zph_table <- zph$table

      idx <- idx + 1
      results_list[[idx]] <- data.frame(
        dis          = dis,
        zph_global_p = zph_table["GLOBAL", "p"],
        zph_shift    = safe_zph_p(zph_table, "shift"),
        zph_disturb  = safe_zph_p(zph_table, "disturbance"),
        zph_chrono   = safe_zph_p(zph_table, "chronotype"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results_list) > 0) {
    fwrite(do.call(rbind, results_list), output_file)
  }
}

# ---------- 4. 三次分析 ----------
# 4.1 all
pred_all <- read_pred_bart("bart/allpred.csv")
analyze_zph(pred_all, "allbart2dis_chro.csv", "allbart2dis_zph.csv",
            time_data, covar_data, lifestyle_data, chro_data,
            icd_files, blood_collection, end_date)

# 4.2 g1
pred_g1 <- read_pred_bart("bart/g1pred.csv")
analyze_zph(pred_g1, "g1bart2dis_chro.csv", "g1bart2dis_zph.csv",
            time_data, covar_data, lifestyle_data, chro_data,
            icd_files, blood_collection, end_date)

# 4.3 g2
pred_g2 <- read_pred_bart("bart/g2pred.csv")
analyze_zph(pred_g2, "g2bart2dis_chro.csv", "g2bart2dis_zph.csv",
            time_data, covar_data, lifestyle_data, chro_data,
            icd_files, blood_collection, end_date)
