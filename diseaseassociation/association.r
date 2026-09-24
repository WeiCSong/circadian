# ============================================================
# 与 ICD 疾病的横断面关联（lm）和纵向关联（Cox）。
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(survival)
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

# ---------- 2. 辅助函数 ----------
# rank-based inverse normal transformation
rank_qnorm <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# 从 xgboost 输出读取（已含 shift / disturbance 列）
read_pred_xgb <- function(file) {
  p <- fread(file, data.table = FALSE)
  data.frame(
    ID          = p[[1]],
    shift       = rank_qnorm(p[["shift"]]),
    disturbance = rank_qnorm(p[["disturbance"]]),
    stringsAsFactors = FALSE
  )
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

# ---------- 3. 主分析函数 ----------
# 对每个 ICD 文件、每个疾病：
#   (1) lm:  shift / disturbance ~ disease + 协变量
#   (2) cox: Surv(time, status) ~ shift + disturbance + chronotype + 协变量
analyze_dis <- function(pred, output_prefix,
                        time_data, covar_data, lifestyle_data, chro_data,
                        icd_files, blood_collection, end_date) {


  pred_matched <- pred[match(time_data[, 1], pred$ID), ]


  covar_sub     <- covar_data[match(time_data[, 1], covar_data[, 1]), 2:16]
  lifestyle_sub <- lifestyle_data[match(time_data[, 1], lifestyle_data[, 1]), -1]
  chro_sub      <- chro_data[match(time_data[, 1], chro_data[, 1]), 2]

  lm_results  <- list()
  cox_results <- list()
  idx_lm  <- 0
  idx_cox <- 0

  for (path in icd_files) {
    diag <- fread(path, data.table = FALSE)

    for (dis in colnames(diag)[-1]) {
      diagnosis_date <- ymd(diag[match(time_data[, 1], diag[, 1]), dis])

      status <- ifelse(!is.na(diagnosis_date) & diagnosis_date > blood_collection, 1, 0)
      dtime  <- ifelse(status == 1,
                       as.numeric(difftime(diagnosis_date, blood_collection, units = "days")),
                       as.numeric(difftime(end_date, blood_collection, units = "days")))

      former_diag <- which(!is.na(diagnosis_date) & diagnosis_date <= blood_collection)
      dtime[former_diag]  <- NA
      status[former_diag] <- NA

      # ===== (1) 横断面：lm =====
      df_lm <- data.frame(
        diag         = 0,
        covar_sub,
        lifestyle_sub,
        shift        = pred_matched$shift,
        disturbance  = pred_matched$disturbance,
        chronotype   = chro_sub
      )
      df_lm[former_diag, "diag"] <- 1
      df_lm <- df_lm[!is.na(df_lm$diag) & !is.na(df_lm$shift), ]
      df_lm[is.na(df_lm)] <- 0

      # 只在病例数 > 10 时拟合（避免过拟合 / 无解）
      if (sum(df_lm$diag == 1) > 10) {
        # 从预测变量中排除 shift / disturbance / chronotype
        pred_vars <- df_lm %>% select(-shift, -disturbance, -chronotype)

        fit_shift <- lm(shift ~ ., data = data.frame(shift = df_lm$shift, pred_vars))
        fit_dist  <- lm(disturbance ~ .,
                        data = data.frame(disturbance = df_lm$disturbance, pred_vars))

        coef_shift <- summary(fit_shift)$coefficients["diag", ]
        coef_dist  <- summary(fit_dist)$coefficients["diag", ]

        idx_lm <- idx_lm + 1
        lm_results[[idx_lm]] <- data.frame(
          dis    = dis,
          metric = c("shift", "disturbance"),
          rbind(coef_shift, coef_dist),
          stringsAsFactors = FALSE
        )
      }

      # ===== (2) 纵向：Cox 模型 =====
      df_cox <- data.frame(
        status       = status,
        time         = dtime,
        covar_sub,
        lifestyle_sub,
        shift        = pred_matched$shift,
        disturbance  = pred_matched$disturbance,
        chronotype   = chro_sub
      )
      df_cox <- df_cox[!is.na(df_cox$status) & !is.na(df_cox$time) &
                       !is.na(df_cox$shift)  & !is.na(df_cox$chronotype), ]
      df_cox[is.na(df_cox)] <- 0

      cox_model <- coxph(Surv(time, status) ~ ., data = df_cox)
      coef_all  <- summary(cox_model)$coefficients
      coef_sel  <- coef_all[c("shift", "disturbance", "chronotype"), ]

      idx_cox <- idx_cox + 1
      cox_results[[idx_cox]] <- data.frame(
        dis    = dis,
        metric = c("shift", "disturbance", "chronotype"),
        coef_sel,
        stringsAsFactors = FALSE
      )
    }
  }

  # 一次性写出，带表头
  if (length(lm_results) > 0) {
    fwrite(do.call(rbind, lm_results),
           paste0(output_prefix, "_lm.csv"))
  }
  if (length(cox_results) > 0) {
    fwrite(do.call(rbind, cox_results),
           paste0(output_prefix, "_cox.csv"))
  }
}

# ---------- 4. 三次分析 ----------
# 4.1 bart all
pred_all <- read_pred_bart("bart/allpred.csv")
analyze_dis(pred_xgb, "alllin2dis_chro",
            time_data, covar_data, lifestyle_data, chro_data,
            icd_files, blood_collection, end_date)

# 4.2 bart g1
pred_g1 <- read_pred_bart("bart/g1pred.csv")
analyze_dis(pred_g1, "g1bart2dis_chro",
            time_data, covar_data, lifestyle_data, chro_data,
            icd_files, blood_collection, end_date)

# 4.3 bart g2
pred_g2 <- read_pred_bart("bart/g2pred.csv")
analyze_dis(pred_g2, "g2bart2dis_chro",
            time_data, covar_data, lifestyle_data, chro_data,
            icd_files, blood_collection, end_date)
