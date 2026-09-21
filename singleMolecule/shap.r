# ============================================================
# shift/disturbance -> xgboost -> shap
# ============================================================

library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(caret)
library(lightgbm)
library(shapviz)
library(ggplot2)

# ---------- 1. 读取数据 ----------
load("allid.RData")
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

# ---------- 4. 构造特征矩阵 X（样本 × 特征）----------
match_by_id <- function(df, ids) {
  df[match(ids, df[[1]]), , drop = FALSE]
}

cell_aligned  <- match_by_id(cell_data,  time_data$id)
chem_aligned  <- match_by_id(chem_data,  time_data$id)
metab_aligned <- match_by_id(metab_data, time_data$id)

# 去掉第一列 id，只保留特征
X_all <- cbind(
  cell_aligned[,  -1, drop = FALSE],
  chem_aligned[,  -1, drop = FALSE],
  metab_aligned[, -1, drop = FALSE]
)
rownames(X_all) <- time_data$id
X_all <- as.matrix(X_all)

# 列名合法化
colnames(X_all) <- gsub("[()\\-]", ".", colnames(X_all))

# ---------- 5. 读取聚类标签 ----------
load("cluster_labels.RData")
g1_features <- names(cluster_labels)[cluster_labels == 1]
g2_features <- names(cluster_labels)[cluster_labels == 2]

# 对齐特征名（若 cluster_labels 的 names 未做合法化，需同步替换）
g1_features <- gsub("[()\\-]", ".", g1_features)
g2_features <- gsub("[()\\-]", ".", g2_features)

# 三组特征矩阵
X_g1 <- X_all[, intersect(g1_features, colnames(X_all)), drop = FALSE]
X_g2 <- X_all[, intersect(g2_features, colnames(X_all)), drop = FALSE]

# ---------- 6. 定义辅助函数 ----------
# rank-qnorm 归一化
rank_qnorm <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# 从 BART 预测文件读取 shift 和 disturbance
read_bart_pred <- function(file) {
  pred <- fread(file, data.table = FALSE)
  # 按列名取，避免位置索引出错
  # 期望列：ID, pred, true, sd, interval
  shift       <- pred$pred - pred$true          # 预测偏差
  disturbance <- pred$interval                   # 后验区间宽度（不确定性）

  data.frame(
    ID          = pred$ID,
    shift       = rank_qnorm(shift),
    disturbance = rank_qnorm(disturbance),
    stringsAsFactors = FALSE
  )
}

# 用 lightgbm 拟合 + SHAP 解释
run_lgb_shap <- function(X, Y, tag, out_dir = "shap") {
  # 对齐 X 和 Y（按 rownames(X) 匹配 pred 的 ID）
  Y <- Y[match(rownames(X), Y$ID), ]
  keep <- !is.na(Y$value)
  X_use <- X[keep, , drop = FALSE]
  y_use <- Y$value[keep]

  dtrain <- lgb.Dataset(data = as.matrix(X_use), label = y_use)

  params <- list(
    objective    = "regression",
    metric       = "mse",
    boosting     = "gbdt",
    num_leaves   = 100,
    learning_rate = 0.05,
    num_threads  = 19
  )

  model <- lgb.train(params, dtrain, nrounds = 1000)

  predictions <- predict(model, as.matrix(X_use))
  perf <- postResample(predictions, y_use)

  # SHAP
  shap_obj <- shapviz(
    model,
    X_pred = as.matrix(X_use),
    predict_function = function(model, newdata) {
      predict(model, newdata)   # lightgbm 不需要 type="response"
    }
  )

  save(shap_obj, file = file.path(out_dir, paste0(tag, "_shap.RData")))

  imp_plot <- sv_importance(shap_obj,
                            kind        = "beeswarm",
                            max_display = 20,
                            show_plot   = FALSE)
  ggsave(file.path(out_dir, paste0(tag, "_importance.png")),
         plot = imp_plot, width = 8, height = 10, dpi = 300)

  list(model = model, perf = perf, shap = shap_obj)
}

# ---------- 7. 主流程：6 个组合 ----------
# 组合：all / g1 / g2 × shift / disturbance
dir.create("shap", showWarnings = FALSE)

# 7.1 读取三组 BART 预测
bart_all <- read_bart_pred("bart/allpred.csv")
bart_g1  <- read_bart_pred("bart/g1pred.csv")
bart_g2  <- read_bart_pred("bart/g2pred.csv")

# 7.2 定义组合表
combos <- list(
  list(tag = "all_shift",       X = X_all, Y = bart_all, col = "shift"),
  list(tag = "all_disturbance", X = X_all, Y = bart_all, col = "disturbance"),
  list(tag = "g1_shift",        X = X_g1,  Y = bart_all,  col = "shift"),
  list(tag = "g1_disturbance",  X = X_g1,  Y = bart_all,  col = "disturbance"),
  list(tag = "g2_shift",        X = X_g2,  Y = bart_all,  col = "shift"),
  list(tag = "g2_disturbance",  X = X_g2,  Y = bart_all,  col = "disturbance")
)

# 7.3 逐个运行
results <- list()
for (combo in combos) {
  message("Running: ", combo$tag)

  Y_df <- data.frame(ID = combo$Y$ID, value = combo$Y[[combo$col]])
  res <- run_lgb_shap(combo$X, Y_df, tag = combo$tag)
  results[[combo$tag]] <- res$perf

  message("  Performance: ",
          paste(names(res$perf), round(res$perf, 4), sep = "=", collapse = ", "))
}

# ---------- 8. 保存性能指标 ----------
perf_table <- do.call(rbind, lapply(names(results), function(tag) {
  data.frame(tag = tag, t(results[[tag]]), stringsAsFactors = FALSE)
}))
fwrite(perf_table, "shap/lgb_performance.csv")
