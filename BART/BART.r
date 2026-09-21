# ============================================================
# all/early/late molecules -> BART 
# ============================================================

library(data.table)
library(stochtree)
library(lubridate)

# ---------- 1. 读取聚类标签 ----------
load("cluster_labels.RData")

g1_features <- names(cluster_labels)[cluster_labels == 1]
g2_features <- names(cluster_labels)[cluster_labels == 2]

# ---------- 2. 读取残差矩阵 ----------
adj_cell  <- fread("adjcell.txt.gz",  data.table = FALSE)
adj_chem  <- fread("adjchem.txt.gz",  data.table = FALSE)
adj_metab <- fread("adjmetab.txt.gz", data.table = FALSE)

# 列名去掉 fread 自动加的 "X" 前缀
strip_X <- function(df) {
  colnames(df) <- gsub("^X", "", colnames(df))
  df
}
adj_cell  <- strip_X(adj_cell)
adj_chem  <- strip_X(adj_chem)
adj_metab <- strip_X(adj_metab)

# 第一列是特征名，作为行名；其余为样本
to_feature_by_sample <- function(df) {
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  t(df)   # 转置为 样本 × 特征
}
adj_cell  <- to_feature_by_sample(adj_cell)
adj_chem  <- to_feature_by_sample(adj_chem)
adj_metab <- to_feature_by_sample(adj_metab)

# ---------- 3. 合并三类组学 ----------
adj_blood <- cbind(adj_cell, adj_chem, adj_metab)
rm(adj_cell, adj_chem, adj_metab)
gc()

# 列名合法化（括号、减号替换为点）
colnames(adj_blood) <- gsub("[()\\-]", ".", colnames(adj_blood))

# ---------- 4. 读取采血时间 ----------
time_data <- fread("bloodtime_participant.csv", data.table = FALSE)
time_data <- time_data[, c(1, 2, 6)]
colnames(time_data) <- c("id", "datetime", "fasting")

parsed_times <- ymd_hms(time_data$datetime)
time_data$hour <- hour(parsed_times) +
                  minute(parsed_times) / 60 +
                  second(parsed_times) / 3600

# 建立 id -> hour 的映射
id_to_hour <- setNames(time_data$hour, time_data$id)

# ---------- 5. 定义 BART 拟合 + 预测函数 ----------
fit_bart <- function(feature_matrix, id_to_hour,
                     num_gfr = 100, num_burnin = 500, num_mcmc = 1200,
                     num_threads = 12, seed = 1234,
                     interval_level = 0.90) {

  # 对齐时间
  sample_ids <- rownames(feature_matrix)
  y <- id_to_hour[sample_ids]
  keep <- !is.na(y)
  feature_matrix <- feature_matrix[keep, , drop = FALSE]
  y <- y[keep]

  # 缺失值用列中位数填补
  X <- as.matrix(feature_matrix)
  X[] <- apply(X, 2, function(col) {
    col[is.na(col)] <- median(col, na.rm = TRUE)
    col
  })

  bart_model <- bart(
    X_train = X,
    y_train = y,
    num_gfr = num_gfr,
    num_burnin = num_burnin,
    num_mcmc = num_mcmc,
    mean_forest_params = list(
      min_samples_leaf = 10,
      alpha            = 0.95,
      max_depth        = 7,
      num_trees        = 70,
      beta             = 2.0,
      num_features_subsample = 40
    ),
    general_params = list(
      num_threads      = num_threads,
      keep_gfr         = FALSE,
      cutpoint_grid_size = 80,
      keep_burnin      = FALSE,
      random_seed      = seed
    )
  )

  # 后验预测抽样（每个样本一个 draw）
  post_preds <- sampleBARTPosteriorPredictive(
    model_object         = bart_model,
    X                    = X,
    num_draws_per_sample = 1
  )

  pred_mean <- apply(post_preds, 1, mean)

  # 均值森林后验
  mf_posterior <- predict(
    object = bart_model,
    X      = X,
    type   = "posterior",
    terms  = "mean_forest"
  )
  pred_sd <- apply(mf_posterior, 1, sd)

  # 后验区间宽度
  intervals <- computeBARTPosteriorInterval(
    model_object = bart_model,
    terms        = "mean_forest",
    X            = X,
    level        = interval_level
  )
  interval_width <- intervals$upper - intervals$lower

  data.frame(
    ID        = rownames(X),
    pred      = pred_mean,
    true      = y,
    sd        = pred_sd,
    interval  = interval_width,
    stringsAsFactors = FALSE
  )
}

# ---------- 6. 三次拟合 ----------


# 6.1 all
pred_all <- fit_bart(adj_blood, id_to_hour)
fwrite(pred_all, "bart/allpred.csv")

# 6.2 early-peak
pred_g1 <- fit_bart(adj_blood[, g1_features, drop = FALSE], id_to_hour)
fwrite(pred_g1, "bart/g1pred.csv")

# 6.3 late-peak
pred_g2 <- fit_bart(adj_blood[, g2_features, drop = FALSE], id_to_hour)
fwrite(pred_g2, "bart/g2pred.csv")
