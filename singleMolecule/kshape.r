library(dtwclust)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

gam_cell <- fread("gam_cell.csv",data.table=F)
gam_chem <- fread("gam_chem.csv",data.table=F)
gam_metab <- fread("gam_metab.csv",data.table=F)

gam_dat <- data.frame(gam_cell,gam_chem[,-1],gam_metab[,-1])
curves_matrix <- gam_dat %>% select(-time) %>%  as.matrix() %>% t() 

curves_list <- split(curves_matrix, row(curves_matrix))
set.seed(1234)   # 保证可重复

ks_results <- tsclust(
  series   = curves_list,
  type     = "partitional",
  k        = 2:7,
  distance = "sbd",                  
  centroid = "shape",                
  preproc  = zscore,                 
  control  = partitional_control(
    nrep     = 15,                   
    iter.max = 300
  ),
  seed     = 1234,
  trace    = TRUE
)

metrics <- sapply(ks_results, function(x) {
  cvi(x, 
      type = c("Sil", "CH", "D"), 
      log.base = 10)
})

metrics_table <- as.data.frame(t(metrics))
colnames(metrics_table) <- c("Silhouette", "Calinski_Harabasz", "Dunn")

best_index <- 5
best_ks <- ks_results[[best_index]]

pdf("kshape.pdf",height=5,width=15)
plot(best_ks, type = "series", 
     main = paste("k-Shape centroid (k =", best_k, ")"),
     xlab = "time", 
     ylab = "Z-score")
dev.off()

cluster_labels <- best_ks@cluster
save(cluster_labels,file="cluster_labels.RData")
