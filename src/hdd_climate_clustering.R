# Province-level HDD profile clustering for Türkiye
# -------------------------------------------------
# Input:  data/hdd_daily_matrix.csv
# Output: figures/ and outputs/ tables

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(cluster)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", args[grep("--file=", args)])
project_dir <- if (length(script_file) > 0) normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE) else getwd()
input_file  <- file.path(project_dir, "data", "hdd_daily_matrix.csv")
output_dir  <- file.path(project_dir, "outputs")
figure_dir  <- file.path(project_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


hdd_wide <- read_csv(input_file, show_col_types = FALSE)
province_names <- names(hdd_wide)[names(hdd_wide) != "date"]

# Rows are provinces and columns are daily HDD observations.
X <- t(as.matrix(hdd_wide[, province_names]))
rownames(X) <- province_names

# Row-wise z-score: clustering captures the shape of each province's winter HDD profile,
# not only the absolute coldness level.
Xz <- t(scale(t(X)))
Xz[is.na(Xz)] <- 0

D <- dist(Xz, method = "euclidean")
hc_single <- hclust(D, method = "single")
hc_complete <- hclust(D, method = "complete")

calc_ch <- function(X, cl) {
  n <- nrow(X)
  k <- length(unique(cl))
  overall <- colMeans(X)
  W <- 0
  B <- 0
  for (g in unique(cl)) {
    Xg <- X[cl == g, , drop = FALSE]
    ng <- nrow(Xg)
    mu <- colMeans(Xg)
    W <- W + sum((Xg - matrix(mu, ng, ncol(X), byrow = TRUE))^2)
    B <- B + ng * sum((mu - overall)^2)
  }
  (B / (k - 1)) / (W / (n - k))
}

k_grid <- 2:10
validation <- tibble(
  k = k_grid,
  calinski_harabasz = sapply(k_grid, function(k) calc_ch(Xz, cutree(hc_complete, k = k))),
  average_silhouette = sapply(k_grid, function(k) mean(silhouette(cutree(hc_complete, k = k), D)[, "sil_width"]))
)
write_csv(validation, file.path(output_dir, "cluster_validation_metrics.csv"))

main_k <- validation$k[which.max(validation$calinski_harabasz)]
cl_k3 <- cutree(hc_complete, k = 3)
cl_k5 <- cutree(hc_complete, k = 5)

set.seed(123)
km3 <- kmeans(Xz, centers = 3, nstart = 50)$cluster
km5 <- kmeans(Xz, centers = 5, nstart = 50)$cluster

cluster_assignments <- tibble(
  province = province_names,
  hierarchical_k3 = as.integer(cl_k3),
  hierarchical_k5 = as.integer(cl_k5),
  kmeans_k3 = as.integer(km3),
  kmeans_k5 = as.integer(km5),
  mean_hdd = rowMeans(X),
  min_hdd = apply(X, 1, min),
  max_hdd = apply(X, 1, max)
)
write_csv(cluster_assignments, file.path(output_dir, "hdd_cluster_assignments.csv"))

summarise_clusters <- function(labels) {
  tibble(province = province_names, cluster = labels) %>%
    group_by(cluster) %>%
    summarise(
      n_provinces = n(),
      mean_hdd = mean(X[province, ]),
      min_hdd = min(X[province, ]),
      max_hdd = max(X[province, ]),
      provinces = paste(province, collapse = ", "),
      .groups = "drop"
    )
}

write_csv(summarise_clusters(cl_k3), file.path(output_dir, "k3_cluster_summary.csv"))
write_csv(summarise_clusters(cl_k5), file.path(output_dir, "k5_cluster_summary.csv"))

png(file.path(figure_dir, "complete_linkage_dendrogram.png"), width = 1600, height = 800, res = 160)
plot(hc_complete, main = "Complete linkage dendrogram - HDD province profiles", xlab = "", sub = "", cex = 0.55)
dev.off()

png(file.path(figure_dir, "single_linkage_dendrogram.png"), width = 1600, height = 800, res = 160)
plot(hc_single, main = "Single linkage dendrogram - HDD province profiles", xlab = "", sub = "", cex = 0.55)
dev.off()

p_ch <- ggplot(validation, aes(k, calinski_harabasz)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  labs(title = "Calinski-Harabasz criterion - HDD province clusters",
       x = "Number of clusters (k)", y = "Calinski-Harabasz score") +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "calinski_harabasz.png"), p_ch, width = 7, height = 4.5, dpi = 160)

p_sil <- ggplot(validation, aes(k, average_silhouette)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  labs(title = "Average silhouette score - HDD province clusters",
       x = "Number of clusters (k)", y = "Average silhouette score") +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "silhouette_score.png"), p_sil, width = 7, height = 4.5, dpi = 160)

profile <- as.data.frame(Xz) |>
  mutate(cluster = factor(cl_k3), province = province_names) |>
  pivot_longer(cols = -c(cluster, province), names_to = "day", values_to = "z_hdd") |>
  group_by(cluster, day) |>
  summarise(mean_z_hdd = mean(z_hdd), .groups = "drop") |>
  mutate(day_index = as.integer(factor(day, levels = names(as.data.frame(Xz)))))

p_profile <- ggplot(profile, aes(day_index, mean_z_hdd, group = cluster, linetype = cluster)) +
  geom_line(linewidth = 0.55) +
  labs(title = "Mean standardized HDD profile by cluster (k = 3)",
       x = "Day index", y = "Mean standardized HDD", linetype = "Cluster") +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "k3_mean_hdd_profiles.png"), p_profile, width = 8, height = 4.8, dpi = 160)

cat("Main k selected by CH:", main_k, "\n")
print(validation)
print(table(Hierarchical = cl_k3, KMeans = km3))
print(table(Hierarchical = cl_k5, KMeans = km5))
