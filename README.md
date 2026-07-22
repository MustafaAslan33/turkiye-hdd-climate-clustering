# Clustering Turkish Provinces by Winter Heating-Degree-Day Profiles

This repository contains a province-level clustering analysis of winter heating-degree-day (HDD) profiles in Türkiye. The purpose is to group provinces that show similar heating-demand dynamics during the winter season, rather than simply ranking provinces by average coldness.

## Project summary

The analysis uses 151 daily HDD observations for 81 provinces between 1 November 2024 and 31 March 2025. Each province is treated as one time-series profile. The daily profile of every province is standardised using a row-wise z-score, and clustering is performed on the standardised profiles using Euclidean distance.

The main model is complete-linkage hierarchical clustering. Single linkage is used as a diagnostic comparison, while k-means is used as a robustness check.

## Main findings

- The Calinski-Harabasz criterion supports a three-cluster structure.
- The silhouette criterion favours a simpler two-cluster split, but the three-cluster result is more useful for energy-planning interpretation.
- The k = 3 solution separates provinces into western/coastal, transition/inland and eastern-highland HDD profile groups.
- The k = 5 solution provides a more detailed regional view, separating Mediterranean, western, southeastern transition, central/northern inland and eastern-highland zones.
- K-means clustering produces broadly consistent results, supporting the presence of stable HDD profile regions.

## Repository structure

```text
data/
  hdd_daily_matrix.csv
  hdd_daily_matrix.xlsx
  hdd_standardized_profiles.csv

figures/
  01_calinski_harabasz.png
  02_silhouette_score.png
  03_complete_linkage_dendrogram.png
  04_single_linkage_dendrogram.png
  05_k3_mean_hdd_profiles.png
  08_map_hierarchical_k3.png
  09_map_hierarchical_k5.png
  10_map_kmeans_k3.png
  11_map_kmeans_k5.png

outputs/
  hdd_cluster_results.xlsx
  hdd_cluster_assignments.csv
  cluster_validation_metrics.csv
  k3_cluster_summary.csv
  k5_cluster_summary.csv

report/
  HDD_Climate_Clustering_Report_EN.pdf
  HDD_Climate_Clustering_Report_EN.docx

src/
  hdd_climate_clustering.R
```

## Methodology

1. Convert the HDD table into an 81 x 151 province-day matrix.
2. Standardise each province's daily profile with a row-wise z-score.
3. Calculate Euclidean distances between standardised province profiles.
4. Compare single-linkage and complete-linkage hierarchical clustering.
5. Evaluate k = 2 to 10 using the Calinski-Harabasz index and average silhouette score.
6. Use k = 3 as the main clustering solution and k = 5 as a supplementary detailed view.
7. Compare hierarchical clustering outputs with k-means results as a robustness check.

## How to run

Install the required R packages:

```r
install.packages(c("dplyr", "readr", "tidyr", "cluster", "ggplot2"))
```

Then run:

```r
source("src/hdd_climate_clustering.R")
```

The script writes updated outputs to `outputs/` and figures to `figures/`.

## Notes

The clustering is descriptive and exploratory. It should not be interpreted as a full electricity-demand forecast. The results are most useful as regional HDD segments that can feed demand modelling, heating-load analysis or energy-efficiency planning.
