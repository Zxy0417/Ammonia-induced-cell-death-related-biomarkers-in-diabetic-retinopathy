# Integrative data mining and in-silico perturbation reveals UROD-mediated ammonia-induced cell death and immune regulatory networks in diabetic retinopathy
# Data and Code Availability for Biodata Mining
This repository contains the minimal, clean, and reproducible core codes and processed datasets for the manuscript entitled Integrative bioinformatics analysis reveals ammonia-induced cell death-associated molecular patterns and regulatory networks in diabetic retinopathy.
All raw transcriptomic data used in this study are publicly available from the GEO database. No private clinical data or original sequencing files are included here to avoid data leakage.
Environment Requirements
R ≥ 4.2
Required R packages: clusterProfiler, org.Hs.eg.db, limma, DESeq2, GSVA, pROC, glmnet, randomForest, caret, pheatmap, corrplot, ggplot2, CIBERSORT
Repository Content
This repo only retains core reproducible files for open science declaration:
- analysis_clean.R: full analysis script
- gene.txt: intersected core gene list
- AD.csv: ammonia-induced cell death related gene set
- immuno_checkpoint.csv: immune checkpoint gene list
- DR_exp_combined.Rdata: normalized processed expression matrix
