###
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggnewscale)
library(limma)
library(pROC)
library(ROCR)
library(caret)
library(glmnet)
library(e1071)
library(randomForest)
library(forestplot)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(ggsci)
library(sva)
library(DESeq2)
library(GEOquery)
library(ggcorrplot)
set.seed(123)
####################
eSet=getGEO("GSE60436",destdir='.')
?save
save(eSet,file="GSE60436_eSet.Rdata")
load('GSE60436_eSet.Rdata')
a=eSet[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
pd=pd[,c(1,2)]
head(pd)
names(pd) <- c('group','sample')
library(stringr)
pd$group <- c(group=str_split(pd$group,'_',simplify=T)[,1])
pd$group <- ifelse(pd$group=='Retina','control','DR')
group_list=pd[,1]
table(group_list)
annoation <- read.table('GPL6884-11607.txt',sep='\t',fill=TRUE,header = TRUE,check.names = FALSE)
ids=annoation[,c(1,13)]
dat=as.data.frame(dat)
dat1 <- merge(ids,dat,by.x=1,by.y=0)
dat1 =dat1[,-1]
dat1 =dat1[!duplicated(dat1$Symbol),]
rownames(dat1) <- dat1$Symbol
dat1=dat1[,-1]
dat1=na.omit(dat1)
pdf('GSE60436boxplot.pdf', width = 8, height = 7)
boxplot(dat1)
dev.off()
save(dat1,file='GSE60436.Rdata')
save(dat1,pd,group_list,file='GSE60436_dat_group.Rdata')

#
library(data.table)
exp <- fread("GSE102485_raw_counts_GRCh38.p13_NCBI.tsv.gz")
head(exp)
library(stringi)
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
options(connectionObserver = NULL)
library(BiocManager)
name <- bitr(exp$GeneID,fromType='ENTREZID',toType='SYMBOL',OrgDb='org.Hs.eg.db',drop=T)
name$ENTREZID=as.numeric(name$ENTREZID)
exp=right_join(name,exp,by=c('ENTREZID'='GeneID'))
exp=exp[,-1]
exp <- aggregate(.~SYMBOL,FUN=mean,data=exp)
library(GEOquery)
f='GSE102485_eSet.Rdata'
if(!file.exists(f)){
  gset <- getGEO('GSE102485',destdir = ".",
                 AnnotGPL = F,
                 getGPL = F)
  save(gset,file=f)
}
load('GSE102485_eSet.Rdata')
b=gset[[1]]
pd=pData(b)
pd=pd[,c(2,10)]
head(pd)
names(pd) <- c('sample','group')
library(dplyr)
pd1 <- pd %>%
  filter(pd$group %in% c("disease: abbreviations of type II diabetes", "disease: Normal retina","disease: abbreviations of type I diabetes"))
pd1$group=ifelse(pd1$group=="disease: abbreviations of type II diabetes","DR",pd1$group)
pd1$group=ifelse(pd1$group=="disease: abbreviations of type I diabetes","DR",pd1$group)
pd1$group=ifelse(pd1$group=='DR','DR','Control')
print(table(pd1$group))
group_list=pd1[,2]
rownames(exp)<-exp$SYMBOL
exp=exp[,-1]
exp=as.data.frame(t(exp))
exp=merge(pd1,exp,by.x=1,by.y=0)
rownames(exp)<-exp$sample
exp=exp[,-c(1,2)]
setdiff(rownames(exp), pd1$sample)
pd1=merge(pd1,exp,by.x=1,by.y=0)
rownames(pd1)<-pd1$sample
pd1=pd1[,c(1,2)]
exp=as.data.frame(t(exp))
exp1=log(exp+1)
pdf('GSE102485boxplot.pdf', width = 8, height = 7)
boxplot(exp1)
dev.off()
exp1=normalizeBetweenArrays(exp1)
group_list=pd1[,2]
save(exp1,pd1,group_list,file = 'GSE102485_dat_group.Rdata')
exp<-dat1
rm(dat1,pd,group_list)
rm(pd1,group_list)
exp1=as.data.frame(exp1)
comm_gene <- intersect(rownames(exp),rownames(exp1))
exp=exp[comm_gene,]
exp1=exp1[comm_gene,]
pdf('GSE102485exp1.pdf', width = 8, height = 7)
boxplot(exp1)
dev.off()
pdf('GSE60436exp.pdf', width = 8, height = 7)
boxplot(exp)
dev.off()
exp2 = cbind(exp,exp1)
pdf('bindboxlot.pdf', width = 8, height = 7)
boxplot(exp2)
dev.off()
graphics.off()
exp2=normalizeBetweenArrays(exp2)
save(exp2,file = "DR_exp.Rdata")
rm(list = ls())
library(sva)
library(limma)
#ComBat
load("GSE60436_dat_group.Rdata"); g1 <- group_list
load("GSE102485_dat_group.Rdata"); g2 <- group_list
group_list <- c(g1, g2)
gse<-c(rep('GSE60436',9),rep('GSE102485',24))
table(group_list,gse)
load("DR_exp.Rdata")
dat<-exp2
dat[1:4,1:4]
betch<-c(rep('GSE60436',9),rep('GSE102485',24))
design=model.matrix(~group_list)
exp2<- removeBatchEffect(dat,batch = betch,design = design)
dim(exp2) 
exp2=normalizeBetweenArrays(exp2)
graphics.off()
pdf('boxplot2.pdf', width = 8, height = 7)
boxplot(exp2)
dev.off()
save(exp2,group_list,file = "DR_exp_combined.Rdata")
#differetial analysis
library(limma)
library(dplyr)
library(ggplot2)
load("DR_exp_combined.Rdata")
group_list1 <- as.data.frame(group_list)
design<-model.matrix(~0+factor(group_list1$group_list))
colnames(design)<-levels(factor(group_list1$group_list))
rownames(design)<-colnames(exp2)
contrast.matrix<-makeContrasts(DR-Control,levels = design ) 
fit<-lmFit(exp2,design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
DEG<-topTable(fit2,coef = 1,n=Inf)
DEG$regulate<-ifelse(DEG$P.Value<=0.05 & abs(DEG$logFC)>=1,
                     ifelse(DEG$logFC>0, 'Up', 'Down'), 'Stable')
table(DEG$regulate)          
write.csv(data.frame(gene_symbol=rownames(DEG),DEG),file='DEG.csv')
pdf('DEG-VOLCANO.pdf')
library(ggrepel)
DEG$symbol <- rownames(DEG)

top10 <- DEG[order(DEG$P.Value), ][1:10, ]
pdf('DEG-VOLCANO.pdf')
ggplot(DEG, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6, size = 3.5, aes(color = regulate)) +
  ylab('-log10(P.Value)') +
  scale_color_manual(values = c("#2166ac", "darkgrey", "#c51b7d")) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = 'black', lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = 'black', lwd = 0.8) +
  geom_text_repel(data = top10, aes(label = symbol),
                  size = 3.5, box.padding = 0.5, max.overlaps = 10) +
  theme_bw()
while (!is.null(dev.list())) dev.off()
############################################################
# 1. Intersection: DR DEGs x ammonia-induced cell death genes
############################################################
DR <- read.csv("DEG.csv", stringsAsFactors = FALSE)
DR$regulate[DR$regulate == "Stable"] <- NA
DR <- na.omit(DR)

AD <- read.csv("ammonia_icd_gene_set.csv", stringsAsFactors = FALSE)

datalist <- list(DR = DR$gene_symbol, AD = AD$Gene.Symbol)
comm_gene <- Reduce(intersect, list(DR$gene_symbol, AD$Gene.Symbol))
write.table(comm_gene, file = "gene.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

# Venn
pdf("Venn_DR_ammonia.pdf", width = 8, height = 6)
ggvenn(datalist,
       fill_color = c("lightblue", "pink"),
       fill_alpha = 0.5,
       stroke_linetype = "longdash",
       set_name_size = 3,
       text_size = 4)
dev.off()

############################################################
# 2. GSEA: ammonia-induced ICD pathway
############################################################
DEG <- read.csv("DEG.csv", stringsAsFactors = FALSE)

DEG$ENTREZID <- mapIds(org.Hs.eg.db,
  keys = DEG$gene_symbol, column = "ENTREZID",
  keytype = "SYMBOL", multiVals = "first")

geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- geneList[!is.na(names(geneList))]
geneList <- sort(geneList, decreasing = TRUE)

entrez_raw <- mapIds(org.Hs.eg.db,
  keys = unique(trimws(AD$Gene.Symbol)), column = "ENTREZID",
  keytype = "SYMBOL", multiVals = "first")

ammonia_geneset <- data.frame(
  term = "Ammonia_ICD",
  gene = as.character(entrez_raw),
  stringsAsFactors = FALSE
)

set.seed(2024)
gsea_res <- GSEA(
  geneList      = geneList,
  TERM2GENE     = ammonia_geneset,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  minGSSize     = 5,
  maxGSSize     = 1000,
  eps           = 0,
  seed          = 2024,
  verbose       = FALSE
)

gsea_res@result[, c("Description", "NES", "pvalue", "p.adjust", "qvalue")]

pdf("GSEA_ammonia.pdf", width = 8, height = 6)
gseaplot2(gsea_res, geneSetID = "Ammonia_ICD",
  title = "GSEA: Ammonia-induced ICD in DR",
  pvalue_table = TRUE, ES_geom = "line")
dev.off()

############################################################
# 3. ssGSEA pathway score per sample
############################################################
my_ssgsea <- function(expr_mat, gene_sets, alpha = 0.25, normalization = TRUE) {
  require(matrixStats)
  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)
  es <- matrix(0, nrow = length(gene_sets), ncol = n_samples)
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(expr_mat)
  for (i in seq_along(gene_sets)) {
    gs_genes <- intersect(gene_sets[[i]], rownames(expr_mat))
    if (length(gs_genes) == 0) next
    gene_in_set <- rownames(expr_mat) %in% gs_genes
    for (j in seq_len(n_samples)) {
      expr_sample <- expr_mat[, j]
      order_idx <- order(expr_sample, decreasing = TRUE)
      sorted_in_set <- gene_in_set[order_idx]
      p_vals <- expr_sample[order_idx]
      p_vals[sorted_in_set] <- abs(p_vals[sorted_in_set])^alpha
      p_vals[!sorted_in_set] <- 0
      n_total <- sum(p_vals)
      if (n_total == 0) next
      hit_steps <- p_vals / n_total
      miss_steps <- rep(1 / (n_genes - sum(sorted_in_set)), n_genes)
      miss_steps[sorted_in_set] <- 0
      cum_hit <- cumsum(hit_steps)
      cum_miss <- cumsum(miss_steps)
      diff_vec <- cum_hit - cum_miss
      es[i, j] <- diff_vec[which.max(abs(diff_vec))]
    }
  }
  if (normalization && n_samples > 1) {
    for (i in seq_len(nrow(es))) {
      es[i, ] <- (es[i, ] - mean(es[i, ])) / sd(es[i, ])
    }
  }
  return(es)
}

load("DR_exp_combined.Rdata")
genes_use <- intersect(unique(trimws(AD$Gene.Symbol)), rownames(exp2))
gs_list <- list(Ammonia_ICD = genes_use)
ssgsea_scores <- my_ssgsea(as.matrix(exp2), gs_list)
score <- as.numeric(ssgsea_scores["Ammonia_ICD", ])
names(score) <- colnames(exp2)

df_score <- data.frame(Score = score, Group = group_list)

pdf("ssGSEA_score_boxplot.pdf", width = 10, height = 8)
ggplot(df_score, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_classic() +
  labs(title = "ssGSEA: Ammonia-induced ICD pathway activity")
dev.off()

wilcox.test(Score ~ Group, data = df_score)

roc_obj <- roc(df_score$Group, df_score$Score)
pdf("ssGSEA_score_ROC.pdf", width = 10, height = 8)
plot(roc_obj, main = paste0("AUC = ", round(auc(roc_obj), 3)))
dev.off()

############################################################
# 4. GO enrichment on intersection genes
############################################################
rt <- read.table("gene.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(rt) <- "symbol"
entrezIDs <- mget(rt$symbol, org.Hs.egSYMBOL2EG, ifnotfound = NA)
rt$entrezID <- as.character(entrezIDs)
rt <- rt[!is.na(rt$entrezID), ]
gene <- unique(rt$entrezID)

# Prepare fold-change named by ENTREZID
aflogfc <- DEG$logFC
names(aflogfc) <- DEG$ENTREZID

kk_go <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db,
  pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)

GO <- as.data.frame(kk_go)
GO <- GO[GO$pvalue < 0.05, ]
write.table(GO, file = "GO.xls", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("GO_barplot.pdf", width = 9, height = 7)
barplot(kk_go, showCategory = 7, split = "ONTOLOGY", color = "qvalue") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))
dev.off()

pdf("GO_bubble.pdf", width = 9, height = 7)
dotplot(kk_go, showCategory = 7, orderBy = "GeneRatio", color = "qvalue") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))
dev.off()

pdf("GO_cnet.pdf", width = 10, height = 8)
af_go <- setReadable(kk_go, "org.Hs.eg.db", "ENTREZID")
cnetplot(af_go, foldChange = aflogfc, showCategory = 5,
  categorySize = "pvalue", circular = TRUE, colorEdge = TRUE,
  cex_label_category = 0.65, cex_label_gene = 0.6)
dev.off()

############################################################
# 5. KEGG enrichment (NO rm(list=ls()) -- keeps DEG alive)
############################################################
library(pathview)

kk_kegg <- enrichKEGG(gene = gene, organism = "hsa",
  pvalueCutoff = 1, qvalueCutoff = 1)
KEGG <- as.data.frame(kk_kegg)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x)
  paste(rt$symbol[match(strsplit(x, "/")[[1]], rt$entrezID)], collapse = "/")))
KEGG <- KEGG[KEGG$pvalue < 0.05, ]
write.table(KEGG, file = "KEGG.xls", sep = "\t", quote = FALSE, row.names = FALSE)

showNum <- min(20, nrow(KEGG))

pdf("KEGG_barplot.pdf", width = 9, height = 7)
barplot(kk_kegg, showCategory = showNum, color = "pvalue") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))
dev.off()

pdf("KEGG_bubble.pdf", width = 9, height = 7)
dotplot(kk_kegg, showCategory = showNum, orderBy = "GeneRatio", color = "pvalue") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))
dev.off()

pdf("KEGG_cnet.pdf", width = 10, height = 8)
af_kegg <- setReadable(kk_kegg, "org.Hs.eg.db", "ENTREZID")
cnetplot(af_kegg, foldChange = aflogfc, showCategory = 5,
  categorySize = "pvalue", circular = TRUE, colorEdge = TRUE,
  cex_label_category = 0.65, cex_label_gene = 0.6)
dev.off()

############################################################
# 6. Machine learning: LASSO, SVM-RFE, RF, Logistic
############################################################
source("./msvmRFE.R")

gene_df <- read.table("gene.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(gene_df) <- "gene"
exp <- as.data.frame(exp2[gene_df$gene, ])
exp_t <- t(exp)
y <- ifelse(group_list == "Control", 0, 1)

# --- LASSO ---
set.seed(123456789)
x_lasso <- as.matrix(exp_t)
fit_lasso <- glmnet(x_lasso, y, family = "binomial", alpha = 1, lambda = NULL)
cvfit <- cv.glmnet(x_lasso, y, family = "binomial", alpha = 1,
  type.measure = "deviance", nfolds = 5)
lasso_coefs <- coef(cvfit, s = "lambda.min")
lasso_features <- lasso_coefs@Dimnames[[1]][which(lasso_coefs != 0)][-1]
write.csv(data.frame(x = lasso_features), "feature_lasso.csv", row.names = FALSE)

# --- SVM-RFE ---
train_df <- data.frame(y = factor(y), exp_t)
colnames(train_df)[1] <- "group"
input <- train_df
nfold <- 5
nrows <- nrow(input)
folds <- rep(1:nfold, length.out = nrows)[sample(nrows)]
folds <- lapply(folds, function(x) which(folds == x))
results <- lapply(folds, svmRFE.wrap, input, k = 5, halve.above = 20)
top_features <- WriteFeatures(results, input, save = FALSE)
write.csv(top_features, "feature_svm.csv", row.names = FALSE)

FeatSweep <- lapply(1:20, FeatSweep.wrap, results, input)
no.info <- min(prop.table(table(input[, 1])))
errors <- sapply(FeatSweep, function(x) if (is.null(x)) NA else x$error)

min_error_index <- which.min(errors)
SVMRFEgenes <- top_features[1:min_error_index, "FeatureName"]
write.csv(SVMRFEgenes, file = "SVMRFEgenes.csv", row.names = FALSE)

# --- Random Forest ---
set.seed(123456)
exp_rf <- as.data.frame(exp2[gene_df$gene, ])
rownames(exp_rf) <- str_replace_all(rownames(exp_rf), "-", "\\.")
exp_rf <- t(exp_rf)
rf <- randomForest(factor(group_list) ~ ., data = as.data.frame(exp_rf), ntree = 500)
optionTrees <- which.min(rf$err.rate[, 1])
rf2 <- randomForest(factor(group_list) ~ ., data = as.data.frame(exp_rf), ntree = optionTrees)

importance_rf <- as.data.frame(importance(rf2, type = 1))
colnames(importance_rf) <- "importance"
importance_rf$Gene <- str_replace_all(rownames(importance_rf), "\\.", "-")
rfGenes <- importance_rf[order(-importance_rf$importance), ]
write.csv(rfGenes, "rf.importance.csv", row.names = FALSE)

# --- Univariate Logistic ---
exp_glm <- as.data.frame(exp2[gene_df$gene, ])
rownames(exp_glm) <- str_replace_all(rownames(exp_glm), "-", "\\.")
exp_glm <- t(exp_glm)
data_glm <- data.frame(y = factor(y), exp_glm)
colnames(data_glm)[1] <- "group"

Uni_glm_model <- function(var, data) {
  FML <- as.formula(paste0("group ~ ", var))
  glm1 <- glm(FML, family = binomial, data = data)
  glm2 <- summary(glm1)
  OR <- round(exp(coef(glm1)), 2)
  SE <- glm2$coefficients[, 2]
  CI_Lower <- round(exp(coef(glm1) - 1.96 * SE), 2)
  CI_Upper <- round(exp(coef(glm1) + 1.96 * SE), 2)
  CI <- paste0(CI_Lower, "-", CI_Upper)
  P <- signif(glm2$coefficients[, 4], 3)
  data.frame(characteristics = var, OR = OR, CI = CI, P = P)[-1, ]
}

variable_names <- colnames(data_glm)[-1]
Uni_glm <- do.call(rbind, lapply(variable_names, Uni_glm_model, data = data_glm))
Uni_glm$characteristics <- str_replace_all(Uni_glm$characteristics, "\\.", "-")
write.csv(Uni_glm, file = "Uni_glm.csv", row.names = FALSE)

# Forest plot
fp <- Uni_glm
fp <- fp[fp$P < 0.05 & as.numeric(sub("-.*$", "", fp$CI)) > 0, ]
fp$OR_mean_log <- log(as.numeric(fp$OR))
fp$OR_1_log <- log(as.numeric(sub("-.*$", "", fp$CI)))
fp$OR_2_log <- log(as.numeric(sub(".*-", "", fp$CI)))

fp_sorted <- fp[order(fp$OR_mean_log), ]
fp_sorted$Var <- factor(fp_sorted$characteristics, levels = fp_sorted$characteristics)

label_text <- cbind(
  fp_sorted$characteristics,
  paste0(fp_sorted$OR, " (", fp_sorted$CI, ")"),
  ifelse(fp_sorted$P < 0.001, "<0.001",
         ifelse(fp_sorted$P < 0.01, format(fp_sorted$P, digits = 2),
                format(fp_sorted$P, digits = 3)))
)

pdf("forestplot_logistic.pdf", width = 8, height = 8, onefile = FALSE)
forestplot(labeltext = label_text,
  mean = c(NA, fp_sorted$OR_mean_log),
  lower = c(NA, fp_sorted$OR_1_log),
  upper = c(NA, fp_sorted$OR_2_log),
  zero = 0,
  xlim = c(min(fp_sorted$OR_1_log) * 1.1, max(fp_sorted$OR_2_log) * 1.1),
  align = c("l", "c", "c"),
  boxsize = 0.3,
  lineheight = unit(5, "mm"),
  colgap = unit(5, "mm"),
  lwd.zero = 1.5, lwd.ci = 2,
  col = fpColors(box = "red", summary = "blue", lines = "black", zero = "#458B00"),
  xlab = "log(OR)",
  txt_gp = fpTxtGp(ticks = gpar(cex = 0.85), xlab = gpar(cex = 0.8), cex = 0.9),
  lty.ci = "solid",
  title = "Forest plot of logistic regression",
  new_page = TRUE)
dev.off()

write.csv(fp_sorted$characteristics, file = "LogisticGenes.csv", row.names = FALSE)

############################################################
# 7. Intersection of ML methods
############################################################
Lasso <- read.csv("feature_lasso.csv", stringsAsFactors = FALSE)
SVMRFE <- read.csv("SVMRFEgenes.csv", stringsAsFactors = FALSE)
RF <- read.csv("rf.importance.csv", stringsAsFactors = FALSE)
Logistic <- read.csv("LogisticGenes.csv", stringsAsFactors = FALSE)

datalist_ml <- list(Lasso = Lasso$x, RF = RF$Gene, Logistic = Logistic$x)

pdf("ML_intersection_venn.pdf", width = 8, height = 6)
ggvenn(datalist_ml,
  fill_color = c("blue", "yellow", "green"),
  fill_alpha = 0.5,
  stroke_linetype = "longdash",
  set_name_size = 3, text_size = 4)
dev.off()

common_genes <- Reduce(intersect, list(Lasso$x, RF$Gene, Logistic$x))
write.csv(common_genes, "ML_common_genes.csv", row.names = FALSE)

############################################################
# 8. Core genes x ssGSEA score (Spearman)
############################################################
core_genes <- c("PPP1R12C", "UROD", "PPIG", "SLC1A3")
core_expr <- exp2[core_genes, , drop = FALSE]

cor_list <- lapply(core_genes, function(g)
  cor.test(as.numeric(core_expr[g, ]), score, method = "spearman"))
cor_df <- data.frame(
  Gene = core_genes,
  Rho = sapply(cor_list, `[[`, "estimate"),
  Pval = sapply(cor_list, `[[`, "p.value"))
cor_df$FDR <- p.adjust(cor_df$Pval, method = "BH")
print(cor_df)

pdf("UROD_vs_ammonia_score.pdf", width = 8, height = 6)
ggplot(data.frame(Score = score,
                  UROD = as.numeric(core_expr["UROD", ]),
                  Group = group_list),
  aes(x = UROD, y = Score)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_classic() +
  labs(x = "UROD expression", y = "Ammonia-ICD ssGSEA score")
dev.off()

pdf("PPIG_vs_ammonia_score.pdf", width = 8, height = 6)
ggplot(data.frame(Score = score,
                  PPIG = as.numeric(core_expr["PPIG", ]),
                  Group = group_list),
  aes(x = PPIG, y = Score)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_classic() +
  labs(x = "PPIG expression", y = "Ammonia-ICD ssGSEA score")
dev.off()

cor_matrix_genes_score <- cor(t(core_expr), score, method = "spearman")
pdf("core_genes_ssGSEA_heatmap.pdf", width = 8, height = 6)
pheatmap(cor_matrix_genes_score, cluster_rows = TRUE, cluster_cols = FALSE,
  display_numbers = TRUE, number_format = "%.3f")
dev.off()

############################################################
# 9. Internal ROC and expression boxplots
############################################################
load("GSE102485_dat_group.Rdata")
load("GSE60436_dat_group.Rdata")
pd60436 <- pd[, c("sample", "group")]
all_cli <- rbind(pd60436, pd1)
all_cli$group[all_cli$group == "control"] <- "Control"

gene_res <- read.csv("ML_common_genes.csv", stringsAsFactors = FALSE)
exp_roc <- exp2[gene_res$x, ] %>% t() %>% as.data.frame()
exp_roc <- cbind(all_cli[, 2], exp_roc)
colnames(exp_roc)[1] <- "group"

# 5-fold CV logistic ROC
consensus_genes <- c("PPP1R12C", "UROD", "PPIG", "SLC1A3")
cv_data <- exp_roc[, c("group", consensus_genes)]
cv_data$group <- factor(cv_data$group, levels = c("Control", "DR"))

ctrl <- trainControl(method = "cv", number = 5,
  savePredictions = "final", classProbs = TRUE,
  summaryFunction = twoClassSummary)
set.seed(2024)
cv_fit <- train(group ~ ., data = cv_data, method = "glm",
  trControl = ctrl, metric = "ROC")
cv_preds <- cv_fit$pred[order(cv_fit$pred$rowIndex), ]

roc_cv <- roc(response = cv_preds$obs, predictor = cv_preds$DR,
  levels = c("Control", "DR"))
pdf("CV_ROC_4genes.pdf", width = 6, height = 6)
plot(roc_cv, col = "blue", lwd = 2, legacy.axes = TRUE,
  print.auc = TRUE, main = "5-fold CV ROC (4 consensus genes)")
dev.off()

cv_data$PRS_Score <- predict(cv_fit, newdata = cv_data, type = "prob")[, "DR"]
p_density <- ggplot(cv_data, aes(x = PRS_Score, fill = group)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "DR" = "#E64B35")) +
  theme_bw(base_size = 14) +
  labs(x = "PRS score", y = "Density") +
  stat_compare_means(label = "p.format")
ggsave("PRS_density.pdf", p_density, width = 6, height = 5)

# Individual gene ROC (FIXED: SLC1A3 now uses rocSLC1A3)
for (g in consensus_genes) {
  r <- roc(exp_roc$group, exp_roc[[g]], quiet = TRUE)
  pdf(paste0("ROC_", g, ".pdf"), width = 5, height = 5)
  plot(r, col = "red", legacy.axes = TRUE, print.auc = TRUE,
       main = paste0(g, " ROC"))
  legend("bottomright", legend = paste0("AUC = ", round(auc(r), 3)),
         bty = "n")
  dev.off()
}

# Expression boxplots per gene
for (g in consensus_genes) {
  pdf(paste0("Boxplot_", g, ".pdf"), width = 5, height = 5)
  p <- ggplot(exp_roc, aes(x = group, y = .data[[g]], fill = group)) +
    geom_boxplot(alpha = 1) + geom_jitter(width = 0.1, alpha = 0.6) +
    scale_fill_manual(values = c("blue", "red")) +
    theme_classic() +
    stat_compare_means(method = "wilcox.test") +
    labs(y = paste0(g, " expression"))
  print(p)
  dev.off()
}

############################################################
# 10. CIBERSORT immune deconvolution
############################################################
data(LM22, package = "IOBR")
write.table(LM22, file = "LM22.txt", sep = "\t", quote = FALSE,
  row.names = TRUE, col.names = TRUE)
write.table(exp2, file = "exp.txt", sep = "\t", quote = FALSE,
  row.names = TRUE, col.names = TRUE)

res_cibersort <- cibersort("LM22.txt", "exp.txt", perm = 10, QN = FALSE)
save(res_cibersort, file = "res_cibersort.Rdata")

ciber_res <- res_cibersort[, abs(colSums(res_cibersort)) > 0]
ciber_res$group <- all_cli$group

ciber_long <- ciber_res %>%
  select(`P-value`, Correlation, RMSE, everything()) %>%
  pivot_longer(-c(1:4), names_to = "cell_type", values_to = "fraction") %>%
  mutate(cell_type = gsub("_CIBERSORT", "", cell_type),
         cell_type = gsub("_", " ", cell_type))

palette4 <- colorRampPalette(ggsci::pal_npg()(10))(22)
pdf("CIBERSORT_stack.pdf", width = 10, height = 8)
ggplot(ciber_long, aes(x = ID, y = fraction)) +
  geom_bar(stat = "identity", position = "fill", aes(fill = cell_type)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = palette4, name = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom")
dev.off()

pdf("CIBERSORT_boxplot.pdf", width = 10, height = 6)
ggplot(ciber_long, aes(cell_type, fraction, fill = group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  labs(x = NULL, y = "Estimated proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "kruskal.test")
dev.off()

# Core genes x immune cells
exp_core_t <- exp2[gene_res$x, ] %>% t() %>% as.data.frame()
ciber <- ciber_res[, 1:22]
cor_immune <- cor(ciber, exp_core_t)
pdf("CIBERSORT_coregenes_corr.pdf", width = 10, height = 6)
corrplot(t(cor_immune), method = "color",
  col = colorRampPalette(c("#01468b", "white", "#ee0000"))(100),
  addCoef.col = "black", tl.col = "black",
  number.cex = 0.5, tl.cex = 0.7)
dev.off()

############################################################
# 11. Immune checkpoint correlation
############################################################
checkpoint <- read.csv("immuno_checkpoint.csv", stringsAsFactors = FALSE)
checkpoint_exp <- exp2[checkpoint$Symbol, ] %>% na.omit() %>% t() %>% as.data.frame()
cor1 <- cor(checkpoint_exp, exp_core_t)

table3_data <- data.frame(
  Gene = c("UROD", "UROD", "PPIG"),
  Checkpoint = c("LAG3", "CTLA4", "CD40"),
  Spearman_R = c(cor1["LAG3", "UROD"], cor1["CTLA4", "UROD"], cor1["CD40", "PPIG"]),
  Note = c("Immune brake dismantled", "Immune brake dismantled", "Ag presentation blocked")
)
write.csv(table3_data, "Table3_checkpoint_correlation.csv", row.names = FALSE)

pdf("checkpoint_correlation.pdf", width = 10, height = 6)
corrplot(t(cor1), method = "color",
  col = colorRampPalette(c("#01468b", "white", "#ee0000"))(100),
  addCoef.col = "black", tl.col = "black",
  number.cex = 0.5, tl.cex = 0.7)
dev.off()

############################################################
# 12. ssGSEA score x immune cells
############################################################
gsva_score <- score
names(gsva_score) <- colnames(exp2)

immune_frac <- res_cibersort[, 1:22]
common_samples <- intersect(colnames(exp2), rownames(immune_frac))
exp2 <- exp2[, common_samples, drop = FALSE]
immune_frac <- immune_frac[common_samples, , drop = FALSE]
gsva_score <- gsva_score[common_samples]
names(group_list) <- colnames(exp2)
group_list <- group_list[common_samples]

# CIBERSORT differential table
cibersort_diff <- data.frame(Cell = colnames(immune_frac))
cibersort_diff$Control_Mean <- apply(immune_frac, 2, function(x)
  mean(x[group_list == "Control"]))
cibersort_diff$DR_Mean <- apply(immune_frac, 2, function(x)
  mean(x[group_list == "DR"]))
cibersort_diff$Pval <- apply(immune_frac, 2, function(x)
  wilcox.test(x ~ group_list)$p.value)
cibersort_diff$FDR <- p.adjust(cibersort_diff$Pval, method = "BH")
cibersort_diff <- cibersort_diff[order(cibersort_diff$Pval), ]
write.csv(cibersort_diff, "Table2_CIBERSORT_differential.csv", row.names = FALSE)

# Core genes vs pathway score
core_expr <- exp2[core_genes, , drop = FALSE]
gene_path_cor <- apply(core_expr, 1, function(x)
  cor.test(x, gsva_score, method = "spearman"))
gene_path_df <- data.frame(
  Gene = rownames(core_expr),
  Rho = sapply(gene_path_cor, `[[`, "estimate"),
  Pval = sapply(gene_path_cor, `[[`, "p.value"))
gene_path_df$FDR <- p.adjust(gene_path_df$Pval, method = "BH")
write.csv(gene_path_df, "GSVA_spearman_coregenes.csv", row.names = FALSE)

# Selected key immune cells
selected_cells <- c("Monocytes", "Mast cells resting", "Macrophages M1",
                    "T cells CD8", "B cells naive")
combined_mat <- cbind(
  UROD = as.numeric(core_expr["UROD", ]),
  PPIG = as.numeric(core_expr["PPIG", ]),
  GSVA = gsva_score,
  immune_frac[, selected_cells, drop = FALSE])
cor_matrix_final <- cor(combined_mat, method = "spearman")

pdf("CoreGenes_ImmuneCorrelation_heatmap.pdf", width = 8, height = 6)
pheatmap(cor_matrix_final, cluster_rows = TRUE, cluster_cols = TRUE,
  display_numbers = TRUE, number_format = "%.2f",
  color = colorRampPalette(c("steelblue", "white", "#ee0000"))(100),
  main = "Core genes vs key immune cells")
dev.off()

# Scatter: GSVA p-value vs CIBERSORT p-value
gsva_immune_cor <- apply(immune_frac, 2, function(x)
  cor.test(gsva_score, x, method = "spearman"))
gsva_immune_df <- data.frame(
  Cell = colnames(immune_frac),
  Rho = sapply(gsva_immune_cor, `[[`, "estimate"),
  Pval = sapply(gsva_immune_cor, `[[`, "p.value"))

compare_df <- merge(cibersort_diff[, c("Cell", "Pval")],
  gsva_immune_df, by = "Cell")
compare_df$negLog_diff <- -log10(compare_df$Pval.x + 1e-10)
compare_df$negLog_gsva <- -log10(compare_df$Pval.y + 1e-10)

pdf("GSVA_vs_CIBERSORT_scatter.pdf", width = 8, height = 7)
ggplot(compare_df, aes(x = negLog_diff, y = negLog_gsva)) +
  geom_point(aes(color = Pval.y < 0.05), size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("grey80", "#c51b7d"), name = NULL) +
  geom_text_repel(aes(label = Cell), size = 3, max.overlaps = 10) +
  labs(x = "CIBERSORT group difference (-log10 P)",
       y = "GSVA-immune correlation (-log10 P)") +
  theme_bw()
dev.off()

############################################################
# 13. External validation: GSE221521
############################################################
f_matrix <- "GSE221521_series_matrix.txt.gz"
if (!file.exists(f_matrix)) {
  gset <- getGEO("GSE221521", destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
  save(gset, file = "GSE221521_eSet.Rdata")
} else {
  load("GSE221521_eSet.Rdata")
}
b <- gset[[1]]
pd_ext <- pData(b)

group_info <- pd_ext$title
group_full <- case_when(
  grepl("Control", group_info, ignore.case = TRUE) ~ "Control",
  grepl("DM", group_info, ignore.case = TRUE) ~ "DM",
  grepl("DR", group_info, ignore.case = TRUE) ~ "DR",
  TRUE ~ NA_character_)
keep <- group_full %in% c("Control", "DR")
group_ext <- factor(group_full[keep], levels = c("Control", "DR"))
names(group_ext) <- rownames(pd_ext)[keep]

counts_file <- "GSE221521_raw_counts_GRCh38.p13_NCBI.tsv.gz"
raw_counts <- read.table(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
common_samples_ext <- intersect(colnames(raw_counts), names(group_ext))
raw_counts <- raw_counts[, common_samples_ext, drop = FALSE]
group_ext <- group_ext[common_samples_ext]

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
  colData = data.frame(group = group_ext), design = ~ group)
vsd <- vst(dds, blind = FALSE)
expr_norm <- assay(vsd)

rownames_clean <- sub("\\..*", "", rownames(expr_norm))
symbols_ext <- mapIds(org.Hs.eg.db, keys = rownames_clean,
  keytype = "ENSEMBL", column = "SYMBOL")
keep_sym <- !is.na(symbols_ext)
expr_norm <- expr_norm[keep_sym, ]
symbols_ext <- symbols_ext[keep_sym]
rownames(expr_norm) <- symbols_ext
expr_agg <- aggregate(expr_norm, by = list(Gene = rownames(expr_norm)), FUN = mean)
rownames(expr_agg) <- expr_agg$Gene
expr_symbol_ext <- as.matrix(expr_agg[, -1])

# Only 4 consensus genes (FIXED: removed RPE65, HSPA4)
core_genes_ext <- c("PPP1R12C", "UROD", "PPIG", "SLC1A3")
found_genes_ext <- intersect(core_genes_ext, rownames(expr_symbol_ext))
cat("External dataset found:", found_genes_ext, "\n")
expr_final <- expr_symbol_ext[found_genes_ext, , drop = FALSE]
group_ext <- group_ext[colnames(expr_final)]

# Limma differential expression
design_ext <- model.matrix(~ group_ext)
fit_ext <- lmFit(expr_final, design_ext)
fit_ext <- eBayes(fit_ext)
DEG_core_ext <- topTable(fit_ext, coef = 2, number = Inf)
DEG_core_ext$regulate <- ifelse(DEG_core_ext$adj.P.Val < 0.05,
  ifelse(DEG_core_ext$logFC > 0, "Up", "Down"), "Stable")
write.csv(data.frame(Gene = rownames(DEG_core_ext), DEG_core_ext),
  "DEG_core_validation_GSE221521.csv", row.names = FALSE)

# ROC per gene
roc_data_ext <- as.data.frame(t(expr_final))
roc_data_ext$group <- group_ext
auc_summary <- data.frame()

for (g in rownames(expr_final)) {
  r <- roc(roc_data_ext$group, roc_data_ext[[g]],
           levels = c("Control", "DR"), quiet = TRUE)
  auc_val <- round(auc(r), 3)
  ci_val <- round(ci(r), 3)
  auc_summary <- rbind(auc_summary, data.frame(
    Gene = g, AUC = auc_val,
    CI_Lower = ci_val[1], CI_Upper = ci_val[3]))
  pdf(paste0("ROC_", g, "_GSE221521.pdf"), width = 6, height = 6)
  par(mfrow = c(1, 1))
  plot(r, col = "#2166ac", lwd = 2, legacy.axes = TRUE, main = g)
  legend("bottomright", legend = paste0("AUC = ", auc_val),
         bty = "n", text.col = "#c51b7d")
  dev.off()
}
write.csv(auc_summary[order(-auc_summary$AUC), ],
  "AUC_summary_GSE221521.csv", row.names = FALSE)

# PRS: extract weights from training CV model (FIXED: not hardcoded)
coefs_train <- coef(cv_fit$finalModel)
weights_ext <- coefs_train[-1]
# Align gene names
weights_ext <- weights_ext[found_genes_ext]
risk_score <- colSums(expr_final * weights_ext)

roc_prs <- roc(group_ext, risk_score, levels = c("Control", "DR"))
pdf("PRS_ROC_GSE221521.pdf", width = 5, height = 5)
plot(roc_prs, col = "#2166ac", lwd = 2, legacy.axes = TRUE,
  main = "PRS validation in GSE221521")
legend("bottomright",
  legend = paste0("AUC = ", round(auc(roc_prs), 3)),
  bty = "n", text.col = "#c51b7d")
dev.off()

cat("\n===== Bulk analysis complete =====\n")
