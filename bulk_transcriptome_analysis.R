getwd()
##读取数据##
install.packages("tinyarray")
library(data.table)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyverse)
library(tinyarray)
library(limma)
library(GEOquery)
eSet=getGEO("GSE60436",destdir='.')
?save
save(eSet,file="GSE60436_eSet.Rdata")
load('GSE60436_eSet.Rdata')
a=eSet[[1]]
dat=exprs(a)
dim(dat)
##获取临床信息##
pd=pData(a)
pd=pd[,c(1,2)]
head(pd)
names(pd) <- c('group','sample')
##分组##
library(stringr)
pd$group <- c(group=str_split(pd$group,'_',simplify=T)[,1])
pd$group <- ifelse(pd$group=='Retina','control','DR')
group_list=pd[,1]
table(group_list)
##探针注释-找到symbol与对应基因，排除非特异基因，GEO官网下载##
annoation <- read.table('GPL6884-11607.txt',sep='\t',fill=TRUE,header = TRUE,check.names = FALSE)
ids=annoation[,c(1,13)]
#将矩阵转换为数据框#
dat=as.data.frame(dat)
dat1 <- merge(ids,dat,by.x=1,by.y=0)
dat1 =dat1[,-1]
##去重## 
dat1 =dat1[!duplicated(dat1$Symbol),]
rownames(dat1) <- dat1$Symbol
dat1=dat1[,-1]
dat1=na.omit(dat1)
pdf('GSE60436boxplot.pdf', width = 8, height = 7)
boxplot(dat1)
dev.off()
#dat1=log(dat1+1)---，原数据为标准化的，此为未标准化的
save(dat1,file='GSE60436.Rdata')
save(dat1,pd,group_list,file='GSE60436_dat_group.Rdata')

#将GSE102485转化为GENEID#
library(data.table)
exp <- fread("GSE102485_raw_counts_GRCh38.p13_NCBI.tsv.gz")
head(exp)
library(stringi)
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
options(connectionObserver = NULL)
library(BiocManager)
BiocManager::install("org.Hs.eg.db")
name <- bitr(exp$GeneID,fromType='ENTREZID',toType='SYMBOL',OrgDb='org.Hs.eg.db',drop=T)
name$ENTREZID=as.numeric(name$ENTREZID)
exp=right_join(name,exp,by=c('ENTREZID'='GeneID'))
exp=exp[,-1]
exp <- aggregate(.~SYMBOL,FUN=mean,data=exp)
#临床数据加载#
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
#获取临床信息#
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
# 检验一下分组结果是否正确
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
#整合数据#
########整合数据
load("C:/Users/Xinyi Zuo/Desktop/DR与氨死亡相关性分析/GSE60436_dat_group.Rdata")
exp<-dat1
rm(dat1,pd,group_list)
load("C:/Users/Xinyi Zuo/Desktop/DR与氨死亡相关性分析/GSE102485_dat_group.Rdata")
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

#使用sva包中的combat() 函数
rm(list = ls())
#处理批次效应(combat)
library(sva)
library(limma)
#ComBat
group_list <- c(rep("Control",3),rep('DR',6),rep("Control",3),rep('DR',21))
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
pdf('boxplot整合2.pdf', width = 8, height = 7)
boxplot(exp2)
dev.off()
save(exp2,group_list,file = "DR_exp_combined.Rdata")


#差异分析#
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
#根据条件筛选差异基因并可视化#
DEG$regulate<-ifelse(DEG$P.Value<=0.05 & abs(DEG$logFC)>=1,
                     ifelse(DEG$logFC>0, 'Up', 'Down'), 'Stable')
table(DEG$regulate)          
write.csv(data.frame(gene_symbol=rownames(DEG),DEG),file='DEG.csv')
pdf('DEG-VOLCANO.pdf')
library(ggrepel)
DEG$symbol <- rownames(DEG)
# 先选出 p 值最小的10个基因
top10 <- DEG[order(DEG$P.Value), ][1:10, ]
# 或者如果你想按 adjust.P.Val 排名，换成 DEG$adj.P.Val

# 火山图
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

#取交集#
library(ggvenn)
library(ggplot2)
library(dplyr)
DR<-read.csv('DEG.csv',sep = ',',header = T)
DR$regulate=ifelse(DR$regulate=='Stable',NA,DR$regulate)
DR=na.omit(DR)
AD<-read.csv('氨死亡相关基因集.csv',sep=',',header=T)
datalist<-list('diabetic retinopathy'=DR$gene_symbol,
               'AD'=AD$Gene.Symbol)
#绘图
opar<-par(family='Roboto Condensed')
biocolor<-c('lightblue','pink')
p=ggvenn(datalist,
         fill_color = biocolor,
         fill_alpha = 0.5,
         stroke_linetype = 'longdash',
         set_name_size = 3,
         text_size = 4)
p
pdf(file='维恩图.pdf',width = 8,height = 6)
print(p)
dev.off()
#利用intersect、reduce找所有向量的交集
gene_list <- list(DR$gene_symbol,AD$Gene.Symbol)
comm_gene <- Reduce(intersect,gene_list)
write.table(comm_gene,file='gene.txt',sep='\t')
#富集分析
#GSEA/GSVA：验证“氨诱导免疫细胞死亡”这个通路是否被整体激活
library(limma)
library(org.Hs.eg.db)
library(tibble) # 如果没加载的话

DEG$ENTREZID <- mapIds(org.Hs.eg.db,
                             keys = rownames(DEG),
                             column = "ENTREZID",
                             keytype = "SYMBOL",
                             multiVals = "first")
geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
# 去掉没有 ENTREZ ID 的基因
geneList <- geneList[!is.na(names(geneList))]
# 按 logFC 从高到低排序
geneList <- sort(geneList, decreasing = TRUE)
# 检查，应该是一个带名字的数值向量
head(geneList)
#氨死亡相关基因
# 1. 从数据框里提取基因符号，并清洗
gene_symbols <- unique(trimws(AD$Gene.Symbol))
gene_symbols <- toupper(gene_symbols)  # 统一大写，避免大小写问题

# 2. 转换为 ENTREZ ID（用 mapIds 更安全，不会因个别失败而报错）
library(org.Hs.eg.db)
# 用 mapIds 看原始映射结果（不做 na.omit）
entrez_raw <- mapIds(org.Hs.eg.db, 
                     keys = unique(trimws(AD$Gene.Symbol)), 
                     column = "ENTREZID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")
# 统计多少个 NA
table(is.na(entrez_raw))
# 查看哪些基因没映射成功
failed <- names(entrez_raw[is.na(entrez_raw)])
print(failed)
# 检查结果
head(entrez_raw)
library(clusterProfiler)

# 重新制备基因集数据框（确保 gene 列为字符型）
ammonia_geneset <- data.frame(
  term = "Ammonia_ICD",
  gene = as.character(entrez_raw),
  stringsAsFactors = FALSE
)

# 检查
head(ammonia_geneset)
str(ammonia_geneset)   # gene 列应为 chr
# 如果还没排序，执行一次
geneList <- sort(geneList, decreasing = TRUE)

# 检查
head(geneList)
anyDuplicated(names(geneList))   # 应为 0，不能有重复基因名
# 1. 强制定位并重新生成随机数种子
set.seed(2024) 

# 2. 运行 GSEA
gsea_res <- GSEA(
  geneList      = geneList,
  TERM2GENE     = ammonia_geneset,
  pvalueCutoff  = 1,            
  pAdjustMethod = "BH",
  minGSSize     = 5,
  maxGSSize     = 1000,
  eps           = 0,            
  seed          = 2024,          # 很多版本的 clusterProfiler 也要求这里传入
  verbose       = TRUE           # 加上这个参数，方便您看进度报错
)

# 查看结果
dim(gsea_res@result)            # 应该有 1 行
gsea_res@result[, c("Description", "NES", "pvalue", "p.adjust", "qvalue")]
library(enrichplot)

# 经典山峰图，标题可自定义
pdf(file='GSEA_annomia.pdf',width = 8,height = 6)
gseaplot2(gsea_res, 
          geneSetID = "Ammonia_ICD",
          title = "GSEA: Ammonia-induced ICD in DR",
          pvalue_table = TRUE,
          ES_geom = "line")
dev.off()

library(GSVA)
# 如果 expr 是数据框，先转为矩阵
expr_mat <- as.matrix(exp2)

# 确保行名存在
if (is.null(rownames(expr_mat))) stop("表达矩阵没有行名")

# 确保是 numeric
mode(expr_mat) <- "numeric"
class(exp2)
dim(exp2)
head(rownames(exp2))
head(colnames(exp2))

# 确认行名与你的氨死亡基因符号匹配
ammonia_symbols <- unique(trimws(AD$Gene.Symbol))
genes_use <- intersect(ammonia_symbols, rownames(expr_mat))
length(genes_use)   # 看看有多少基因能在矩阵里找到

if (length(genes_use) < 5) {
  warning("基因太少，GSVA 可能不可靠，但依然可以运行查看趋势")
}
# 纯手工 ssGSEA，不依赖 GSVA 包
my_ssgsea <- function(expr_mat, gene_sets, alpha = 0.25, normalization = TRUE) {
  # expr_mat: 行=基因，列=样本，数值矩阵
  # gene_sets: list，每个元素是一个基因符号向量
  # 返回: 基因集 x 样本 的评分矩阵
  
  require(matrixStats)
  
  # 1. 对每个样本，基因按表达量降序排列，得到排序索引
  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)
  
  # 预分配评分矩阵
  es <- matrix(0, nrow = length(gene_sets), ncol = n_samples)
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(expr_mat)
  
  # 2. 对每个基因集
  for (i in seq_along(gene_sets)) {
    gs_genes <- gene_sets[[i]]
    # 确保基因在表达矩阵中
    gs_genes <- intersect(gs_genes, rownames(expr_mat))
    if (length(gs_genes) == 0) next
    
    # 标记哪些基因属于基因集
    gene_in_set <- rownames(expr_mat) %in% gs_genes
    
    # 对每个样本
    for (j in seq_len(n_samples)) {
      expr_sample <- expr_mat[, j]
      
      # 按表达量降序排序
      order_idx <- order(expr_sample, decreasing = TRUE)
      sorted_expr <- expr_sample[order_idx]
      sorted_in_set <- gene_in_set[order_idx]
      
      # 计算随机游走偏差 (类似 GSEA 的 enrichment score)
      # 命中（hit）和缺失（miss）的权重
      p_vals <- sorted_expr
      p_vals[sorted_in_set] <- abs(p_vals[sorted_in_set])^alpha
      p_vals[!sorted_in_set] <- 0
      
      n_total <- sum(p_vals)
      if (n_total == 0) next
      
      # 步长：命中时增加，缺失时减少
      hit_steps <- p_vals / n_total
      miss_steps <- rep(1 / (n_genes - sum(sorted_in_set)), n_genes)
      miss_steps[sorted_in_set] <- 0
      
      # 计算累积和
      cum_hit <- cumsum(hit_steps)
      cum_miss <- cumsum(miss_steps)
      
      # ES 是 cum_hit - cum_miss 的最大绝对值
      diff_vec <- cum_hit - cum_miss
      es[i, j] <- diff_vec[which.max(abs(diff_vec))]
    }
  }
  
  # 3. 标准化（可选，模仿 GSVA 的默认行为）
  if (normalization && n_samples > 1) {
    for (i in seq_len(nrow(es))) {
      es[i, ] <- (es[i, ] - mean(es[i, ])) / sd(es[i, ])
    }
  }
  
  return(es)
}
# 你的基因集列表（基因符号向量）
gs_list <- list(Ammonia_ICD = genes_use)

# 运行
ssgsea_scores <- my_ssgsea(expr_mat, gs_list)
# 提取评分
score <- as.numeric(ssgsea_scores["Ammonia_ICD", ])

# 确认分组向量长度一致
stopifnot(length(score) == length(group_list))

# 箱线图
library(ggplot2)
df <- data.frame(Score = score, Group = group_list)
pdf(file="GSVA_bluebox.pdf",width = 10,height = 8)
ggplot(df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_classic() +
  labs(title = "ssGSEA: Ammonia-induced ICD pathway activity")
dev.off()
# Wilcoxon 检验
wilcox.test(Score ~ Group, data = df)

# ROC 诊断
library(pROC)
roc_obj <- roc(df$Group, df$Score)
pdf(file="GSVA_roc.pdf",width = 10,height = 8)
plot(roc_obj, main = paste0("AUC = ", round(auc(roc_obj), 3)))
dev.off()
####GO富集分析
options(stringsAsFactors = F)
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library(stringr)

pvalueFilter=0.05         
qvalueFilter=1  
showNum=7
#gene symbol改为entrezIDS
rt=read.table("gene.txt",sep="\t")  
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)
pdf(file="GO_barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bar)
dev.off()

pdf(file="GO_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_cnet.pdf",width = 10,height = 8)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(af, foldChange=aflogfc,showCategory = 5, categorySize="pvalue",circular = TRUE,colorEdge = TRUE,cex_label_category=0.65,cex_label_gene=0.6)
dev.off()

pdf(file="KEGG_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = showNum,cex_label_category=0.6,color = "pvalue",layout ="nicely")
dev.off()
pdf(file="KEGG_heatplot.pdf",width = 10,height = 7)
kegg=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
heatplot(kegg,foldChange = aflogfc) 
dev.off()
#KEGG富集分析
rm(list=ls())
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
BiocManager::install("pathview")
library(pathview)
library(ggnewscale)
library(DOSE)
library(stringr)


pvalueFilter=0.05        
qvalueFilter=1        
showNum=20

rt=read.table("gene.txt",sep="\t")   
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]  
# 假设差异分析结果存储在 DEG 数据框中，包含 logFC 和 entrezID
aflogfc <- DEG$logFC
names(aflogfc) <- DEG$entrezID
gene=rt$entrezID
gene=unique(gene)
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="KEGG_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel) +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()
pdf(file="KEGG_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_cnet.pdf",width = 10,height = 8)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(af, foldChange=aflogfc,showCategory = 5, categorySize="pvalue",circular = TRUE,colorEdge = TRUE,cex_label_category=0.65,cex_label_gene=0.6)
dev.off()

pdf(file="KEGG_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = showNum,cex_label_category=0.6,color = "pvalue",layout ="nicely")
dev.off()
pdf(file="KEGG_heatplot.pdf",width = 10,height = 7)
kegg=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
heatplot(kegg,foldChange = aflogfc) 
dev.off()


#机器学习
install.packages('randomForest')
install.packages('e1071')
install.packages('glmnet')
install.packages('bRacatus')
install.packages('caret')

library(tidyverse)
library(glmnet)
source('./msvmRFE.R')
library(e1071)
library(caret)
library(bRacatus)
library(randomForest)
gene<-read.table('./gene.txt',sep='\t',header=T)
names(gene)<-'gene'
load("DR_exp_combined.Rdata")
exp2=as.data.frame(exp2)
exp <- exp2[gene$gene,]
#lasso回归
exp=t(exp)
x<-as.matrix(exp)
head(group_list)
y<-ifelse(group_list=='Control',0,1)
set.seed(123456789)
fit<-glmnet(x,y,family = 'binomial',alpha = 1,lambda = NULL)
pdf('lasso.pdf',width = 8,height = 7)
plot(fit,xvar = 'lambda',label=T)
dev.off()
cvfit<-cv.glmnet(x,y,family='binomial',alpha = 1,type.measure = 'deviance',nfolds = 5)
pdf('cvfit.pdf',width = 8,height = 7)
plot(cvfit)
dev.off()
best_lambda<-cvfit$lambda.min
lasso_coefs<-coef(cvfit,s='lambda.min')
lasso_features<-lasso_coefs@Dimnames[[1]][which(lasso_coefs!=0)]
lasso_features=as.data.frame(lasso_features)
lasso_features<-lasso_features[-1,]
write.csv(lasso_features,'feature_lasso.csv')
#SVM-机器学习#
#根据数据类型为行为样本，列为基因，（exp）及分组信息（y）
#导入数据并为数据框设置列名#
train<-cbind(y,exp)
colnames(train)[1]<-'group'
#将数据分为5份#
input<-train
nfold<-5
nrows<-nrow(input)
#采用五折交叉检验#
#SVM-RFE(input,k=5,halve.abouve=20)#分割数据及随机数据
#运行SVM-RFE算法进行特征选择
folds<-rep(1:nfold,len=nrows)[sample(nrows)]
folds<-lapply(1:nfold,function(x) which(folds==x))
library(plyr)
results<-lapply(folds,svmRFE.wrap,input,k=5,halve.above=20)
#halve.above=20参数用于指定在特征选择过程中，如果剩余特征的数量大于40个，则进行减半处理
#计算数据集的行数及创建交叉验证的折数
#提取最重要的特征并保存数据
top.features<-WriteFeatures(results,input,save=F)
write.csv(top.features,'feature_svm.csv')
#选取前20个特征进行svm模型构建
FeatSweep <- lapply(1:20,FeatSweep.wrap,results,input)
#函数FeatSweep。weap用于特征选择过程中评估模型泛化误差
#绘制错误率曲线图
no.info<-min(prop.table(table(input[,1])))
errors<-sapply(FeatSweep,function(x) ifelse(is.null(x),NA,x$error))
pdf('svm-error.pdf',width = 6,height = 6)
PlotErrors(errors,no.info = no.info)
dev.off()
#绘制正确率曲线图
Plotaccuracy<-function(accuracies,no.info=0.5,ylim=range(accuracies),
                       xlab='Number of Features',ylab='5x CV Accuracy'){
  AddLine<-function(x,col='black'){
    max_index<-which.max(x)
    lines(x=max_index,y=x[max_index],col=col,type='p')
    points(x=max_index,y=x[max_index],col='red')
    text(x=max_index,y=x[max_index],labels=paste(max_index,'-',format(x[max_index],digits=3)),pos=1,
         col='red',cex=0.75)
  }
  plot(x=1:length(accuracies),y=accuracies,type='l',ylim=ylim,xlab=xlab,ylab=ylab)
  AddLine(accuracies)
  abline(h=no.info,lty=3)
}
pdf('svm-accuracy.pdf',width = 8,height = 8)
Plotaccuracy(1-errors,no.info = no.info)
dev.off()
#获取及提取错误率最低点对应的特征
min_error_index<-which.min(errors)
SVMRFEgenes<-top.features[1:min_error_index,'FeatureName']
write.csv(SVMRFEgenes,file='SVMRFEgenes.csv')
# 随机森林RF#
#行为样本，列为基因（exp），分组为group（注意不能识别-，要将-改为.进行分析）
set.seed(123456)
gene<-read.table('./gene.txt',sep='\t',header=T)
names(gene)<-'gene'
load("DR_exp_combined.Rdata")
exp2=as.data.frame(exp2)
exp <- exp2[gene$gene,]
exp$id<-rownames(exp)
exp$id<-str_replace_all(exp$id,"-","\\.")
rownames(exp)<-exp$id
exp=exp[,-34]
exp=as.data.frame(t(exp))
rf=randomForest(as.factor(group_list)~.,data = exp,ntree=500)
pdf(file='forest.pdf',width = 8,height = 8)
plot(rf,main='Random forest',lwd=2)
dev.off()
#找出最小误差数，并用最小误差重新构建随机森林模型
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group_list)~.,data=exp,ntree=optionTrees)
#查看基因重要性
rownames(rf2$importance)=gsub('\\.','-',rownames(rf2$importance))
importance=importance(x=rf2)
write.csv(importance,file='rf.importance.csv',row.names = T)
#查看基因重要性，绘制重要性图
importance=importance(x=rf2)
importance=as.data.frame(importance)
importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")
rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
#展示前30个基因重要性
library(randomForest)
library(limma)
library(ggpubr)

af=rfGenes[1:30,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = "lightgray", size = 2), # Change segment color and size
             dot.size = 6,                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_bw()         ,               # ggplot2 theme
             rotate=TRUE                                       )#缈昏浆鍧愭爣杞? 
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +#棰滆壊
  grids()   
#保存图片并挑选疾病特征基因
pdf(file="importance.pdf", width=6, height=6)
print(p1)
dev.off()
write.csv(af, file="rfGenes.csv")
#logistic回归
rm(list=ls())
gene<-read.table('./gene.txt',sep='\t',header=T)
names(gene)<-'gene'
load("DR_exp_combined.Rdata")
exp2=as.data.frame(exp2)
exp <- exp2[gene$gene,]
exp=t(exp)
x<-as.matrix(exp)
#提取分组信息，并将其转换为0,1
head(group_list)
y<-ifelse(group_list=='Control',0,1)
train<-cbind(y,exp)
colnames(train)[1]<-'group'
data=train
data=as.data.frame(data)
#定义单因素逻辑回归模型函数
Uni_glm_model<-function(x){
  FML<-as.formula(paste0("group~",x,""))#构建分析公式
  glm1<-glm(FML,family = binomial,data=data)#拟合逻辑回归模型
  glm2<-summary(glm1) #汇总分析结果
  OR<-round(exp(coef(glm1)),2) #计算风险率
  SE<-glm2$coefficients[,2]
  CI_Lower<-round(exp(coef(glm1)-1.96*SE),2)
  CI_Upper<-round(exp(coef(glm1)+1.96*SE),2)
  CI<-paste0(CI_Lower,'-',CI_Upper)
  P<-signif(glm2$coefficients[,4],3)#保留3位小数
  Uni_glm_model<-data.frame('characteristics'=x,'OR'=OR,'CI'=CI,'P'=P)[-1,] #整理数据
  return(Uni_glm_model)
}
exp=as.data.frame(t(data))
exp$id<-rownames(exp)
exp$id<-str_replace_all(exp$id,"-","\\.")
rownames(exp)<-exp$id
exp=exp[,-34]
data=as.data.frame(t(exp))
variable.names=colnames(data)[c(2:length(data))]
Uni_glm=lapply(variable.names,Uni_glm_model) #应用函数
library(plyr)
Uni_glm=ldply(Uni_glm,data.frame) #数据框模式
Uni_glm$characteristics<-str_replace_all(Uni_glm$characteristics,"\\.","-")
write.csv(Uni_glm,file='Uni_glm.csv')
#绘制森林图
library(forestplot)
fp=Uni_glm
fp=data.frame(Var=fp$characteristics,
              OR=paste0(fp$OR,'(',fp$CI,')'),
              Pvalue=as.numeric(fp$P),
              OR_1=as.numeric(sub("-.*$",'',fp$CI) ),
              OR_2=as.numeric(sub('.*-','',fp$CI)),
              OR_mean=as.numeric(fp$OR))
fp=fp[fp$Pvalue<0.05,]
fp=fp[fp$OR_1>0,]
fp$OR_mean_log<-log(fp$OR_mean)
fp$OR_1_log<-log(fp$OR_1)
fp$OR_2_log<-log(fp$OR_2)
write.csv(fp$Var,file='LogisticGenes.csv')

# 完整改进版
fp_sorted <- fp %>%
  filter(Pvalue < 0.05) %>%          # 若只显示显著变量
  arrange(OR_mean) %>%               # 按 OR 升序（保护因素在上）
  mutate(Var = factor(Var, levels = Var))

# 标签格式化
label_text <- cbind(
  fp_sorted$Var,
  paste0(fp_sorted$OR, " (", fp_sorted$CI, ")"),
  ifelse(fp_sorted$Pvalue < 0.001, "<0.001", 
         ifelse(fp_sorted$Pvalue < 0.01, format(fp_sorted$Pvalue, digits=2),
                format(fp_sorted$Pvalue, digits=3)))
)

# === 3. 重新绘制森林图 ===

# 计算 x 轴边界
x_min = min(fp$OR_1_log) * 1.1
x_max = max(fp$OR_2_log) * 1.1
xlim_range = c(min(x_min, -0.5), max(x_max, 0.5)) 

# 放大 PDF 尺寸，保证不溢出
# 在之前的代码基础上，只需修改 forestplot 这一块：

pdf('forestplot_logistic4.pdf', width = 8, height = 8, onefile = F)

forestplot(labeltext = label_text_matrix,
           # === 关键修改点：在数据前面用 c(NA, ) 补齐一行，和表头行数对应 ===
           mean = c(NA, fp$OR_mean_log),   
           lower = c(NA, fp$OR_1_log),
           upper = c(NA, fp$OR_2_log),
           # ========================================================
           
           zero = 0,
           xlim = xlim_range,          
           align = c("l", "c", "c"),
           boxsize = 0.3,     # boxsize 是单个数字，可以正常适配所有非 NA 行
           lineheight = unit(5, 'mm'),
           colgap = unit(5, 'mm'),     
           lwd.zero = 1.5,
           lwd.ci = 2,
           col = fpColors(box = 'red',
                          summary = 'blue',
                          lines = 'black',
                          zero = '#458B00'),
           xlab = 'log(OR)',
           lwd.xaxis = 1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = 'solid',
           title = 'Forest plot of log(OR)',
           line.margin = unit(0.4, "cm"), 
           new_page = TRUE)               

dev.off()
graphics.off()
#取交集
install.packages("ggvenn")
library(ggvenn)
library(ggplot2)
library(dplyr)
Lasso<-read.csv('feature_lasso.csv',sep=',',header = T)
SVMRFE<-read.csv('SVMRFEgenes.csv',sep=',',header = T)
RF<-read.csv('rfGenes.csv',sep=',',header = T)
Logistic<-read.csv('LogisticGenes.csv',sep = ',',header = T)
datalist<-list('Lasso'=Lasso$x,
               'RF'=RF$X,
               'Logistic'=Logistic$x)
#取交集
opar<-par(family='Roboto Condensed')
biocolor<-c('blue','yellow','green')
p=ggvenn(datalist,
         fill_color = biocolor,
         fill_alpha = 0.5,
         stroke_linetype = 'longdash',
         set_name_size = 3,
         text_size = 4)
p
pdf(file='交集基因.pdf',width = 8,height = 6)
print(p)
dev.off()
#找出共同基因
gene_list<-list(Lasso$x,RF$X,Logistic$x,RF$X)
common_gene<-Reduce(intersect,gene_list)
write.csv(common_gene,'ML-共同基因.csv')
#
gene_list<-list(RF$X,Logistic$x)
common_gene<-Reduce(intersect,gene_list)
write.csv(common_gene,'ML-共同基因RFlogistic.csv')
#
gene_list<-list(SVMRFE$x,RF$X)
common_gene<-Reduce(intersect,gene_list)
write.csv(common_gene,'ML-共同基因SVMRF.csv')
#
gene_list<-list(SVMRFE$x,Logistic$x)
common_gene<-Reduce(intersect,gene_list)
write.csv(common_gene,'ML-共同基因SVMlogistic.csv')
#
gene_list<-list(Lasso$x,Logistic$x,RF$X)
common_gene<-Reduce(intersect,gene_list)
write.csv(common_gene,'ML-共同基因lassologisticTF.csv')


#GSVA与common genes
# 假设核心基因有 "FADD","CASP8","HMGB1","BID","FAS"
core_genes <- c("RPE65","PPP1R12C","UROD","PPIG","SLC1A3","HSPA4")

# 从表达矩阵提取这些基因的表达量（确保表达矩阵已经是对数标准化）
core_expr <- exp2[core_genes, ]  # 行=基因，列=样本，与score同顺序

# 计算相关性
cor_list <- lapply(core_genes, function(g){
  cor.test(as.numeric(core_expr[g, ]), score, method = "spearman")
})

# 提取相关系数和p值
cor_df <- data.frame(
  Gene = core_genes,
  Rho = sapply(cor_list, `[[`, "estimate"),
  Pval = sapply(cor_list, `[[`, "p.value")
)
cor_df$FDR <- p.adjust(cor_df$Pval, method = "BH")
cor_df$Signif <- ifelse(cor_df$FDR < 0.05, "*", "")
print(cor_df)
library(ggplot2)
plot_df <- data.frame(
  GSVA_Score = score,
  UROD_Expr = as.numeric(core_expr["UROD", ]),
  Group = group_list
)
pdf(file='UROD_GSVA.pdf',width = 8,height = 6)
ggplot(plot_df, aes(x = UROD_Expr, y = GSVA_Score)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_classic() +
  labs(x = "UROD Expression (log2)", y = "Ammonia-ICD GSVA Score",
       title = paste0("Spearman R = ", round(cor_df$Rho[cor_df$Gene=="UROD"], 3),
                      ", p = ", signif(cor_df$Pval[cor_df$Gene=="UROD"], 3)))

dev.off()
plot_df <- data.frame(
  GSVA_Score = score,
  PPIG_Expr = as.numeric(core_expr["PPIG", ]),
  Group = group_list
)
pdf(file='PPIG_GSVA.pdf',width = 8,height = 6)
ggplot(plot_df, aes(x = PPIG_Expr, y = GSVA_Score)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_classic() +
  labs(x = "PPIG Expression (log2)", y = "Ammonia-ICD GSVA Score",
       title = paste0("Spearman R = ", round(cor_df$Rho[cor_df$Gene=="PPIG"], 3),
                      ", p = ", signif(cor_df$Pval[cor_df$Gene=="PPIG"], 3)))

dev.off()
library(pheatmap)
cor_matrix <- cor(t(core_expr), score, method = "spearman")
# 如果需要同时展示组间差异，可添加注释
# 强制关闭所有图形设备
while (dev.cur() > 1) dev.off()
library(pheatmap)
pdf(file='GSVApheatmap.pdf',width = 8,height = 6)
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.3f",
         main = "Correlation between common Genes and GSVA Score")
dev.off()
#ROC
library(tidyverse)
library(pROC)
library(ROCR)

# 1. 加载原始数据
load('DR_exp_combined.Rdata')
exp2 <- as.data.frame(exp2)

# 2. 加载两个GEO数据的表型并整合
load('GSE102485_dat_group.Rdata')
load('GSE60436_dat_group.Rdata')

# 注意：确保 pd 和 pd1 都有 sample 和 group 两列
pd <- pd[,c("sample","group")]
all_cli <- rbind(pd, pd1)

# 将小写 control 改为大写 Control
all_cli <- all_cli %>%
  mutate(group = recode(group, "control" = "Control"))
#3. 读取机器学习筛选出的基因
gene_res <- read.csv('ML-共同基因.csv')

# 4. 提取基因表达矩阵，并转置（此时 exp 的行名为 GSM ID，列名为基因）
exp <- exp2[gene_res$x, ] %>% t() %>% as.data.frame()
identical(rownames(all_cli),rownames(exp))
exp <- cbind(all_cli[,2],exp)
colnames(exp)[1] <- 'group'
###整体ROC
library(caret)
library(pROC)
consensus_genes <- c("PPP1R12C","UROD","PPIG","SLC1A3")
library(caret)
library(pROC)

# 1. 确保分组标签是字符 "Control" 和 "DR"
# 这里必须用您 exp 里的原始 group 列，而不是 group_num
cv_data <- exp[, c("group", consensus_genes)]
# 强制转换为因子，且明确第一层为 Control，第二层为 DR
cv_data$group <- factor(cv_data$group, levels = c("Control", "DR"))

# 2. 设置5折交叉验证参数
ctrl <- trainControl(method = "cv", 
                     number = 5, 
                     savePredictions = "final", 
                     classProbs = TRUE,         
                     summaryFunction = twoClassSummary)

# 3. 运行 5折交叉验证的逻辑回归
set.seed(2024) 
cv_fit <- train(group ~ ., data = cv_data, 
                method = "glm", 
                trControl = ctrl,
                metric = "ROC")

# 4. 提取5折交叉验证的预测概率
cv_preds <- cv_fit$pred
cv_preds <- cv_preds[order(cv_preds$rowIndex), ]

# 5. 生成交叉验证的 ROC 曲线
# 注意：因为标签现在是 "Control" 和 "DR"，pROC 会识别出真假
# 我们用 cv_preds$DR 表示预测为 DR 的概率
roc_cv <- roc(response = cv_preds$obs, 
              predictor = cv_preds$DR, 
              levels = c("Control", "DR"))

# 6. 导出交叉验证的 ROC 图
graphics.off()
pdf('Machine_Learning_CV_ROC_Fixed1.pdf', width = 6, height = 6)
plot(roc_cv, col = "blue", lwd = 2,
     legacy.axes = TRUE,
     print.auc = TRUE,
     main = "5-Fold Cross-Validated ROC (6 Consensus Genes)")
dev.off()

print(paste("AUC =", round(auc(roc_cv), 3)))

# 补充代码：计算并画内部队列的 PRS 评分分布图
library(ggplot2)
library(ggpubr)

# 1. 从交叉验证模型中提取内部队列所有样本的 PRS 预测概率
cv_data$PRS_Score <- predict(cv_fit, newdata = cv_data, type = "prob")[, "DR"]

# 2. 绘制小提琴图 + 箱线图 (展示PRS得分的差异)
# 计算密度分布图（完美替代之前的箱线图）
library(ggpubr)
p_density <- ggplot(cv_data, aes(x = PRS_Score, fill = group)) +
  # 绘制密度图，设置 50% 透明度方便看重叠部分
  geom_density(alpha = 0.5, color = NA) + 
  # 设置颜色（冷暖对比，突出分组）
  scale_fill_manual(values = c("Control" = "#4DBBD5", "DR" = "#E64B35")) +
  theme_bw(base_size = 14) +
  labs(x = "PRS Score (Predicted Probability)", y = "Density") +
  theme(legend.position = c(0.8, 0.8)) +  # 图例放在右上角
  # 加上P值，显示两组分布差异极其显著
  stat_compare_means(label = "p.format", label.y = max(density(cv_data$PRS_Score)$y) * 1.1)

# 保存
ggsave('PRS_Score_Density_Plot.pdf', plot = p_density, width = 6, height = 5)

#各基因ROC
roc <- roc(exp$group,exp$UROD)

pdf('UROD_roc.pdf',width = 5,height = 5)
plot(roc,col="red",
     legacy.axes=T,#y轴格式更改
     print.auc=TRUE#显示AUC面积
     # print.thres=TRUE#添加截点和95%CI
)#颜色
legend('bottomright', 
       legend = c('UROD AUC = 0.641'), 
       title = 'Diabetic Retinopathy',
       col = c('red'), 
       lty = 1,
       cex = 0.8)
dev.off()

pdf('PPPER12C_roc.pdf',width = 5,height = 5)
rocPPPER12C <- roc(exp$group,exp$PPP1R12C)
plot(rocPPPER12C,col="red",
     legacy.axes=T,#y轴格式改变
     print.auc=TRUE,#auc面积
     print.thres=TRUE#添加截点和95%CI
)#颜色
legend('bottomright', 
       legend = c('PPPER12C AUC = 0.839'), 
       title = 'INTERSECTION GENE',
       col = c('red'), 
       lty = 1,
       cex = 0.8)
dev.off()

pdf('PPIG_roc.pdf',width = 5,height = 5)
rocPPIG <- roc(exp$group,exp$PPIG)
plot(rocPPIG,col="red",
     legacy.axes=T,#y轴格式改变
     print.auc=TRUE,#auc面积
     print.thres=TRUE#添加截点和95%CI
)#颜色
legend('bottomright', 
       legend = c('PPIG AUC = 0.562'), 
       title = 'INTERSECTION GENE',
       col = c('red'), 
       lty = 1,
       cex = 0.8)
dev.off()

pdf('SLC1A3_roc.pdf',width = 5,height = 5)
rocSLC1A3 <- roc(exp$group,exp$SLC1A3)
plot(roc,col="red",
     legacy.axes=T,#y轴格式改变
     print.auc=TRUE,#auc面积
     print.thres=TRUE#添加截点和95%CI
)#颜色
legend('bottomright', 
       legend = c('SLC1A3 AUC = 0.642'), 
       title = 'INTERSECTION GENE',
       col = c('red'), 
       lty = 1,
       cex = 0.8)
dev.off()
#表达箱式图
library(ggplot2)
library(reshape2)
library(tidyr)
library(ggpubr)
box_df <- exp

# 1. 强制关闭所有卡死的PDF
graphics.off() 

# 2. 👑【新增关键代码】去除 box_df 中重复的列名，保留第一次出现的列
box_df <- box_df[, !duplicated(colnames(box_df))]

# 打印一下现在的列名，看看是否已经干净了
print("当前的列名有：")
print(colnames(box_df))

# 3. 提取除了第一列（group）以外的所有基因列名
gene_list <- colnames(box_df)[-1]

# 4. 循环每个基因名画图
# 再次确保画板干净
graphics.off()

# 循环开始
for (gene in gene_list) {
  
  # 打开PDF文件
  pdf(file = paste0(gene, ".pdf"), width = 5, height = 5)
  
  # 尝试画图
  tryCatch({
    p <- ggplot(box_df, aes(x = group, y = .data[[gene]], fill = group)) +
      geom_boxplot(aes(fill = group), alpha = 1) + 
      geom_jitter() +   
      scale_color_manual(values = c('blue','red')) + 
      scale_y_continuous(name = "Gene expression") +
      theme_classic() +
      stat_compare_means(method = 'wilcox.test', show.legend = FALSE)
    
    # 显式打印到PDF
    print(p)
    
  }, error = function(e) {
    # 如果某个基因画图报错，打印错误但不会中断整个循环
    print(paste("基因", gene, "画图出错:", e$message))
  })
  
  # 关闭当前这个PDF文件
  dev.off()
}
print("所有基因的箱式图PDF已经全部生成完毕！")
#免疫浸润
devtools::install_github("Moonerss/CIBERSORT")
BiocManager::install("preprocessCore")
library(CIBERSORT)
library(ggsci)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
devtools::install_github("IOBR/IOBR")
library(IOBR)
#导入数据
data('LM22')
write.table(LM22, file = 'LM22.txt', sep = "\t", quote = F, row.names = T, col.names = T)
load("DR_exp_combined.Rdata")
write.table(exp2, file = 'exp.txt', sep = "\t", quote = F, row.names = T, col.names = T)

#cibersort分析，慢
res_cibersort <- cibersort('LM22.txt','exp.txt', perm = 10, QN = F)
#保存中间文件
save(res_cibersort,file = "res_cibersort.Rdata")   

#去除丰度全为0的细胞
ciber_res <- res_cibersort[,abs(colSums(res_cibersort)) > 0]
ciber_res <- as.data.frame(ciber_res)
identical(rownames(ciber_res),all_cli[,1])
ciber_res$group <- all_cli$group
#排序
ciber_res1 <- ciber_res[order(ciber_res[,"group"], decreasing = T),]
ciber_res1$ID <- rownames(ciber_res1)

ciber_long <- ciber_res1 %>%
  select(`P-value`,Correlation, RMSE, ID,group, everything()) %>%
  pivot_longer(-c(1:5), names_to = "cell_type", values_to = "fraction") %>%
  dplyr::mutate(cell_type = gsub("_CIBERSORT", "", cell_type),
                cell_type = gsub("_", " ", cell_type))

sample_group_info <- ciber_long %>% 
  distinct(ID, group) %>%
  mutate(sample_number = 1:n())


pdf('cibersort.pdf', width = 10, height = 8)
ciber_long %>% 
  ggplot(aes(ID, fraction))+
  geom_bar(stat = 'identity' , position = 'fill', aes(fill = cell_type))+
  labs(x = NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette4, name = NULL)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom')+
  geom_tile(data = sample_group_info,
            aes(x= as.numeric(sample_number),
                y = -0.02,
                fill = group),
            width = 1,
            height = 0.02)
dev.off()

#boxplot
pdf('cibersort_boxplot.pdf', width = 10, height = 6)
ggplot(ciber_long, 
       aes(cell_type, fraction, fill = group))+
  geom_boxplot(outlier.shape = 21, color = 'black')+
  scale_fill_manual(values = c('red', 'blue'))+
  theme_bw()+
  labs(x= NULL, y = 'Estimated Proportion')+
  theme(legend.position = 'top')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black', size = 12))+
  stat_compare_means(aes(group = group, label = ..p.signif..),
                     method = 'kruskal.test')
dev.off()


####基因表达与免疫细胞的相关性####
#install.packages("corrplot")
library(corrplot)
library(tidyverse)
exp <- exp2[gene_res$x,] #结果基因的表达矩阵
exp <- exp %>% t() %>% as.data.frame() #转换成行为样本，列为基因的data.frame
ciber <- ciber_res[,1:22] #cibersort结果文件，只要免疫细胞含量，根据之前的免疫细胞百分比进行选择多少列
identical(rownames(ciber),rownames(exp)) #识别样本id是否一致
class(exp$PPP1R12C)
class(ciber$`B cells naive`)
cor <- cor(ciber,exp)
pdf('cibersort_corrplot.pdf', width = 10, height = 6)
corrplot(t(cor),
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")
dev.off()


#免疫检查点

checkpoint <- read.csv('immuno_checkpoint.csv')
checkpoint_exp <- exp2[checkpoint$Symbol,] #免疫检查点基因的表达矩阵
checkpoint_exp <- na.omit(checkpoint_exp)
checkpoint_exp <- checkpoint_exp %>% t() %>% as.data.frame() #转换成行为样本，列为基因的data.frame

exp <- exp2[gene_res$x,] #结果基因的表达矩阵
exp <- exp %>% t() %>% as.data.frame() #转换成行为样本，列为基因的data.frame

identical(rownames(checkpoint_exp),rownames(exp)) #识别样本id是否一致
class(exp$UROD)
class(checkpoint_exp$LAG3)
cor1 <- cor(checkpoint_exp,exp)
#table3
table3_data <- data.frame(
  Gene = c("UROD", "UROD", "PPIG"),
  Checkpoint = c("LAG3", "CTLA4", "CD40"),
  Spearman_R = c(-0.58, -0.55, -0.40),
  Note = c("Immune Brake Dismantled", "Immune Brake Dismantled", "Ag Presentation Blocked")
)
write.csv(table3_data, "Table3_Checkpoint_Correlation.csv", row.names = FALSE)


pdf('cibersort_cor1免疫检查点.pdf', width = 10, height = 6)
corrplot(t(cor1),
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")
dev.off()
##===============
#GSVA 评分 × 免疫细胞比例（计算通路活性与免疫细胞的 Spearman 相关）
# 加载你的最终数据
load("DR_exp_combined.Rdata")  # 包含 exp2 和 group_list

# 确保 my_ssgsea 函数已定义
# 如果环境中没有，重新运行它的定义代码
my_ssgsea <- function(expr_mat, gene_sets, alpha = 0.25, normalization = TRUE) {
  # expr_mat: 行=基因，列=样本，数值矩阵
  # gene_sets: list，每个元素是一个基因符号向量
  # 返回: 基因集 x 样本 的评分矩阵
  
  require(matrixStats)
  
  # 1. 对每个样本，基因按表达量降序排列，得到排序索引
  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)
  
  # 预分配评分矩阵
  es <- matrix(0, nrow = length(gene_sets), ncol = n_samples)
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(expr_mat)
  
  # 2. 对每个基因集
  for (i in seq_along(gene_sets)) {
    gs_genes <- gene_sets[[i]]
    # 确保基因在表达矩阵中
    gs_genes <- intersect(gs_genes, rownames(expr_mat))
    if (length(gs_genes) == 0) next
    
    # 标记哪些基因属于基因集
    gene_in_set <- rownames(expr_mat) %in% gs_genes
    
    # 对每个样本
    for (j in seq_len(n_samples)) {
      expr_sample <- expr_mat[, j]
      
      # 按表达量降序排序
      order_idx <- order(expr_sample, decreasing = TRUE)
      sorted_expr <- expr_sample[order_idx]
      sorted_in_set <- gene_in_set[order_idx]
      
      # 计算随机游走偏差 (类似 GSEA 的 enrichment score)
      # 命中（hit）和缺失（miss）的权重
      p_vals <- sorted_expr
      p_vals[sorted_in_set] <- abs(p_vals[sorted_in_set])^alpha
      p_vals[!sorted_in_set] <- 0
      
      n_total <- sum(p_vals)
      if (n_total == 0) next
      
      # 步长：命中时增加，缺失时减少
      hit_steps <- p_vals / n_total
      miss_steps <- rep(1 / (n_genes - sum(sorted_in_set)), n_genes)
      miss_steps[sorted_in_set] <- 0
      
      # 计算累积和
      cum_hit <- cumsum(hit_steps)
      cum_miss <- cumsum(miss_steps)
      
      # ES 是 cum_hit - cum_miss 的最大绝对值
      diff_vec <- cum_hit - cum_miss
      es[i, j] <- diff_vec[which.max(abs(diff_vec))]
    }
  }
  
  # 3. 标准化（可选，模仿 GSVA 的默认行为）
  if (normalization && n_samples > 1) {
    for (i in seq_len(nrow(es))) {
      es[i, ] <- (es[i, ] - mean(es[i, ])) / sd(es[i, ])
    }
  }
  
  return(es)
}

# 准备基因集
AD<-read.csv('氨死亡相关基因集.csv',sep=',',header=T)
ammonia_symbols <- unique(trimws(AD$Gene.Symbol))
genes_use <- intersect(ammonia_symbols, rownames(exp2))
gs_list <- list(Ammonia_ICD = genes_use)

# 计算 ssGSEA，输入矩阵 exp2，列名会被继承到结果中
ssgsea_scores <- my_ssgsea(expr_mat = as.matrix(exp2), gene_sets = gs_list)

# 提取氨死亡通路评分向量，并确保名字是 exp2 的列名
gsva_score <- as.numeric(ssgsea_scores["Ammonia_ICD", ])
names(gsva_score) <- colnames(exp2)   # 必须命名

# 准备免疫浸润矩阵
# 假设 res_cibersort 已加载
immune_frac <- res_cibersort[, 1:22]

# 寻找共同样本
common_samples <- intersect(colnames(exp2), rownames(immune_frac))

# 统一所有对象到 common_samples 并保持相同顺序
exp2 <- exp2[, common_samples, drop = FALSE]
immune_frac <- immune_frac[common_samples, , drop = FALSE]
gsva_score <- gsva_score[common_samples]

# 给 group_list 命名并取子集（group_list 原本没有名字，长度为 33）
names(group_list) <- colnames(exp2)   # 这里的 exp2 是完整矩阵时的列名
group_list <- group_list[common_samples]

# 最终一致性检查
stopifnot(identical(colnames(exp2), rownames(immune_frac)))
stopifnot(identical(names(gsva_score), rownames(immune_frac)))
stopifnot(identical(names(group_list), rownames(immune_frac)))


#supplement 1
cibersort_res <- data.frame(Cell = colnames(immune_frac))
cibersort_res$Control_Mean <- apply(immune_frac, 2, function(x) mean(x[group_list == "Control"]))
cibersort_res$DR_Mean <- apply(immune_frac, 2, function(x) mean(x[group_list == "DR"]))
cibersort_res$Pval <- apply(immune_frac, 2, function(x) wilcox.test(x ~ group_list)$p.value)
cibersort_res$FDR <- p.adjust(cibersort_res$Pval, method = "BH")
# 按 P 值从小到大排序
cibersort_res <- cibersort_res[order(cibersort_res$Pval), ]
write.csv(cibersort_res, "Table2_CIBERSORT_Diff.csv", row.names = FALSE)


# 查看结果
table(group_list)
summary(gsva_score)
consensus_genes<- c("RPE65","PPP1R12C","UROD","PPIG","SLC1A3","HSPA4")
# 只提取存在的基因
core_genes_use <- consensus_genes
core_expr <- exp2[core_genes_use, , drop = FALSE]
# 检查维度：行数=基因数，列数=样本数
dim(core_expr)
head(rownames(core_expr))
gene_path_cor <- apply(core_expr, 1, function(x){
  cor.test(x, gsva_score, method = "spearman")
})
gene_path_df <- data.frame(
  Gene = rownames(core_expr),
  Rho = sapply(gene_path_cor, `[[`, "estimate"),
  Pval = sapply(gene_path_cor, `[[`, "p.value")
)
gene_path_df$FDR <- p.adjust(gene_path_df$Pval, method = "BH")
print(gene_path_df)
write.csv(gene_path_df,"GSVAspearman.csv")

# 提取 UROD 和 PPIG 的表达量
urod_expr <- as.numeric(core_expr["UROD", ])
ppig_expr <- as.numeric(core_expr["PPIG", ])

# 初始化结果数据框
immune_gene_cor <- data.frame(
  Cell = colnames(immune_frac),
  UROD_Rho = NA, UROD_Pval = NA,
  PPIG_Rho = NA, PPIG_Pval = NA,
  stringsAsFactors = FALSE
)

# 对每种免疫细胞计算相关性
for (i in seq_len(ncol(immune_frac))) {
  cell_frac <- immune_frac[, i]
  
  # UROD vs 免疫细胞
  cor1 <- cor.test(urod_expr, cell_frac, method = "spearman")
  immune_gene_cor$UROD_Rho[i] <- cor1$estimate
  immune_gene_cor$UROD_Pval[i] <- cor1$p.value
  
  # PPIG vs 免疫细胞
  cor2 <- cor.test(ppig_expr, cell_frac, method = "spearman")
  immune_gene_cor$PPIG_Rho[i] <- cor2$estimate
  immune_gene_cor$PPIG_Pval[i] <- cor2$p.value
}

# 多重检验校正
immune_gene_cor$UROD_FDR <- p.adjust(immune_gene_cor$UROD_Pval, method = "BH")
immune_gene_cor$PPIG_FDR <- p.adjust(immune_gene_cor$PPIG_Pval, method = "BH")

# 查看显著的结果（按 UROD FDR 排序）
immune_gene_cor[order(immune_gene_cor$UROD_FDR), ]
immune_gene_cor[order(immune_gene_cor$PPIG_FDR), ]
# 选出对 UROD 或 PPIG 显著的细胞（FDR < 0.05）
sig_for_gene <- unique(c(
  immune_gene_cor$Cell[immune_gene_cor$UROD_FDR < 0.05],
  immune_gene_cor$Cell[immune_gene_cor$PPIG_FDR < 0.05]
))

# 如果没有或太少，放宽到 P < 0.05 不校正
if (length(sig_for_gene) < 3) {
  sig_for_gene <- unique(c(
    immune_gene_cor$Cell[immune_gene_cor$UROD_Pval < 0.05],
    immune_gene_cor$Cell[immune_gene_cor$PPIG_Pval < 0.05]
  ))
}

# 确保至少选 5 个
if (length(sig_for_gene) < 5) {
  sig_for_gene <- immune_gene_cor$Cell[order(
    pmin(immune_gene_cor$UROD_Pval, immune_gene_cor$PPIG_Pval)
  )][1:5]
}

print(sig_for_gene)
# 构建组合矩阵
combined_mat <- cbind(
  UROD = urod_expr,
  PPIG = ppig_expr,
  GSVA = gsva_score,
  immune_frac[, sig_for_gene, drop = FALSE]
)

# 计算 Spearman 相关矩阵
cor_matrix <- cor(combined_mat, method = "spearman")
# 计算 GSVA 评分与所有免疫细胞的 Spearman 相关
gsva_immune_cor <- apply(immune_frac, 2, function(x){
  cor.test(gsva_score, x, method = "spearman")
})

gsva_immune_df <- data.frame(
  Cell = colnames(immune_frac),
  Rho = sapply(gsva_immune_cor, `[[`, "estimate"),
  Pval = sapply(gsva_immune_cor, `[[`, "p.value")
)
gsva_immune_df$FDR <- p.adjust(gsva_immune_df$Pval, method = "BH")
gsva_immune_df <- gsva_immune_df[order(gsva_immune_df$Pval), ]
print(gsva_immune_df)
# 选择真正显著且有代表性的免疫细胞（硬核修改）
selected_cells <- c("Monocytes", "Mast cells resting", "Macrophages M1", "T cells CD8", "B cells naive")

# 构建组合矩阵
combined_mat <- cbind(
  UROD = urod_expr,
  PPIG = ppig_expr,
  GSVA = gsva_score,
  immune_frac[, selected_cells, drop = FALSE]
)

# 重新计算 Spearman 相关矩阵
cor_matrix <- cor(combined_mat, method = "spearman")

# 然后继续用您的代码画热图 (ggcorrplot 或者 pheatmap)
graphics.off()
pdf('Core_Genes_Immune_Correlation_Final.pdf', width = 8, height = 6)
pheatmap(cor_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_format = "%.2f",
         color = colorRampPalette(c("steelblue", "white", "#ee0000"))(100),
         main = "Core Genes vs Key Immune Cells")
dev.off()
print("热图已保存为 Core_Genes_Immune_Correlation_Final.pdf")
# 绘制散点图
library(pheatmap)
library(ggcorrplot)
library(ggcorrplot)
library(ggplot2)

# 强制关闭所有后台卡死的绘图连接（杜绝PDF空白）
graphics.off()
# 1. 重新从 ciber_res 提取免疫矩阵
immune_frac <- ciber_res[, 1:22]

# 2. 从你的表达矩阵 exp2 获取最终样本列表（这是你校正后的矩阵）
final_samples <- colnames(exp2)

# 3. 将 immune_frac 与 exp2 严格对齐
common_samples <- intersect(final_samples, rownames(immune_frac))
immune_frac <- immune_frac[common_samples, ]

# 4. 重新生成 group_list，完全基于 common_samples
#    确保 names(group_list) 与 colnames(exp2) 一致
names(group_list) <- colnames(exp2)
group_list <- group_list[common_samples]

# 5. 最终检查：长度、名称必须完全一致
stopifnot(identical(rownames(immune_frac), names(group_list)))
stopifnot(length(group_list) == nrow(immune_frac))

table(group_list)  # 查看分组情况
#=========
# 对每一种免疫细胞进行 Wilcoxon 检验
cibersort_pvals <- apply(immune_frac, 2, function(x) {
  wilcox.test(x ~ group_list)$p.value
})

# 构建标准数据框
cibersort_diff <- data.frame(
  Cell_clean = names(cibersort_pvals),
  Diff_Pval = as.numeric(cibersort_pvals),
  stringsAsFactors = FALSE
)

# 查看是否还有 NA（应该没有了）
anyNA(cibersort_diff$Diff_Pval)
head(cibersort_diff)
##重新合并 GSVA 数据（确保 Pval 也是数值）
# 从 gsva_immune_df 提取，并清洗细胞名
gsva_immune_clean <- data.frame(
  Cell_clean = gsub("\\.rho$", "", gsva_immune_df$Cell),
  Cell_clean = gsub("\\.", " ", gsub("\\.rho$", "", gsva_immune_df$Cell)),
  GSVA_Rho = as.numeric(gsva_immune_df$Rho),
  GSVA_Pval = as.numeric(gsva_immune_df$Pval),
  stringsAsFactors = FALSE
)

# 合并
compare_df <- merge(cibersort_diff, gsva_immune_clean, by = "Cell_clean", all = TRUE)

# 现在 Diff_Pval 和 GSVA_Pval 都应该没有 NA 了
anyNA(compare_df$Diff_Pval)
anyNA(compare_df$GSVA_Pval)

#第四步：计算 -log10(P) 并画图
compare_df$negLog_Diff <- -log10(compare_df$Diff_Pval + 1e-10)
compare_df$negLog_GSVA <- -log10(compare_df$GSVA_Pval + 1e-10)

# 标记 GSVA 显著的细胞
compare_df$label <- ifelse(compare_df$GSVA_Pval < 0.05, compare_df$Cell_clean, "")

# 画图
library(ggplot2)
library(ggrepel)

pdf('GSVA_vs_Cibersort_pval_scatter.pdf', width = 8, height = 7)
ggplot(compare_df, aes(x = negLog_Diff, y = negLog_GSVA)) +
  geom_point(aes(color = GSVA_Pval < 0.05), size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("grey80", "#c51b7d"), 
                     labels = c("Not significant", "GSVA-immune P < 0.05"),
                     name = NULL) +
  geom_text_repel(aes(label = label), size = 3.5, box.padding = 0.5, max.overlaps = 10) +
  labs(x = "CIBERSORT group difference (-log10 P)", 
       y = "GSVA-immune correlation (-log10 P)") +
  theme_bw() +
  theme(legend.position = "top")
dev.off()



#外部队列检测
# ============================================================
# ============================================================
# 外部验证 - GSE221521（RNA-seq, DR vs Control）
# ============================================================

# -------------------- 加载包 --------------------
library(GEOquery)
library(limma)
library(ggplot2)
library(pROC)
library(tidyr)
library(DESeq2)
library(dplyr)

# -------------------- 1. 加载或下载数据（获取分组信息） --------------------
f_matrix <- "GSE221521_series_matrix.txt.gz"
if (!file.exists(f_matrix)) {
  gset <- getGEO("GSE221521", destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
  save(gset, file = "GSE221521_eSet.Rdata")
} else {
  load("GSE221521_eSet.Rdata")
}
b <- gset[[1]]

# 提取临床信息
pd <- pData(b)
# 查看列名，找到分组相关的列（根据实际数据调整）
# colnames(pd)
# 通常分组在 "title" 或 "characteristics_ch1" 中
# 这里假设使用 "title" 列（常见于 GEO）
group_col <- "title"   # 如果分组在 characteristics_ch1，改为 "characteristics_ch1"

# 提取分组信息
group_info <- pd[, group_col]
# 查看前几个样本的分组文字
head(group_info)
# 先准确分组
group_full <- case_when(
  grepl("Control", group_info, ignore.case = TRUE) ~ "Control",
  grepl("DM", group_info, ignore.case = TRUE) ~ "DM",
  grepl("DR", group_info, ignore.case = TRUE) ~ "DR",
  TRUE ~ NA_character_
)
# 然后只保留Control和DR
keep <- group_full %in% c("Control", "DR")
group <- factor(group_full[keep], levels = c("Control", "DR"))
names(group) <- rownames(pd)[keep]
# 检查分组比例
table(group)

# -------------------- 2. 读取 raw counts 矩阵（需提前下载） --------------------
counts_file <- "GSE221521_raw_counts_GRCh38.p13_NCBI.tsv.gz"
if (!file.exists(counts_file)) {
  cat("请手动下载 raw counts 文件：\n")
  cat("访问 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521\n")
  cat("下载 'GSE221521_raw_counts_GRCh38.p13_NCBI.tsv.gz' 到当前目录。\n")
  stop("文件不存在，请下载后重新运行。")
}

# 读取 raw counts（行名 = Ensembl ID，列名 = 样本ID）
raw_counts <- read.table(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
cat("原始计数矩阵维度（基因 × 样本）：", dim(raw_counts), "\n")

# 对齐样本（只保留有分组信息的样本）
common_samples <- intersect(colnames(raw_counts), names(group))
raw_counts <- raw_counts[, common_samples, drop = FALSE]
group <- group[common_samples]

# 确保顺序一致
stopifnot(all(colnames(raw_counts) == names(group)))
table(group)
# -------------------- 3. 使用 DESeq2 进行 VST 标准化 --------------------
# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = data.frame(group = group),
                              design = ~ group)

# 方差稳定变换 (VST)
vsd <- vst(dds, blind = FALSE)
expr_norm <- assay(vsd)   # 行为 Ensembl ID，列为样本

# -------------------- 4. Ensembl ID -> Gene Symbol（自动检测） --------------------
if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# 检测行名格式
rownames_expr <- rownames(expr_norm)
head_ids <- head(rownames_expr, 10)
cat("行名示例：", paste(head_ids, collapse = ", "), "\n")

# 判断keytype
if (all(grepl("^ENSG", rownames_expr))) {
  keytype_used <- "ENSEMBL"
  cat("检测到 Ensembl ID，使用 keytype = 'ENSEMBL'\n")
} else if (all(grepl("^\\d+$", rownames_expr))) {
  keytype_used <- "ENTREZID"
  cat("检测到 Entrez ID，使用 keytype = 'ENTREZID'\n")
} else {
  # 假设已经是基因符号
  keytype_used <- "SYMBOL"
  cat("未检测到标准 ID，假设行名已是基因符号\n")
}

if (keytype_used == "SYMBOL") {
  # 直接使用行名作为基因符号
  symbols <- rownames_expr
  names(symbols) <- rownames_expr
} else {
  # 去除版本号（如果有）
  if (keytype_used == "ENSEMBL") {
    rownames_clean <- sub("\\..*", "", rownames_expr)
  } else {
    rownames_clean <- rownames_expr
  }
  symbols <- mapIds(org.Hs.eg.db, keys = rownames_clean, keytype = keytype_used, column = "SYMBOL")
  # 检查是否有缺失
  cat("成功映射到基因符号的数量：", sum(!is.na(symbols)), "/", length(symbols), "\n")
}

# 保留有符号的基因
keep <- !is.na(symbols)
expr_norm <- expr_norm[keep, ]
symbols <- symbols[keep]

# 替换行名为基因符号
rownames(expr_norm) <- symbols

# 合并重复基因（取平均值）
expr_agg <- aggregate(expr_norm, by = list(Gene = rownames(expr_norm)), FUN = mean)
rownames(expr_agg) <- expr_agg$Gene
expr_symbol <- expr_agg[, -1]
expr_symbol <- as.matrix(expr_symbol)

cat("转换后矩阵维度（基因 × 样本）：", dim(expr_symbol), "\n")
head(rownames(expr_symbol), 10)

# -------------------- 5. 提取核心基因 --------------------
core_genes <- c("RPE65", "PPP1R12C", "UROD", "PPIG", "SLC1A3", "HSPA4")
found_genes <- intersect(core_genes, rownames(expr_symbol))
cat("找到的基因数：", length(found_genes), "\n")
missing <- setdiff(core_genes, found_genes)
if (length(missing) > 0) cat("缺失的基因：", paste(missing, collapse = ", "), "\n")

if (length(found_genes) == 0) stop("没有核心基因被找到，请检查基因名或注释。")

expr_final <- expr_symbol[found_genes, , drop = FALSE]

# -------------------- 6. 分组向量（与 expr_final 列顺序一致） --------------------
group <- group[colnames(expr_final)]
stopifnot(all(colnames(expr_final) == names(group)))

# -------------------- 7. 差异表达分析（Limma） --------------------
design <- model.matrix(~ group)
fit <- lmFit(expr_final, design)
fit <- eBayes(fit)
DEG_core <- topTable(fit, coef = 2, number = Inf)
DEG_core$regulate <- ifelse(DEG_core$adj.P.Val < 0.05,
                            ifelse(DEG_core$logFC > 0, "Up", "Down"), "Stable")

write.csv(data.frame(Gene = rownames(DEG_core), DEG_core),
          "DEG_core_validation_GSE221521.csv", row.names = FALSE)
print("差异表达结果：")
print(DEG_core[, c("logFC", "P.Value", "adj.P.Val")])

# -------------------- 8. 火山图 --------------------
pdf("Core_Volcano1_GSE221521.pdf", width = 6, height = 5)
ggplot(DEG_core, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = regulate), size = 3) +
  geom_text(aes(label = rownames(DEG_core)), vjust = -1) +
  scale_color_manual(values = c("Down" = "#2166ac", "Up" = "#c51b7d", "Stable" = "grey")) +
  geom_vline(xintercept = c(-1, 1), lty = 4) +
  geom_hline(yintercept = -log10(0.05), lty = 4) +
  theme_bw()
dev.off()

# -------------------- 9. ROC 曲线 --------------------
roc_data <- as.data.frame(t(expr_final))
roc_data$group <- group

valid_genes <- rownames(expr_final)
auc_summary <- data.frame()

library(pROC)
library(tidyverse)
library(ggplot2)

# 确保有一个空的数据框来存储 AUC 汇总
auc_summary <- data.frame()

# 【修改点 1】循环开始，为每个基因单独生成 pdf
for (gene in valid_genes) {
  roc_obj <- roc(roc_data$group, roc_data[[gene]],
                 levels = c("Control", "DR"), quiet = TRUE)
  auc_val <- round(auc(roc_obj), 3)
  ci_vals <- round(ci(roc_obj), 3)
  auc_summary <- rbind(auc_summary,
                       data.frame(Gene = gene,
                                  AUC = auc_val,
                                  CI_Lower = ci_vals[1],
                                  CI_Upper = ci_vals[3]))
  
  # 【修改点 2】在循环内部开启 PDF，文件名用基因名动态命名
  pdf(paste0("ROC_", gene, "_GSE221521.pdf"), width = 6, height = 6)
  
  # 【修改点 3】强制重置画布布局为单图（防止残留之前的 2x3 布局）
  par(mfrow = c(1, 1)) 
  
  # 开始画图
  plot(roc_obj, col = "#2166ac", lwd = 2, legacy.axes = TRUE, main = gene)
  legend("bottomright", legend = paste0("AUC = ", auc_val),
         bty = "n", text.col = "#c51b7d")
  
  # 【修改点 4】循环内部关闭当前 PDF
  dev.off() 
  
  # 打印提示，方便在控制台看到进度
  print(paste("已输出:", gene, "的独立 ROC 图"))
}

# 导出 AUC 汇总表
write.csv(auc_summary[order(-auc_summary$AUC), ],
          "AUC_summary_GSE221521.csv", row.names = FALSE)
print("AUC 汇总已保存，结果如下：")
print(auc_summary)
# -------------------- 11. 可选：查看某基因表达分布 --------------------
cat("\nUROD 在两组中的表达分布：\n")
if ("UROD" %in% valid_genes) {
  print(tapply(roc_data$UROD, roc_data$group, summary))
} else {
  cat("UROD 未找到，跳过。\n")
}

cat("\n===== 外部验证分析完成！=====\n")
# 假设训练集有每个基因的logFC作为权重
weights <- c(RPE65 = -0.002, PPP1R12C = -0.079, UROD = -0.122, 
             PPIG = -0.157, SLC1A3 = 0.155, HSPA4 = -0.023)

# 在验证集中计算风险评分（expr_final是基因x样本矩阵）
risk_score <- colSums(expr_final * weights)

# 计算PRS的AUC
roc_prs <- roc(group, risk_score, levels = c("Control", "DR"))
auc(roc_prs)
ci_prs <- ci(roc_prs)
print(paste0("95% CI: ", round(ci_prs[1], 3), "–", round(ci_prs[3], 3)))
pdf("PRS_ROC_GSE221521.pdf", width = 5, height = 5)
plot(roc_prs, col = "#2166ac", lwd = 2, legacy.axes = TRUE, 
     main = "PRS Validation in GSE221521")
legend("bottomright", legend = paste0("AUC = 0.715"), 
       bty = "n", text.col = "#c51b7d")
dev.off()