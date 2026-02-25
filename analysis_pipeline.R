if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# 1. 加载包
library(DESeq2)

# 2. 读取数据 
data <- read.table("E:/GSE313754Samples.txt/GSE313754_FPKMs_allSamples.txt", header=TRUE, row.names=1)

"E:/GSE313754Samples.txt/GSE313754_FPKMs_allSamples.txt"
# 3. 设置分组信息 (根据 GEO 描述：4个对照组 vs 4个 MTHFD2 敲低组)
colData <- data.frame(
  row.names = colnames(data),
  condition = factor(c(rep("Control", 4), rep("Knockdown", 4)))
)

# 4. 运行 DESeq2
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# 5. 获取结果
res <- results(dds)
write.csv(res, "analysis_results.csv")

# 正式安装火山图工具
BiocManager::install("EnhancedVolcano")
# 安装并加载绘图包
if (!require("EnhancedVolcano")) BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# 绘制火山图
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'MTHFD2 Knockdown vs Control',
                pCutoff = 10e-5,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0)
ggsave("results/volcano_plot.png", width = 8, height = 8)
# 第一步：正式安装热图工具
install.packages(
  "pheatmap"
)

# 安装并加载热图包
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# 提取差异最显著的前50个基因
select <- order(res$padj)[1:50]
ntd <- normTransform(dds) # 转换数据，让颜色更好看
plot_data <- assay(ntd)[select, ]

# 绘制热图
pheatmap(plot_data, 
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         cluster_cols=TRUE, 
         annotation_col=colData,
         main = "Top 50 Differentially Expressed Genes")
