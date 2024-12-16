# Load necessary packages
library(edgeR)
library(ggplot2)
library(dplyr)

# Set file paths
# Load data
countdata <- read.csv("data/EPmerge_remove_afterPCA.csv", header = T, row.names = "Geneid", stringsAsFactors = F, check.names = F)
metadata <- read.csv("data/metaEP.csv", header = T, row.names = "Sample", stringsAsFactors = T, check.names = F)

# Manual filtering to match DESeq2's standard
filter <- rowSums(countdata >= 20) >= 90
countdata_filtered <- countdata[filter, ]
nrow(countdata) - nrow(countdata_filtered)

conditions <- factor(metadata$Condition)
conditions <- relevel(conditions, ref = "WT") # Set WT as the reference group

# TMM normalization using edgeR
y <- DGEList(counts = countdata_filtered, group = conditions)

# Perform TMM normalization and convert to CPM (Counts Per Million)
y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)

# Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm), function(i) {
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i, ])), conditions)
  p <- wilcox.test(gene ~ conditions, data)$p.value
  return(p)
})
fdr <- p.adjust(pvalues, method = "BH")

# Calculate log2 fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[, c(which(conditions == conditionsLevel[1]))]
dataCon2 <- count_norm[, c(which(conditions == conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2) / rowMeans(dataCon1))

# Output results based on FDR threshold
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
write.csv(outRst, file = "resultEP/Wilcoxon.csv", row.names = TRUE)

# Prepare data for plotting
outRst$Group <- "STABLE"
outRst$Group[outRst$FDR < fdrThres & outRst$log2foldChange > 1] <- "UP"
outRst$Group[outRst$FDR < fdrThres & outRst$log2foldChange < -1] <- "DOWN"
outRst$LogP <- -log10(outRst$FDR)

# Create volcano plot
p <- ggscatter(
  outRst, x = "log2foldChange", y = "LogP",
  color = "Group", 
  palette = c('#2f5688', '#BBBBBB', "#CC0000"), # Down, Stable, Up
  size = 2, alpha = 0.3,
  title = "Volcano plot EP vs WT (Wilcoxon Test)",
  xlab = "log2 Fold Change",
  ylab = "-log10(FDR)"
) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlim(-3.8, 3.7) +
  scale_color_manual(
    values = c('#2f5688', '#BBBBBB', "#CC0000"),
    labels = c(
      paste("Downregulated genes:", sum(outRst$Group == "DOWN")),
      paste("No change:", sum(outRst$Group == "STABLE")),
      paste("Upregulated genes:", sum(outRst$Group == "UP"))
    )
  ) +
  theme(
    legend.position = c(0.88, 0.9),
    legend.text = element_text(size = 7),
    legend.title = element_blank()
  )

print(p)

# Save plot
pdf("resultEP/volcanoEP_vs_WT_Wilcoxon.pdf", width = 8, height = 8)
print(p)
dev.off()