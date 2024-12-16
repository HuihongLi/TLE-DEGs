library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggthemes)
library(ggrepel)

# Load data
countdata <- read.csv("data/EPmerge_remove_afterPCA.csv", header = T, row.names = "Geneid", stringsAsFactors = F, check.names = F)
metadata <- read.csv("data/metaEP.csv", header = T, row.names = "Sample", stringsAsFactors = T, check.names = F)

# Manual filtering to match DESeq2's standard
filter <- rowSums(countdata >= 20) >= 90
countdata_filtered <- countdata[filter, ]
nrow(countdata) - nrow(countdata_filtered)

# Create DGEList object
group <- factor(metadata$Condition)
group <- relevel(group, ref = "WT")
dge <- DGEList(counts = countdata_filtered, group = group)

# Normalize
dge <- calcNormFactors(dge)

# Voom transformation
design <- model.matrix(~ group) # Automatically treats the first level as the reference (WT)
v <- voom(dge, design, plot = TRUE)

# Fit the linear model
fit <- lmFit(v, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract results for EP vs WT (WT is the reference)
res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
res$Geneid <- rownames(res)
res <- res %>%
  mutate(
    Group = case_when(
      (adj.P.Val < 0.05 & logFC > 1) ~ "UP",
      (adj.P.Val < 0.05 & logFC < -1) ~ "DOWN",
      TRUE ~ "STABLE"
    ),
    LogP = -log10(adj.P.Val)
  )

# Write results to file
write.csv(res, file = "resultEP/limma.csv", row.names = TRUE)

# Create volcano plot
p <- ggscatter(
  res, x = "logFC", y = "LogP",
  color = "Group", 
  palette = c('#2f5688', '#BBBBBB', "#CC0000"),
  size = 2, alpha = 0.3,
  title = "Volcano plot EP vs WT (limma)",
  xlab = "log2FoldChange",
  ylab = "-log10(Adjust P-value)"
) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlim(-3.8, 3.7) +
  scale_color_manual(
    values = c('#2f5688', '#BBBBBB', "#CC0000"),
    labels = c(
      paste("Downregulated genes:", sum(res$Group == "DOWN")),
      paste("No change:", sum(res$Group == "STABLE")),
      paste("Upregulated genes:", sum(res$Group == "UP"))
    )
  ) +
  theme(
    legend.position = c(0.88, 0.9),
    legend.text = element_text(size = 7),
    legend.title = element_blank()
  )

print(p)

# Save plot
pdf("resultEP/volcanoEP_vs_WT_limma.pdf", width = 8, height = 8)
print(p)
dev.off()
