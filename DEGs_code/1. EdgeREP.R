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
group <- relevel(group, ref = "WT") # Set WT as the reference group
dge <- DGEList(counts = countdata_filtered, group = group)

# Normalize
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)

# Fit model and test for differential expression
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2) # Compare EP against WT (WT is the reference)

# Extract results
res <- topTags(lrt, n = Inf)$table
res$Geneid <- rownames(res)
res <- res %>%
  mutate(
    Group = case_when(
      (FDR < 0.05 & logFC > 1) ~ "UP",
      (FDR < 0.05 & logFC < -1) ~ "DOWN",
      TRUE ~ "STABLE"
    ),
    LogP = -log10(FDR)
  )

# Write results to file
write.csv(res, file = "resultEP/edgeR.csv", row.names = TRUE)

# Create volcano plot
p <- ggscatter(
  res, x = "logFC", y = "LogP",
  color = "Group", 
  palette = c('#2f5688', '#BBBBBB', "#CC0000"),
  size = 2, alpha = 0.3,
  title = "Volcano plot EP vs WT",
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
pdf("resultEP/volcanoEP_vs_WT.pdf", width = 8, height = 8)
print(p)
dev.off()
