library(ggplot2)
library(dplyr)
library(DESeq2)
library(rlog)
library(pheatmap)

countdata <- read.csv("data/EPmerge_remove_afterPCA.csv", header = T, row.names = "Geneid", stringsAsFactors = F,check.names = F)
metadata<- read.csv("data/metaEP.csv", header = T, row.names = "Sample", stringsAsFactors = T,check.names = F)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~Condition)

# filter
filter <- rowSums(counts(dds) >= 20) >= 90
dds.filtered <- dds[filter, ]
nrow(counts(dds))-nrow(counts(dds.filtered))
vst.deseq <- vst(dds.filtered, blind = FALSE) #transformed daa
vst1 <- as.data.frame(assays(vst.deseq))
vst1$group <- NULL
vst1$group_name <- NULL
write.csv(vst1, file = "resultEP/norm_EP.csv", row.names = T)

#Get the result
dds <- DESeq(dds.filtered)
resultsNames(dds)
#Get the result 
res <- results(dds,tidy=TRUE, contrast=c("Condition","EP","WT"))
head(res)
res <- tbl_df(res)
res %>% arrange(padj) %>% head()
head(as.data.frame(res))
write.csv(res, file = "resultEP/DESEQEP.csv", row.names = T)
# remove NA
#res <- read.csv("C:/Users/InkMu/Project/DESEQADR.csv", header = T, row.names = "row", stringsAsFactors = F,check.names = F)
res$Group = "STABLE"
#logFC_t <- with(res,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange))) 
logFC_t <- 1
res$Group[which((res$padj < 0.05) & (res$log2FoldChange > logFC_t)) ] = "UP"
res$Group[which((res$padj < 0.05) & (res$log2FoldChange < -logFC_t)) ] = "DOWN"
table(res$Group)
res$LogP <- -log10(res$padj)
res %>% 
  filter(padj<0.05 & abs(log2FoldChange) > logFC_t ) %>% 
  write.csv("resultEP/sigresultsEP.csv")
res <- data.frame(res)
res$Group <- as.factor(res$Group)
library(ggpubr)
library(ggthemes)
library(ggrepel)
p <- ggscatter(res, x="log2FoldChange", y="LogP", 
               color = "Group", 
               palette = c('#2f5688','#BBBBBB',"#CC0000"), 
               size = 2,
               alpha = 0.3,  # 设置点的透明度
               title = "Volcano plot EP",
               font.label = 10,
               repel = T,
               xlab = "log2FoldChange",
               ylab = "-log10(Adjust P-value)") + 
  theme_base() + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed") + 
  geom_vline(xintercept = c(-logFC_t,logFC_t), linetype="dashed") +
  xlim(-3.8, 3.7)

p <- p + scale_color_manual(values=c('#2f5688','#BBBBBB',"#CC0000"), 
                            labels = c(paste("Downregulated genes:", table(res$Group)[1]), 
                                       paste("No change:",table(res$Group)[2]),
                                       paste("Upregulated genes:", table(res$Group)[3])))
p <- p + theme(legend.position=c(0.2, 0.87),legend.text = element_text(size = 7), legend.title=element_blank())
print(p)
# 打印图形
pdf("resultEP/volcanoEP.pdf", width=8, height=8)
print(p)
dev.off()


