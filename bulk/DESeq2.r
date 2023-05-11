#
mat <- mat[, c(4:6,1:3)]

# build condition
condition <- factor(c(rep("PBS",3),rep("DA",3)), levels = c("PBS","DA"))
colData <- as.data.frame(condition)
rownames(colData) <- colnames(mat) 

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = colData,
                              design = ~ condition ,
                              tidy = FALSE)

dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

res <- res[order(res$padj),]
head(res)

#write.table(as.data.frame(res), "DEG_DESeq2.csv", col.names = T, row.names = T, quote = FALSE, sep = ",")


