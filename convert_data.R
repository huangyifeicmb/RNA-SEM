library(lavaan)
df1 <- read.csv("CD4_Exon_PROseq_expression_20180925.csv")
df2 <- read.csv("Gene_features_20180925.csv")
# print(head(df1))
# print(head(df2))
df <- merge(df1, df2, by="X")
df$H1_Exon <- log2(df$H1_Exon)
df$H2_Exon <- log2(df$H2_Exon)
df$H1_PROseq <- log2(df$H1_PROseq)
df$H2_PROseq <- log2(df$H2_PROseq)

df <- df[complete.cases(df), ]
df <- df[is.finite(rowSums(df[, 2:5])), ]
print(head(df))
write.csv(df[, c(1:5, 10:13)], "test_data.csv", quote=FALSE, row.names=FALSE)
