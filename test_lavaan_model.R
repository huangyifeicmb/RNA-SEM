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

df <- df[is.finite(rowSums(df[, 2:5])), ]
# df <- df[complete.cases(df), ]
print(head(df))
model <- '
# obs model
RNA =~ 1 * H1_Exon + 1 * H2_Exon
trans =~ 1 * H1_PROseq + 1 * H2_PROseq
H1_Exon ~~ res_exon * H1_Exon
H2_Exon ~~ res_exon * H2_Exon
H1_PROseq ~~ res_pro * H1_PROseq
H2_PROseq ~~ res_pro * H2_PROseq
RNA ~ 1 * trans + U + G + C + Intron_number
trans ~ U + G + C + Intron_number

RNA ~ 1
trans ~ 1
H1_Exon ~ 0 * 1
H2_Exon ~ 0 * 1
H1_PROseq ~ 0 * 1
H2_PROseq ~ 0 * 1
'

# fit <- cfa(model, data = df, std.lv = TRUE)
fit <- sem(model, data = df)
print(summary(fit))
print(coef(fit))
# print(summary(fit, estimates = TRUE, fit.measures = TRUE))
