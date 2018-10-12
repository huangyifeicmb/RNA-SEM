#!/usr/bin/env Rscript
library(lavaan)
args = commandArgs(trailingOnly=TRUE)

# build SEM
build_model <- function(var_names){
    # print(var_names)
    pro_var = grep("_PROseq$", var_names, perl=TRUE, value=TRUE)
    rna_var = grep("_Exon$", var_names, perl=TRUE, value=TRUE)
    covariate = setdiff(var_names, union(pro_var, rna_var))
    # print(var_names)
    # print(pro_var)
    # print(rna_var)
    # print(covariate)

    model_str <- "# Measurement model\n"
    for (x in pro_var){
        new_line = paste("latent_trans_ =~ 1 *", x)
        model_str = paste(model_str, new_line, sep="\n")
        # print(new_line)
        new_line = paste(x, "~~ residual_trans_ *", x)
        model_str = paste(model_str, new_line, sep="\n")
        new_line = paste(x, "~ 0*1\n")
        model_str = paste(model_str, new_line, sep="\n")
    }

    for (x in rna_var){
        new_line = paste("latent_RNA_ =~ 1 *", x)
        model_str = paste(model_str, new_line, sep="\n")
        # print(new_line)
        new_line = paste(x, "~~ residual_RNA_ *", x)
        model_str = paste(model_str, new_line, sep="\n")
        new_line = paste(x, "~ 0*1\n")
        model_str = paste(model_str, new_line, sep="\n")
    }

    model_str = paste(model_str, "# Regression model\n", sep="\n")

    new_line = paste("latent_trans_ ~", paste(covariate, collapse = ' + '))
    model_str = paste(model_str, new_line, sep="\n")
    
    new_line = paste("latent_RNA_ ~ 1 * latent_trans_ +", paste(covariate, collapse = ' + '))
    model_str = paste(model_str, new_line, sep="\n")

    model_str = paste(model_str, "latent_trans_ ~ 1", sep="\n")
    model_str = paste(model_str, "latent_RNA_ ~ 1", sep="\n")
    model_str = paste(model_str, "latent_trans_ ~~ latent_trans_", sep="\n")
    model_str = paste(model_str, "latent_RNA_ ~~ latent_RNA_", sep="\n")
    model_str = paste(model_str, "latent_RNA_ ~~ 0 * latent_trans_", sep="\n")

    return(model_str)
}



df <- read.csv(args[1], row.names=1, header=TRUE)
# print(head(df))
model_str <- build_model(names(df))
cat("Structural equation model:\n\n")
cat(model_str, "\n")

fit <- sem(model_str, data = df, verbose=FALSE)
cat("\n--------------------------------------\n")
cat("\nModel summary:\n\n")
print(summary(fit))
para <- parameterEstimates(fit)
para <- para[para$op == "~" & para$rhs != "latent_trans_", ]
para$response <- para$lhs
para$response[para$response == "latent_trans_"] <- "Transcription_rate"
para$response[para$response == "latent_RNA_"] <- "Half_life"
para$covariate <- para$rhs

cat("\n--------------------------------------\n")
cat("\nModel fitting:\n\n")
print(fitMeasures(fit))
write.csv(para[, c("response", "covariate", "est", "pvalue", "z", "se", "ci.lower", "ci.upper")], quote=FALSE, row.names=FALSE, file=args[2])
