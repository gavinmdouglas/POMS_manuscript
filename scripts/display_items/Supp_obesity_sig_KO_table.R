rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/")

ERP002061_out <- readRDS(file = "ERP002061_POMS_out.rds")[["ko"]]

ERP003612_out <- readRDS(file = "ERP003612_POMS_out.rds")[["ko"]]


ERP002061_out_sig_df <- ERP002061_out$df[which(ERP002061_out$df$multinomial_corr < 0.25), ]
ERP002061_out_sig_df$func <- rownames(ERP002061_out_sig_df)
rownames(ERP002061_out_sig_df) <- NULL
ERP002061_out_sig_df$dataset <- "Obesity 1"

ERP003612_out_sig_df <- ERP003612_out$df[which(ERP003612_out$df$multinomial_corr < 0.25), ]
ERP003612_out_sig_df$func <- rownames(ERP003612_out_sig_df)
rownames(ERP003612_out_sig_df) <- NULL 
ERP003612_out_sig_df$dataset <- "Obesity 2"

combined_obesity_sig_df <- rbind(ERP002061_out_sig_df, ERP003612_out_sig_df)

combined_obesity_sig_df_fdr0.15 <- combined_obesity_sig_df[which(combined_obesity_sig_df$multinomial_corr < 0.15), ]

# Excluded many KOs under lenient cut-off:
nrow(combined_obesity_sig_df) - nrow(combined_obesity_sig_df_fdr0.15)

write.table(x = combined_obesity_sig_df,
            file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_KOs_RAW.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
