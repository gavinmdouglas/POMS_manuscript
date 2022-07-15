rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/")

Almeida_2019_out <- readRDS("combined_output.rds")

ERP002061_out <-  Almeida_2019_out$ERP002061$kos

ERP003612_out <- Almeida_2019_out$ERP003612$kos


ERP002061_out_sig_df <- ERP002061_out$results[which(ERP002061_out$results$multinomial_corr < 0.25), ]
ERP002061_out_sig_df$func <- rownames(ERP002061_out_sig_df)
rownames(ERP002061_out_sig_df) <- NULL
ERP002061_out_sig_df$dataset <- "Obesity 1"

ERP003612_out_sig_df <- ERP003612_out$results[which(ERP003612_out$results$multinomial_corr < 0.25), ]
ERP003612_out_sig_df$func <- rownames(ERP003612_out_sig_df)
rownames(ERP003612_out_sig_df) <- NULL 
ERP003612_out_sig_df$dataset <- "Obesity 2"

combined_obesity_sig_df <- rbind(ERP002061_out_sig_df, ERP003612_out_sig_df)

combined_obesity_sig_df_fdr0.15 <- combined_obesity_sig_df[which(combined_obesity_sig_df$multinomial_corr < 0.15), ]

# Excluded many KOs above this cut-off:
nrow(combined_obesity_sig_df) - nrow(combined_obesity_sig_df_fdr0.15)

write.table(x = combined_obesity_sig_df_fdr0.15,
            file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_KOs_RAW.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
