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

ERP002061_out_pathways <- Almeida_2019_out$ERP002061$pathways

ERP003612_out_pathways <- Almeida_2019_out$ERP003612$pathways


ERP002061_out_pathways_sig_df <- ERP002061_out_pathways$results[which(ERP002061_out_pathways$results$multinomial_corr < 0.25), ]
ERP002061_out_pathways_sig_df$func <- rownames(ERP002061_out_pathways_sig_df)
rownames(ERP002061_out_pathways_sig_df) <- NULL
ERP002061_out_pathways_sig_df$dataset <- "Obesity 1"
ERP002061_out_pathways_sig_df$func_type <- "Pathway"

ERP003612_out_pathways_sig_df <- ERP003612_out_pathways$results[which(ERP003612_out_pathways$results$multinomial_corr < 0.25), ]
ERP003612_out_pathways_sig_df$func <- rownames(ERP003612_out_pathways_sig_df)
rownames(ERP003612_out_pathways_sig_df) <- NULL 
ERP003612_out_pathways_sig_df$dataset <- "Obesity 2"
ERP003612_out_pathways_sig_df$func_type <- "Pathway"


ERP002061_out_modules <- Almeida_2019_out$ERP002061$modules

ERP003612_out_modules <- Almeida_2019_out$ERP003612$modules


ERP002061_out_modules_sig_df <- ERP002061_out_modules$results[which(ERP002061_out_modules$results$multinomial_corr < 0.25), ]
ERP002061_out_modules_sig_df$func <- rownames(ERP002061_out_modules_sig_df)
rownames(ERP002061_out_modules_sig_df) <- NULL
ERP002061_out_modules_sig_df$dataset <- "Obesity 1"
ERP002061_out_modules_sig_df$func_type <- "Module"

ERP003612_out_modules_sig_df <- ERP003612_out_modules$results[which(ERP003612_out_modules$results$multinomial_corr < 0.25), ]
ERP003612_out_modules_sig_df$func <- rownames(ERP003612_out_modules_sig_df)
rownames(ERP003612_out_modules_sig_df) <- NULL 
ERP003612_out_modules_sig_df$dataset <- "Obesity 2"
ERP003612_out_modules_sig_df$func_type <- "Module"

combined_obesity_sig_df <- rbind(ERP002061_out_pathways_sig_df, ERP003612_out_pathways_sig_df,
                                 ERP002061_out_modules_sig_df, ERP003612_out_modules_sig_df)

write.table(x = combined_obesity_sig_df,
            file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_higher_levels_RAW.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
