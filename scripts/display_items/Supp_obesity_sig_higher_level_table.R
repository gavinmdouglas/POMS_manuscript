rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/")

ERP002061_out_pathways <- readRDS(file = "ERP002061_POMS_out.rds")[["pathways"]]

ERP003612_out_pathways <- readRDS(file = "ERP003612_POMS_out.rds")[["pathways"]]


ERP002061_out_pathways_sig_df <- ERP002061_out_pathways$df[which(ERP002061_out_pathways$df$multinomial_corr < 0.25), ]
ERP002061_out_pathways_sig_df$func <- rownames(ERP002061_out_pathways_sig_df)
rownames(ERP002061_out_pathways_sig_df) <- NULL
ERP002061_out_pathways_sig_df$dataset <- "Obesity 1"
ERP002061_out_pathways_sig_df$func_type <- "Pathway"

ERP003612_out_pathways_sig_df <- ERP003612_out_pathways$df[which(ERP003612_out_pathways$df$multinomial_corr < 0.25), ]
ERP003612_out_pathways_sig_df$func <- rownames(ERP003612_out_pathways_sig_df)
rownames(ERP003612_out_pathways_sig_df) <- NULL 
ERP003612_out_pathways_sig_df$dataset <- "Obesity 2"
ERP003612_out_pathways_sig_df$func_type <- "Pathway"


ERP002061_out_modules <- readRDS(file = "ERP002061_POMS_out.rds")[["modules"]]

ERP003612_out_modules <- readRDS(file = "ERP003612_POMS_out.rds")[["modules"]]


ERP002061_out_modules_sig_df <- ERP002061_out_modules$df[which(ERP002061_out_modules$df$multinomial_corr < 0.25), ]
ERP002061_out_modules_sig_df$func <- rownames(ERP002061_out_modules_sig_df)
rownames(ERP002061_out_modules_sig_df) <- NULL
ERP002061_out_modules_sig_df$dataset <- "Obesity 1"
ERP002061_out_modules_sig_df$func_type <- "Module"

ERP003612_out_modules_sig_df <- ERP003612_out_modules$df[which(ERP003612_out_modules$df$multinomial_corr < 0.25), ]
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
