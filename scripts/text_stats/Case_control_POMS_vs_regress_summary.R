rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/")

POMS_obesity_ko_sig <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_KOs_RAW.tsv",
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

POMS_obesity_higher_sig <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_higher_levels_RAW.tsv",
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")


regress_Almeida_out <- readRDS("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_regress_specificity_output/combined_output.rds")

regress_Almeida_out$ERP002061$kos$func <- rownames(regress_Almeida_out$ERP002061$kos)
regress_Almeida_out$ERP003612$kos$func <- rownames(regress_Almeida_out$ERP003612$kos)

regress_Almeida_out$ERP002061$pathways$func <- rownames(regress_Almeida_out$ERP002061$pathways)
regress_Almeida_out$ERP003612$pathways$func <- rownames(regress_Almeida_out$ERP003612$pathways)

regress_Almeida_out$ERP002061$modules$func <- rownames(regress_Almeida_out$ERP002061$modules)
regress_Almeida_out$ERP003612$modules$func <- rownames(regress_Almeida_out$ERP003612$modules)

regress_obesity_ko_results_ERP002061_sig <- regress_Almeida_out$ERP002061$kos[which(regress_Almeida_out$ERP002061$kos$BH < 0.25), ]
regress_obesity_ko_results_ERP003612_sig <- regress_Almeida_out$ERP003612$kos[which(regress_Almeida_out$ERP003612$kos$BH < 0.25), ]

POMS_obesity1_sig_kos <- POMS_obesity_ko_sig[which(POMS_obesity_ko_sig$dataset == "Obesity 1"), ]
length(which(POMS_obesity1_sig_kos$func %in% regress_obesity_ko_results_ERP002061_sig$func))

POMS_obesity2_sig_kos <- POMS_obesity_ko_sig[which(POMS_obesity_ko_sig$dataset == "Obesity 2"), ]
length(which(POMS_obesity2_sig_kos$func %in% regress_obesity_ko_results_ERP003612_sig$func))


regress_obesity_module_results_ERP002061_sig <- regress_Almeida_out$ERP002061$modules[which(regress_Almeida_out$ERP002061$modules$BH < 0.25), ]
regress_obesity_module_results_ERP003612_sig <- regress_Almeida_out$ERP003612$modules[which(regress_Almeida_out$ERP003612$modules$BH < 0.25), ]

POMS_obesity1_sig_modules <- POMS_obesity_higher_sig[which(POMS_obesity_higher_sig$dataset == "Obesity 1" & POMS_obesity_higher_sig$func_type == "Module"), ]
length(which(POMS_obesity1_sig_modules$func %in% regress_obesity_module_results_ERP002061_sig$func))

POMS_obesity1_sig_modules$func[which(POMS_obesity1_sig_modules$func %in% regress_obesity_module_results_ERP002061_sig$func)]


POMS_obesity2_sig_modules <- POMS_obesity_higher_sig[which(POMS_obesity_higher_sig$dataset == "Obesity 2" & POMS_obesity_higher_sig$func_type == "Module"), ]
length(which(POMS_obesity2_sig_modules$func %in% regress_obesity_module_results_ERP003612_sig$func))

POMS_obesity2_sig_modules$func[which(POMS_obesity2_sig_modules$func %in% regress_obesity_module_results_ERP003612_sig$func)]




regress_obesity_pathway_results_ERP002061_sig <- regress_Almeida_out$ERP002061$pathways[which(regress_Almeida_out$ERP002061$pathways$BH < 0.25), ]
regress_obesity_pathway_results_ERP003612_sig <- regress_Almeida_out$ERP003612$pathways[which(regress_Almeida_out$ERP003612$pathways$BH < 0.25), ]

POMS_obesity1_sig_pathways <- POMS_obesity_higher_sig[which(POMS_obesity_higher_sig$dataset == "Obesity 1" & POMS_obesity_higher_sig$func_type == "Pathway"), ]
length(which(POMS_obesity1_sig_pathways$func %in% regress_obesity_pathway_results_ERP002061_sig$func))

POMS_obesity2_sig_pathways <- POMS_obesity_higher_sig[which(POMS_obesity_higher_sig$dataset == "Obesity 2" & POMS_obesity_higher_sig$func_type == "Pathway"), ]
length(which(POMS_obesity2_sig_pathways$func %in% regress_obesity_pathway_results_ERP003612_sig$func))

