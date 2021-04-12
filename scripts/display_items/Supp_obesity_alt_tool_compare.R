rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(ggVennDiagram)

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/")
source("../../scripts/POMS_manuscript_functions.R")

almeida_DA_out <- readRDS(file = "Almeida_2019_DA_tool_output.rds")

almeida_POMS_out_ERP002061 <- readRDS("Almeida_2019_POMS_output/ERP002061_POMS_out.rds")

ERP002061_sig_KOs <- list()
ERP002061_sig_KOs[["POMS"]] <- rownames(almeida_POMS_out_ERP002061$df)[which(almeida_POMS_out_ERP002061$df$multinomial_corr < 0.25)]
ERP002061_sig_KOs[["ALDEx2"]] <- rownames(almeida_DA_out$ERP002061$aldex2)[which(almeida_DA_out$ERP002061$aldex2$BH_corr_p < 0.05)]
ERP002061_sig_KOs[["limma-voom"]] <- rownames(almeida_DA_out$ERP002061$limma.voom)[which(almeida_DA_out$ERP002061$limma.voom$BH_corr_p < 0.05)]
ERP002061_sig_KOs[["Wilcoxon test"]] <- rownames(almeida_DA_out$ERP002061$wilcoxon.relab)[which(almeida_DA_out$ERP002061$wilcoxon.relab$BH_corr_p < 0.05)]

ERP002061_venn <- ggVennDiagram(ERP002061_sig_KOs) +
                                scale_fill_gradient(name = "Count", low = "white", high = "firebrick3") +
                                ggtitle("Obesity 1 dataset sig. KOs") +
                                theme(plot.title = element_text(hjust = 0.5),
                                      plot.margin = unit(c(0, 0, 0, 0), "mm"))


almeida_POMS_out_ERP003612 <- readRDS("Almeida_2019_POMS_output/ERP003612_POMS_out.rds")

ERP003612_sig_KOs <- list()
ERP003612_sig_KOs[["POMS"]] <- rownames(almeida_POMS_out_ERP003612$df)[which(almeida_POMS_out_ERP003612$df$multinomial_corr < 0.25)]
ERP003612_sig_KOs[["ALDEx2"]] <- rownames(almeida_DA_out$ERP003612$aldex2)[which(almeida_DA_out$ERP003612$aldex2$BH_corr_p < 0.05)]
ERP003612_sig_KOs[["limma-voom"]] <- rownames(almeida_DA_out$ERP003612$limma.voom)[which(almeida_DA_out$ERP003612$limma.voom$BH_corr_p < 0.05)]
ERP003612_sig_KOs[["Wilcoxon test"]] <- rownames(almeida_DA_out$ERP003612$wilcoxon.relab)[which(almeida_DA_out$ERP003612$wilcoxon.relab$BH_corr_p < 0.05)]

ERP003612_venn <- ggVennDiagram(ERP003612_sig_KOs) +
                                scale_fill_gradient(name = "Count", low = "white", high = "firebrick3") +
                                ggtitle("Obesity 2 dataset sig. KOs") +
                                theme(plot.title = element_text(hjust = 0.5),
                                      plot.margin = unit(c(0, 0, 0, 0), "mm"))

almeida_venn_plot <- plot_grid(ERP002061_venn,
                               ERP003612_venn,
                               nrow = 1, ncol = 2, labels = c('a', 'b'))

ggsave(filename = "../../display_items/Supp_obesity_alt_tool_compare.pdf",
       plot = almeida_venn_plot,
       width = 10, height = 5, dpi = 600, device = "pdf")

ggsave(filename = "../../display_items/Supp_obesity_alt_tool_compare.png",
       plot = almeida_venn_plot,
       width = 10, height = 5, dpi = 600, device = "png")