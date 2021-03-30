rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(ggVennDiagram)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")
source("../scripts/POMS_manuscript_functions.R")

DA_tool_BY_sig_features <- function(DA_tool_out, correction="BY") {
  DA_tool_BY <- list()
  DA_tool_BY[["ALDEx2"]] <- rownames(DA_tool_out$aldex2)[which(p.adjust(DA_tool_out$aldex2$wi.ep, correction) < 0.05)]
  DA_tool_BY[["DESeq2"]] <- rownames(DA_tool_out$deseq2)[which(p.adjust(DA_tool_out$deseq2$pvalue, correction) < 0.05)]
  DA_tool_BY[["Limma-Voom"]] <- rownames(DA_tool_out$limma.voom)[which(p.adjust(DA_tool_out$limma.voom$P.Value, correction) < 0.05)]
  DA_tool_BY[["Wilcoxon (relab.)"]] <- rownames(DA_tool_out$wilcoxon.relab)[which(p.adjust(DA_tool_out$wilcoxon.relab$wilcox_p, correction) < 0.05)]
  DA_tool_BY[["Wilcoxon (USCG-norm.)"]] <- rownames(DA_tool_out$wilcoxon.musicc)[which(p.adjust(DA_tool_out$wilcoxon.musicc$wilcox_p, correction) < 0.05)]
  return(DA_tool_BY)
}

DA_out <- readRDS(file = "ref.based_ex_KOs_DA_out.rds")

DA_out <- DA_out[-grep("wilcoxon.musicc", names(DA_out))]

settings <- c("mu0.1_nu_0.5\"", "mu0.1_nu_0.9\"", "mu0.1_nu_0.99\"", "mu0.01_nu_0.99\"")
representative_KOs <- c("K00480", "K01845", "K06077")
tools <- c("aldex2", "limma.voom", "wilcoxon.relab")

DA_split_out <- list()

for(setting in settings) {
 
  matching_ids <- grep(setting, names(DA_out), value=TRUE)
  
  for(ko in representative_KOs) {

    ko_matching_ids <- grep(ko, matching_ids, value=TRUE)
    
    setting_clean <- gsub("\"", "", setting)
    DA_split_out[[paste(setting_clean, ko, sep=" ")]] <- list()
    
    for(tool in tools) {
      
        tool_ko_matching_id <- grep(tool, ko_matching_ids, value=TRUE)
      
        DA_split_out[[paste(setting_clean, ko, sep=" ")]][[tool]] <- DA_out[[tool_ko_matching_id]]
        
    }
    
  }
  
}

DA_out_features_BY_0.05 <- lapply(DA_split_out, DA_tool_BY_sig_features)

POMS_out <- readRDS("ref.genome_sim_summaries/ref.genome_sim_ex_KOs.rds")


representative_KOs <- c("K00480", "K01845", "K06077")



for(setting in settings) {
  for(ko in representative_KOs) {
    category <- paste(setting, ko)
    
    DA_out_features_BY_0.05[[category]][["POMS"]] <- names(which(p.adjust(POMS_out[[setting]][[ko]]$pseudo_null_p, "BY") < 0.05))
    
  }
}

names(DA_out_features_BY_0.05) <- gsub("mu0.1_nu_0.5 ", "Setting 1, ", names(DA_out_features_BY_0.05))
names(DA_out_features_BY_0.05) <- gsub("mu0.1_nu_0.9 ", "Setting 2, ", names(DA_out_features_BY_0.05))
names(DA_out_features_BY_0.05) <- gsub("mu0.1_nu_0.99 ", "Setting 3, ", names(DA_out_features_BY_0.05))
names(DA_out_features_BY_0.05) <- gsub("mu0.01_nu_0.99 ", "Setting 4, ", names(DA_out_features_BY_0.05))

venn_plots <- list()

for(category in names(DA_out_features_BY_0.05)) {
  venn_plots[[category]] <- ggVennDiagram(DA_out_features_BY_0.05[[category]]) +
                                            scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
                                            ggtitle(category) +
                                            theme(plot.title = element_text(hjust = 0.5),
                                                  plot.margin=unit(c(0, 0, 0, 0), "mm"))
}

low prop sig KO example: K00480. Ranks: 12, 10.5, 9, 6; Num. sig (cutoff 1): 221, 193, 208, 224
encoded in 540/3000 genomes

mid prop sig KO example: K01845. Ranks: 57, 52, 51, 52; Num sig (cutoff 1): 1407, 1393, 1396, 1366
encoded in 2292/3000 genomes

high prop sig KO example: K06077. Ranks: 92, 87, 86, 87; Num sig (cutoff 1): 3416, 3187, 3132, 3136
encoded in 828/3000 genomes


ERP002061_venn <- ggVennDiagram(almeida_DA_out_features_BY_0.05$ERP002061) +
                                scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
                                ggtitle("Primary obesity dataset sig. KOs") +
                                theme(plot.title = element_text(hjust = 0.5),
                                      plot.margin=unit(c(0, 0, 0, 0), "mm"))


almeida_POMS_out_ERP003612 <- readRDS("Almeida_2019_POMS_output/ERP003612_pseudo.null_sig.rds")
almeida_DA_out_features_BY_0.05$ERP003612$POMS <- c(almeida_POMS_out_ERP003612$KO_up, almeida_POMS_out_ERP003612$KO_down)
almeida_DA_out_features_BY_0.05$ERP003612[["Wilcoxon (USCG-norm.)"]] <- NULL
almeida_DA_out_features_BY_0.05$ERP003612[["DESeq2"]] <- NULL
ERP003612_venn <- ggVennDiagram(almeida_DA_out_features_BY_0.05$ERP003612) +
  scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
  ggtitle("Secondary obesity dataset sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5))


almeida_POMS_out_ERP012177 <- readRDS("Almeida_2019_POMS_output/ERP012177_pseudo.null_sig.rds")
almeida_DA_out_features_BY_0.05$ERP012177$POMS <- c(almeida_POMS_out_ERP012177$KO_up, almeida_POMS_out_ERP012177$KO_down)
almeida_DA_out_features_BY_0.05$ERP012177[["Wilcoxon (USCG-norm.)"]] <- NULL
almeida_DA_out_features_BY_0.05$ERP012177[["DESeq2"]] <- NULL
ERP012177_venn <- ggVennDiagram(almeida_DA_out_features_BY_0.05$ERP012177) +
                                scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
                                ggtitle("Colorectal cancer dataset sig. KOs") +
                                theme(plot.title = element_text(hjust = 0.5))



bottom_row <- plot_grid(ggplot() + theme_void(), ERP012177_venn, ggplot() + theme_void(),
                        rel_widths = c(0.25, 5, 0.25), nrow=1, ncol=3, labels=c('', 'c', ''), label_x = 0.225, label_y=1.05)
top_row <- plot_grid(ERP002061_venn, ERP003612_venn, labels=c('a', 'b'))

almeida_venn_plot <- plot_grid(top_row, bottom_row,nrow=2, ncol=1)

ggsave(filename = "../figures/almeida_DA_tool_venn_plot.pdf", plot = almeida_venn_plot,
       width = 10, height=8)



