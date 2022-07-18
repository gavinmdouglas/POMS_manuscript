rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(ggbeeswarm)

setwd("/home/gavin/github_repos/POMS_manuscript/data/results")

Almeida_POMS_out <- readRDS("Almeida_2019_POMS_output/combined_output.rds")

names(Almeida_POMS_out)[which(names(Almeida_POMS_out) == "ERP002061")] <- "Obesity 1"
names(Almeida_POMS_out)[which(names(Almeida_POMS_out) == "ERP003612")] <- "Obesity 2"
names(Almeida_POMS_out)[which(names(Almeida_POMS_out) == "ERP012177")] <- "Colorectal cancer"

TARA_POMS_out <- list()
TARA_POMS_out[["Mean salinity"]] <- readRDS("TARA_POMS_out.rds")[["Mean_Salinity"]]
TARA_POMS_out[["PO4"]] <- readRDS("TARA_POMS_out.rds")[["PO4"]]


Almeida_regress_out <- readRDS("Almeida_2019_regress_specificity_output/combined_output.rds")

names(Almeida_regress_out)[which(names(Almeida_regress_out) == "ERP002061")] <- "Obesity 1"
names(Almeida_regress_out)[which(names(Almeida_regress_out) == "ERP003612")] <- "Obesity 2"
names(Almeida_regress_out)[which(names(Almeida_regress_out) == "ERP012177")] <- "Colorectal cancer"

TARA_regress_out <- list()
TARA_regress_out[["Mean salinity"]] <- readRDS("TARA_regress_out.rds")[["Mean_Salinity"]]
TARA_regress_out[["PO4"]] <- readRDS("TARA_regress_out.rds")[["PO4"]]


Almeida_DA_out <- readRDS("Almeida_2019_DA_tool_output.rds")
names(Almeida_DA_out)[which(names(Almeida_DA_out) == "ERP002061")] <- "Obesity 1"
names(Almeida_DA_out)[which(names(Almeida_DA_out) == "ERP003612")] <- "Obesity 2"
names(Almeida_DA_out)[which(names(Almeida_DA_out) == "ERP012177")] <- "Colorectal cancer"

TARA_spearman_out <- readRDS("TARA_spearman_out.rds")
names(TARA_spearman_out)[which(names(TARA_spearman_out) == "Mean_Salinity")] <- "Mean salinity"

Almeida_sig_func_other_rankings_all <- data.frame(dataset = NA, func = NA, tool = NA, abs_rank = NA, percentile_rank = NA)
Almeida_sig_func_other_rankings_all <- Almeida_sig_func_other_rankings_all[-1, ]

for (dataset in names(Almeida_POMS_out)) {
  
  for (func_type in names(Almeida_POMS_out[[dataset]])) {
    
    #dataset <- "Obesity 1"
    #func_type <- "kos"
    
    sig_POMS_funcs <- rownames(Almeida_POMS_out[[dataset]][[func_type]]$results[which(Almeida_POMS_out[[dataset]][[func_type]]$results$multinomial_corr < 0.25), ])
    
    if (length(sig_POMS_funcs) > 0) {
    
      num_other_tools <- length(Almeida_DA_out[[dataset]][[func_type]]) + 1
      
      sig_func_other_rankings <- data.frame(matrix(NA, nrow = num_other_tools * length(sig_POMS_funcs), ncol = 6))
      colnames(sig_func_other_rankings) <- c("dataset", "func_type", "func", "tool", "abs_rank", "percentile_rank")
      
      current_row <- 1
      
      if (func_type == "kos") {
       DA_func_type <- "ko" 
      } else {
        DA_func_type <- func_type
      }
      
      regress_abs_ranks <- rank(Almeida_regress_out[[dataset]][[func_type]]$BH)
      names(regress_abs_ranks) <- rownames(Almeida_regress_out[[dataset]][[func_type]])
      
      regress_percent_ranks <- (regress_abs_ranks / length(regress_abs_ranks)) * 100
      
      for (sig_func in sig_POMS_funcs) {
        
        sig_func_other_rankings[current_row, c("dataset",  "func_type", "func", "tool")] <- c(dataset, func_type, sig_func, "regress") 
        
        if (sig_func %in% names(regress_abs_ranks)) {
          sig_func_other_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(regress_abs_ranks[sig_func],
                                                                                      regress_percent_ranks[sig_func])
        } else {
          sig_func_other_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(NA, NA)
        }
        
        current_row <- current_row + 1
        
      }
      
      for (DA_tool in names(Almeida_DA_out[[dataset]][[DA_func_type]])) {
        
        DA_tool_abs_ranks <- rank(Almeida_DA_out[[dataset]][[DA_func_type]][[DA_tool]]$BH_corr_p)
        names(DA_tool_abs_ranks) <- rownames(Almeida_DA_out[[dataset]][[DA_func_type]][[DA_tool]])
        
        DA_tool_percent_ranks <- (DA_tool_abs_ranks / length(DA_tool_abs_ranks)) * 100
        
        for (sig_func in sig_POMS_funcs) {
          
          sig_func_other_rankings[current_row, c("dataset",  "func_type", "func", "tool")] <- c(dataset, func_type, sig_func, DA_tool) 
          
          if (sig_func %in% names(DA_tool_abs_ranks)) {
            sig_func_other_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(DA_tool_abs_ranks[sig_func],
                                                                                     DA_tool_percent_ranks[sig_func])
          } else {
            sig_func_other_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(NA, NA)
          }
          
          current_row <- current_row + 1
    
        }
        
      }
      
      Almeida_sig_func_other_rankings_all <- rbind(Almeida_sig_func_other_rankings_all, sig_func_other_rankings)
    }
  }
}

Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$tool == "regress"), "tool"] <- "Phylo. regress.\n(specificity)"
Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$tool == "aldex2"), "tool"] <- "ALDEx2"
Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$tool == "deseq2"), "tool"] <- "DESeq2"
Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$tool == "limma.voom"), "tool"] <- "limma-voom"
Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$tool == "wilcoxon.musicc"), "tool"] <- "Wilcoxon (corr.)"
Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$tool == "wilcoxon.relab"), "tool"] <- "Wilcoxon (relab.)"

Almeida_sig_func_other_rankings_all$tool <- factor(Almeida_sig_func_other_rankings_all$tool,
                                                   levels = c("Phylo. regress.\n(specificity)",
                                                              "ALDEx2",
                                                              "DESeq2",
                                                              "limma-voom",
                                                              "Wilcoxon (corr.)",
                                                              "Wilcoxon (relab.)"))

Almeida_sig_func_other_rankings_ko <- Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$func_type == "kos"), ]
Almeida_ko_panel <- ggplot(Almeida_sig_func_other_rankings_ko, aes(x = tool, y = abs_rank)) +
                            geom_boxplot(outlier.colour = NA, fill = "grey90") +
                            geom_beeswarm(cex = 3) +
                            facet_grid(. ~ dataset) +
                            theme_bw() +
                            ggtitle("Human case-control: KEGG orthologs") +
                            xlab("") +
                            ylab("Rank") +
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                  legend.position = c(0.45, 0.7),
                                  legend.key = element_rect(fill = "white", colour = "white"),
                                  legend.box.background = element_rect(colour = "black"),
                                  legend.background = element_rect(colour = "light grey"),
                                  plot.title = element_text(hjust = 0.5))
                            
Almeida_sig_func_other_rankings_pathways <- Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$func_type == "pathways"), ]
Almeida_pathways_panel <- ggplot(Almeida_sig_func_other_rankings_pathways, aes(x = tool, y = abs_rank)) +
                                #geom_boxplot(outlier.colour = NA, fill = "grey90") +
                                geom_beeswarm(cex = 3) +
                                facet_grid(. ~ dataset) +
                                theme_bw() +
                                ggtitle("Human case-control: KEGG pathways") +
                                xlab("") +
                                ylab("Rank") +
                                ylim(c(0, 60)) +
                                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                      legend.position = c(0.45, 0.7),
                                      legend.key = element_rect(fill = "white", colour = "white"),
                                      legend.box.background = element_rect(colour = "black"),
                                      legend.background = element_rect(colour = "light grey"),
                                      plot.title = element_text(hjust = 0.5))

Almeida_sig_func_other_rankings_modules <- Almeida_sig_func_other_rankings_all[which(Almeida_sig_func_other_rankings_all$func_type == "modules"), ]
Almeida_modules_panel <- ggplot(Almeida_sig_func_other_rankings_modules, aes(x = tool, y = abs_rank)) +
                                geom_boxplot(outlier.colour = NA, fill = "grey90") +
                                geom_beeswarm(cex = 3) +
                                facet_grid(. ~ dataset) +
                                theme_bw() +
                                ggtitle("Human case-control: KEGG modules") +
                                xlab("") +
                                ylab("Rank") +
                                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                      legend.position = c(0.45, 0.7),
                                      legend.key = element_rect(fill = "white", colour = "white"),
                                      legend.box.background = element_rect(colour = "black"),
                                      legend.background = element_rect(colour = "light grey"),
                                      plot.title = element_text(hjust = 0.5))


# Similar, but for TARA
TARA_sig_func_rankings_all <- data.frame(dataset = NA, func_type = NA,  func = NA, tool = NA, abs_rank = NA, percentile_rank = NA)
TARA_sig_func_rankings_all <- TARA_sig_func_rankings_all[-1, ]


for (dataset in names(TARA_POMS_out)) {
  
  for (func_type in names(TARA_POMS_out[[dataset]])) {
    
    sig_POMS_funcs <- rownames(TARA_POMS_out[[dataset]][[func_type]]$results[which(TARA_POMS_out[[dataset]][[func_type]]$results$multinomial_corr < 0.25), ])
    
    if (length(sig_POMS_funcs) > 0) {
      
      sig_func_TARA_rankings <- data.frame(matrix(NA, nrow = length(sig_POMS_funcs) * 2, ncol = 6))
      colnames(sig_func_TARA_rankings) <- c("dataset", "func_type", "func", "tool", "abs_rank", "percentile_rank")
      
      current_row <- 1
      
      regress_abs_ranks <- rank(TARA_regress_out[[dataset]][[func_type]]$BH)
      names(regress_abs_ranks) <- rownames(TARA_regress_out[[dataset]][[func_type]])
      
      regress_percent_ranks <- (regress_abs_ranks / length(regress_abs_ranks)) * 100
      
      for (sig_func in sig_POMS_funcs) {
        
        sig_func_TARA_rankings[current_row, c("dataset",  "func_type", "func", "tool")] <- c(dataset, func_type, sig_func, "regress") 
        
        if (sig_func %in% names(regress_abs_ranks)) {
          sig_func_TARA_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(regress_abs_ranks[sig_func],
                                                                                      regress_percent_ranks[sig_func])
        } else {
          sig_func_TARA_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(NA, NA)
        }
        
        current_row <- current_row + 1
        
      }
      
      spearman_tool_abs_ranks <- rank(TARA_spearman_out[[dataset]][[func_type]]$p_corr)
      names(spearman_tool_abs_ranks) <- rownames(TARA_spearman_out[[dataset]][[func_type]])
      
      spearman_tool_percent_ranks <- (spearman_tool_abs_ranks / length(spearman_tool_abs_ranks)) * 100
      
      for (sig_func in sig_POMS_funcs) {
        
        sig_func_TARA_rankings[current_row, c("dataset",  "func_type", "func", "tool")] <- c(dataset, func_type, sig_func, "spearman") 
        
        if (sig_func %in% names(spearman_tool_abs_ranks)) {
          sig_func_TARA_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(spearman_tool_abs_ranks[sig_func],
                                                                                         spearman_tool_percent_ranks[sig_func])
        } else {
          sig_func_TARA_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(NA, NA)
        }
        
        current_row <- current_row + 1
        
      }

      TARA_sig_func_rankings_all <- rbind(TARA_sig_func_rankings_all, sig_func_TARA_rankings)
    }
  }
}


TARA_sig_func_rankings_all[which(TARA_sig_func_rankings_all$dataset == "PO4"), "dataset"] <- "Phosphate"
TARA_sig_func_rankings_all[which(TARA_sig_func_rankings_all$func_type == "ko"), "func_type"] <- "KEGG orthologs"
TARA_sig_func_rankings_all[which(TARA_sig_func_rankings_all$func_type == "module"), "func_type"] <- "KEGG modules"
TARA_sig_func_rankings_all[which(TARA_sig_func_rankings_all$func_type == "pathway"), "func_type"] <- "KEGG pathways"

TARA_sig_func_rankings_all$func_type <- factor(TARA_sig_func_rankings_all$func_type, levels = c("KEGG orthologs",
                                                                                      "KEGG pathways",
                                                                                      "KEGG modules"))

TARA_sig_func_rankings_all[which(TARA_sig_func_rankings_all$tool == "regress"), "tool"] <- "Phylo. regress.\n(specificity)"
TARA_sig_func_rankings_all[which(TARA_sig_func_rankings_all$tool == "spearman"), "tool"] <- "Spearman cor."

TARA_sig_func_rankings_all$tool <- factor(TARA_sig_func_rankings_all$tool,
                                                   levels = c("Phylo. regress.\n(specificity)",
                                                              "Spearman cor."))

TARA_sig_func_rankings_all_panel <- ggplot(TARA_sig_func_rankings_all, aes(x = tool, y = abs_rank, colour = dataset)) +
                                            geom_beeswarm(cex = 3) +
                                            theme_bw() +
                                            ggtitle("Tara oceans: All categories") +
                                            xlab("") +
                                            ylab("Rank") +
                                            scale_colour_manual(name = "Environmental\nfactor", values = c("black", "grey")) +
                                            facet_wrap(. ~ func_type, scales = "free_y") +
                                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                                  legend.key = element_rect(fill = "white", colour = "white"),
                                                  legend.box.background = element_rect(colour = "black"),
                                                  legend.background = element_rect(colour = "light grey"),
                                                  plot.title = element_text(hjust = 0.5))




POMS_sig_vs_other_rank_plot <- plot_grid(TARA_sig_func_rankings_all_panel, Almeida_ko_panel, Almeida_modules_panel, Almeida_pathways_panel,
                                    ncol = 2, nrow = 2, labels = c('a', 'b', 'c', 'd'))

ggsave(filename = "../../display_items/Supp_POMS_sig_vs_other_rank.pdf",
       plot = POMS_sig_vs_other_rank_plot,
       width = 11, height = 7, device = "pdf", dpi = 600)


ggsave(filename = "../../display_items/Supp_POMS_sig_vs_other_rank.png",
       plot = POMS_sig_vs_other_rank_plot,
       width = 11, height = 7, device = "png", dpi = 300)

