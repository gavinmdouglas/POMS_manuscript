rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)

setwd("/home/gavin/github_repos/POMS_manuscript/data/results")

Almeida_POMS_out <- list()
Almeida_POMS_out[["Obesity 1"]] <- readRDS("Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
Almeida_POMS_out[["Obesity 2"]] <- readRDS("Almeida_2019_POMS_output/ERP003612_POMS_out.rds")
Almeida_POMS_out[["Colorectal cancer"]] <- readRDS("Almeida_2019_POMS_output/ERP012177_POMS_out.rds")

TARA_POMS_out <- list()
TARA_POMS_out[["Mean salinity"]] <- readRDS("TARA_POMS_out.rds")[["Mean_Salinity"]]
TARA_POMS_out[["PO4"]] <- readRDS("TARA_POMS_out.rds")[["PO4"]]

Almeida_DA_out <- readRDS("Almeida_2019_DA_tool_output.rds")
names(Almeida_DA_out)[which(names(Almeida_DA_out) == "ERP002061")] <- "Obesity 1"
names(Almeida_DA_out)[which(names(Almeida_DA_out) == "ERP003612")] <- "Obesity 2"
names(Almeida_DA_out)[which(names(Almeida_DA_out) == "ERP012177")] <- "Colorectal cancer"

TARA_spearman_out <- readRDS("TARA_spearman_out.rds")
names(TARA_spearman_out)[which(names(TARA_spearman_out) == "Mean_Salinity")] <- "Mean salinity"

Almeida_sig_func_DA_rankings_all <- data.frame(dataset = NA, func = NA, tool = NA, abs_rank = NA, percentile_rank = NA)
Almeida_sig_func_DA_rankings_all <- Almeida_sig_func_DA_rankings_all[-1, ]

for (dataset in names(Almeida_POMS_out)) {
  
  for (func_type in names(Almeida_POMS_out[[dataset]])) {
    
    #dataset <- "Obesity 1"
    #func_type <- "ko"
    
    sig_funcs <- rownames(Almeida_POMS_out[[dataset]][[func_type]]$df[which(Almeida_POMS_out[[dataset]][[func_type]]$df$multinomial_corr < 0.25), ])
    
    if (length(sig_funcs) > 0) {
    
      num_DA_tools <- length(Almeida_DA_out[[dataset]][[func_type]])
      
      sig_func_DA_rankings <- data.frame(matrix(NA, nrow = num_DA_tools * length(sig_funcs), ncol = 6))
      colnames(sig_func_DA_rankings) <- c("dataset", "func_type", "func", "tool", "abs_rank", "percentile_rank")
      
      current_row <- 1
      
      for (DA_tool in names(Almeida_DA_out[[dataset]][[func_type]])) {
        
        DA_tool_abs_ranks <- rank(Almeida_DA_out[[dataset]][[func_type]][[DA_tool]]$BH_corr_p)
        names(DA_tool_abs_ranks) <- rownames(Almeida_DA_out[[dataset]][[func_type]][[DA_tool]])
        
        DA_tool_percent_ranks <- (DA_tool_abs_ranks / length(DA_tool_abs_ranks)) * 100
        
        for (sig_func in sig_funcs) {
          
          sig_func_DA_rankings[current_row, c("dataset",  "func_type", "func", "tool")] <- c(dataset, func_type, sig_func, DA_tool) 
          
          if (sig_func %in% names(DA_tool_abs_ranks)) {
            sig_func_DA_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(DA_tool_abs_ranks[sig_func],
                                                                                     DA_tool_percent_ranks[sig_func])
          } else {
            sig_func_DA_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(NA, NA)
          }
          
          current_row <- current_row + 1
    
        }
        
      }
      
      Almeida_sig_func_DA_rankings_all <- rbind(Almeida_sig_func_DA_rankings_all, sig_func_DA_rankings)
    }
  }
}

Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$tool == "aldex2"), "tool"] <- "ALDEx2"
Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$tool == "deseq2"), "tool"] <- "DESeq2"
Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$tool == "limma.voom"), "tool"] <- "limma-voom"
Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$tool == "wilcoxon.musicc"), "tool"] <- "Wilcoxon (MUSiCC-corr.)"
Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$tool == "wilcoxon.relab"), "tool"] <- "Wilcoxon (relab.)"



Almeida_sig_func_DA_rankings_ko <- Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$func_type == "ko"), ]
Almeida_ko_panel <- ggplot(Almeida_sig_func_DA_rankings_ko, aes(x = tool, y = abs_rank)) +
                            geom_boxplot() +
                            facet_grid(. ~ dataset) +
                            theme_bw() +
                            ggtitle("KEGG orthologs") +
                            xlab("") +
                            ylab("Rank") +
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                  legend.position = c(0.45, 0.7),
                                  legend.key = element_rect(fill = "white", colour = "white"),
                                  legend.box.background = element_rect(colour = "black"),
                                  legend.background = element_rect(colour = "light grey"),
                                  plot.title = element_text(hjust = 0.5))
                            
Almeida_sig_func_DA_rankings_pathways <- Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$func_type == "pathways"), ]
Almeida_pathways_panel <- ggplot(Almeida_sig_func_DA_rankings_pathways, aes(x = tool, y = abs_rank)) +
  geom_boxplot() +
  facet_grid(. ~ dataset) +
  theme_bw() +
  ggtitle("KEGG pathways") +
  xlab("") +
  ylab("Rank") +
  ylim(c(0, 60)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = c(0.45, 0.7),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_rect(colour = "light grey"),
        plot.title = element_text(hjust = 0.5))

Almeida_sig_func_DA_rankings_modules <- Almeida_sig_func_DA_rankings_all[which(Almeida_sig_func_DA_rankings_all$func_type == "modules"), ]
Almeida_modules_panel <- ggplot(Almeida_sig_func_DA_rankings_modules, aes(x = tool, y = abs_rank)) +
  geom_boxplot() +
  facet_grid(. ~ dataset) +
  theme_bw() +
  ggtitle("KEGG modules") +
  xlab("") +
  ylab("Rank") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = c(0.45, 0.7),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_rect(colour = "light grey"),
        plot.title = element_text(hjust = 0.5))

TARA_sig_func_spearman_rankings_all <- data.frame(dataset = NA, func = NA, tool = NA, abs_rank = NA, percentile_rank = NA)
TARA_sig_func_spearman_rankings_all <- TARA_sig_func_spearman_rankings_all[-1, ]


for (dataset in names(TARA_POMS_out)) {
  
  for (func_type in names(TARA_POMS_out[[dataset]])) {
    
    sig_funcs <- rownames(TARA_POMS_out[[dataset]][[func_type]]$df[which(TARA_POMS_out[[dataset]][[func_type]]$df$multinomial_corr < 0.25), ])
    
    if (length(sig_funcs) > 0) {
      
      sig_func_spearman_rankings <- data.frame(matrix(NA, nrow = length(sig_funcs), ncol = 5))
      colnames(sig_func_spearman_rankings) <- c("dataset", "func_type", "func", "abs_rank", "percentile_rank")
      
      current_row <- 1
      
      spearman_tool_abs_ranks <- rank(TARA_spearman_out[[dataset]][[func_type]]$p_corr)
      names(spearman_tool_abs_ranks) <- rownames(TARA_spearman_out[[dataset]][[func_type]])
      
      spearman_tool_percent_ranks <- (spearman_tool_abs_ranks / length(spearman_tool_abs_ranks)) * 100
      
      for (sig_func in sig_funcs) {
        
        sig_func_spearman_rankings[current_row, c("dataset",  "func_type", "func")] <- c(dataset, func_type, sig_func) 
        
        if (sig_func %in% names(spearman_tool_abs_ranks)) {
          sig_func_spearman_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(spearman_tool_abs_ranks[sig_func],
                                                                                         spearman_tool_percent_ranks[sig_func])
        } else {
          sig_func_spearman_rankings[current_row, c("abs_rank", "percentile_rank")] <- c(NA, NA)
        }
        
        current_row <- current_row + 1
        
      }

      TARA_sig_func_spearman_rankings_all <- rbind(TARA_sig_func_spearman_rankings_all, sig_func_spearman_rankings)
    }
  }
}


TARA_sig_func_spearman_rankings_ko <- TARA_sig_func_spearman_rankings_all[which(TARA_sig_func_spearman_rankings_all$func_type == "ko"), ]

TARA_sig_func_ranks_ko_panel <- ggplot(TARA_sig_func_spearman_rankings_ko, aes(x = dataset, y = abs_rank)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("TARA oceans KEGG orthologs") +
  xlab("") +
  ylab("Rank") +
  facet_grid(. ~ func_type, scales = "free_y") +
  theme(legend.position = c(0.45, 0.7),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_rect(colour = "light grey"),
        plot.title = element_text(hjust = 0.5))

plot_grid(Almeida_ko_panel, Almeida_pathways_panel, Almeida_modules_panel, ncol = 2, nrow = 2)


# sig_abun_vs_enrich_plot <- ,
#                                      rel_widths = c(0.42, 0.42, 0.16))
# 
# ggsave(filename = "../../display_items/Supp_obesity_compare_enrichment_vs_relabun.pdf",
#        plot = sig_abun_vs_enrich_plot,
#        width = 8.5, height = 3.5, device = "pdf", dpi = 600)
# 
# 
# ggsave(filename = "../../display_items/Supp_obesity_compare_enrichment_vs_relabun.png",
#        plot = sig_abun_vs_enrich_plot,
#        width = 8.5, height = 3.5, device = "png", dpi = 300)
# 
