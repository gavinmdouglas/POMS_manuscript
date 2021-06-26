rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

# Genome prevalence plot
nu0.50_abun_tab <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.50_sigma1.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
nu0.65_abun_tab <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.65_sigma1.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
nu0.80_abun_tab <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.80_sigma1.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
nu0.95_abun_tab <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.95_sigma1.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
MAG.based_abun_tab <- readRDS("../MAG.based_simulations/MAG.based_prepped_abun.rds")

genome_prev <- data.frame(dataset = c(rep("nu=0.50", nrow(nu0.50_abun_tab)),
                                      rep("nu=0.65", nrow(nu0.65_abun_tab)),
                                      rep("nu=0.80", nrow(nu0.80_abun_tab)),
                                      rep("nu=0.95", nrow(nu0.95_abun_tab)),
                                      rep("MAG-based", nrow(MAG.based_abun_tab))),
                          
                          prev = c((rowSums(nu0.50_abun_tab > 0) / ncol(nu0.50_abun_tab)) * 100,
                                   (rowSums(nu0.65_abun_tab > 0) / ncol(nu0.65_abun_tab)) * 100,
                                   (rowSums(nu0.80_abun_tab > 0) / ncol(nu0.80_abun_tab)) * 100,
                                   (rowSums(nu0.95_abun_tab > 0) / ncol(nu0.95_abun_tab)) * 100,
                                   (rowSums(MAG.based_abun_tab > 0) / ncol(MAG.based_abun_tab)) * 100),
                          
                          fill = c(rep("Reference genomes",  nrow(nu0.50_abun_tab) +  nrow(nu0.65_abun_tab) + nrow(nu0.80_abun_tab) + nrow(nu0.95_abun_tab)),
                                   rep("MAGs", nrow(MAG.based_abun_tab))))

genome_prev$dataset <- factor(genome_prev$dataset, levels = c("nu=0.50", "nu=0.65", "nu=0.80", "nu=0.95", "MAG-based"))
genome_prev$fill <- factor(genome_prev$fill, levels = c("Reference genomes", "MAGs"))

genome_prev_boxplots <- ggplot(genome_prev, aes(x = dataset, y = prev, fill = fill)) +
                                geom_boxplot() +
                                ylim(0, 100) +
                                ylab("Genome prevalence across samples") +
                                xlab("Abundance table sim. approach") +
                                scale_fill_manual(name = "Dataset", values = c("black", "light blue")) +
                                theme_bw() +
                                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                      legend.position = c(0.45, 0.7),
                                      legend.key = element_rect(fill = "white", colour = "white"),
                                      legend.box.background = element_rect(colour = "black"),
                                      legend.background = element_rect(colour = "light grey"))




# Ranking plot
nu0.50_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.50.rds")
nu0.65_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.65.rds")
nu0.80_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.80.rds")
nu0.95_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.95.rds")
MAG.based_sim_func_summary <- readRDS(file = "../MAG.based_simulations/simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")


rankings <- data.frame(dataset = c(rep("nu=0.50", nrow(nu0.50_sim_func_summary) * 2),
                                      rep("nu=0.65", nrow(nu0.65_sim_func_summary) * 2),
                                      rep("nu=0.80", nrow(nu0.80_sim_func_summary) * 2),
                                      rep("nu=0.95", nrow(nu0.95_sim_func_summary) * 2),
                                      rep("MAG-based", nrow(MAG.based_sim_func_summary) * 2)),
                       
                       Method = c(rep("Wilcoxon test", nrow(nu0.50_sim_func_summary)),
                                  rep("POMS", nrow(nu0.50_sim_func_summary)),
                                    
                                  rep("Wilcoxon test", nrow(nu0.65_sim_func_summary)),
                                  rep("POMS", nrow(nu0.65_sim_func_summary)),
                                  
                                  rep("Wilcoxon test", nrow(nu0.80_sim_func_summary)),
                                  rep("POMS", nrow(nu0.80_sim_func_summary)),
                                  
                                  rep("Wilcoxon test", nrow(nu0.95_sim_func_summary)),
                                  rep("POMS", nrow(nu0.95_sim_func_summary)),
                       
                                 rep("Wilcoxon test", nrow(MAG.based_sim_func_summary)),
                                 rep("POMS", nrow(MAG.based_sim_func_summary))),
                          
                          rank = c(nu0.50_sim_func_summary$wilcoxon.musicc_rank_0.05,
                                   nu0.50_sim_func_summary$POMS_rank_0.05,
                                   
                                   nu0.65_sim_func_summary$wilcoxon.musicc_rank_0.05,
                                   nu0.65_sim_func_summary$POMS_rank_0.05,
                                   
                                   nu0.80_sim_func_summary$wilcoxon.musicc_rank_0.05,
                                   nu0.80_sim_func_summary$POMS_rank_0.05,
                                   
                                   nu0.95_sim_func_summary$wilcoxon.musicc_rank_0.05,
                                   nu0.95_sim_func_summary$POMS_rank_0.05,
                                   
                                   MAG.based_sim_func_summary$wilcoxon.musicc_rank_0.05,
                                   MAG.based_sim_func_summary$POMS_rank_0.05))

rankings$dataset <- factor(rankings$dataset, levels = c("nu=0.50", "nu=0.65", "nu=0.80", "nu=0.95", "MAG-based"))

rankings_boxplots <- ggplot(rankings, aes(x = dataset, y = rank, fill = Method)) +
                                              geom_boxplot() +
                                              ylab("Focal gene rank") +
                                              xlab("Abundance table sim. approach") +
                                              scale_fill_manual(name = "Method", values = c("steelblue4", "springgreen4")) +
                                              theme_bw() +
                                              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                                    legend.position = c(0.157, 0.7),
                                                    legend.key = element_rect(fill = "white", colour = "white"),
                                                    legend.box.background = element_rect(colour = "black"),
                                                    legend.background = element_rect(colour = "light grey"))
                                                                                      
genome_prev_and_ranking_plot <- plot_grid(genome_prev_boxplots, rankings_boxplots, nrow = 1, ncol = 2, labels = c('a', 'b'), rel_widths = c(1, 2))


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_ref.based_abun_and_ranking_POMS_vs_wilcoxon.pdf",
       plot = genome_prev_and_ranking_plot,
       device = "pdf",
       width = 8,
       height = 4,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_ref.based_abun_and_ranking_POMS_vs_wilcoxon.png",
       plot = genome_prev_and_ranking_plot,
       device = "png",
       width = 8,
       height = 4,
       dpi = 300)


# Summary stats

wilcox.test(nu0.50_sim_func_summary$POMS_rank_0.05, nu0.50_sim_func_summary$wilcoxon.musicc_rank_0.05)
wilcox.test(nu0.65_sim_func_summary$POMS_rank_0.05, nu0.65_sim_func_summary$wilcoxon.musicc_rank_0.05)
wilcox.test(nu0.80_sim_func_summary$POMS_rank_0.05, nu0.80_sim_func_summary$wilcoxon.musicc_rank_0.05)
wilcox.test(nu0.95_sim_func_summary$POMS_rank_0.05, nu0.95_sim_func_summary$wilcoxon.musicc_rank_0.05)

cor.test(nu0.95_sim_func_summary$wilcoxon.musicc_rank_0.05, nu0.95_sim_func_summary$num_focal_pos_mags)
