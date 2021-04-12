rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

nu0.50_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.50.rds")
nu0.65_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.65.rds")
nu0.80_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.80.rds")
nu0.95_sim_func_summary <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_nu0.95.rds")

ref.based_ranks_and_genome_nums_wilcoxon_nu0.50 <- ggplot(nu0.50_sim_func_summary, aes(y = num_focal_pos_mags, x = wilcoxon.musicc_rank_0.001)) +
                                                          geom_point(colour = "springgreen4", size = 2) +
                                                          theme_bw() +
                                                          xlab("") +
                                                          ylab("No. genomes encoding focal gene") +
                                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                          ggtitle("Wilcoxon test (nu = 0.50)") +
                                                          xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_wilcoxon_nu0.50_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_wilcoxon_nu0.50, type = "histogram", size = 10)


ref.based_ranks_and_genome_nums_POMS_nu0.50 <- ggplot(nu0.50_sim_func_summary, aes(y = num_focal_pos_mags, x = POMS_rank_0.25)) +
                                                      geom_point(colour = "steelblue4", size = 2) +
                                                      theme_bw() +
                                                      xlab("Focal gene rank") +
                                                      ylab("No. genomes encoding focal gene") +
                                                      theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                      ggtitle("POMS (nu = 0.50)") +
                                                      xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_POMS_nu0.50_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_POMS_nu0.50, type = "histogram", size = 10)





ref.based_ranks_and_genome_nums_wilcoxon_nu0.65 <- ggplot(nu0.65_sim_func_summary, aes(y = num_focal_pos_mags, x = wilcoxon.musicc_rank_0.001)) +
                                                          geom_point(colour = "springgreen4", size = 2) +
                                                          theme_bw() +
                                                          xlab("") +
                                                          ylab("") +
                                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                          ggtitle("Wilcoxon test (nu = 0.65)") +
                                                          xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_wilcoxon_nu0.65_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_wilcoxon_nu0.65, type = "histogram", size = 10)


ref.based_ranks_and_genome_nums_POMS_nu0.65 <- ggplot(nu0.65_sim_func_summary, aes(y = num_focal_pos_mags, x = POMS_rank_0.25)) +
                                                      geom_point(colour = "steelblue4", size = 2) +
                                                      theme_bw() +
                                                      xlab("Focal gene rank") +
                                                      ylab("") +
                                                      theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                      ggtitle("POMS (nu = 0.65)") +
                                                      xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_POMS_nu0.65_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_POMS_nu0.65, type = "histogram", size = 10)



ref.based_ranks_and_genome_nums_wilcoxon_nu0.80 <- ggplot(nu0.80_sim_func_summary, aes(y = num_focal_pos_mags, x = wilcoxon.musicc_rank_0.001)) +
                                                          geom_point(colour = "springgreen4", size = 2) +
                                                          theme_bw() +
                                                          xlab("") +
                                                          ylab("") +
                                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                          ggtitle("Wilcoxon test (nu = 0.80)") +
                                                          xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_wilcoxon_nu0.80_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_wilcoxon_nu0.80, type = "histogram", size = 10)


ref.based_ranks_and_genome_nums_POMS_nu0.80 <- ggplot(nu0.80_sim_func_summary, aes(y = num_focal_pos_mags, x = POMS_rank_0.25)) +
                                                      geom_point(colour = "steelblue4", size = 2) +
                                                      theme_bw() +
                                                      xlab("Focal gene rank") +
                                                      ylab("") +
                                                      theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                      ggtitle("POMS (nu = 0.80)") +
                                                      xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_POMS_nu0.80_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_POMS_nu0.80, type = "histogram", size = 10)




ref.based_ranks_and_genome_nums_wilcoxon_nu0.95 <- ggplot(nu0.95_sim_func_summary, aes(y = num_focal_pos_mags, x = wilcoxon.musicc_rank_0.001)) +
                                                          geom_point(colour = "springgreen4", size = 2) +
                                                          theme_bw() +
                                                          xlab("") +
                                                          ylab("") +
                                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                          ggtitle("Wilcoxon test (nu = 0.95)") +
                                                          xlim(0, 3100) + ylim(0, 500)

ref.based_ranks_and_genome_nums_wilcoxon_nu0.95_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_wilcoxon_nu0.95, type = "histogram", size = 10)


ref.based_ranks_and_genome_nums_POMS_nu0.95 <- ggplot(nu0.95_sim_func_summary, aes(y = num_focal_pos_mags, x = POMS_rank_0.25)) +
                                                      geom_point(colour = "steelblue4", size = 2) +
                                                      theme_bw() +
                                                      xlab("Focal gene rank") +
                                                      ylab("") +
                                                      theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 12)) +
                                                      ggtitle("POMS (nu = 0.95)") +
                                                      xlim(0, 3100) +
                                                      ylim(0, 500)
                                                      

ref.based_ranks_and_genome_nums_POMS_nu0.95_marginal <- ggMarginal(ref.based_ranks_and_genome_nums_POMS_nu0.95, type = "histogram", size = 10)


ref.based_ranks_and_genome_nums_plots <- plot_grid(ref.based_ranks_and_genome_nums_wilcoxon_nu0.50_marginal,
                                                   ref.based_ranks_and_genome_nums_wilcoxon_nu0.65_marginal,
                                                   ref.based_ranks_and_genome_nums_wilcoxon_nu0.80_marginal,
                                                   ref.based_ranks_and_genome_nums_wilcoxon_nu0.95_marginal,
                                                   ref.based_ranks_and_genome_nums_POMS_nu0.50_marginal,
                                                   ref.based_ranks_and_genome_nums_POMS_nu0.65_marginal,
                                                   ref.based_ranks_and_genome_nums_POMS_nu0.80_marginal,
                                                   ref.based_ranks_and_genome_nums_POMS_nu0.95_marginal,
                                                   nrow = 2,
                                                   ncol = 4,
                                                   labels = c('a', 'b', 'c', 'd', 'e', 'f' ,'g', 'h'))


ggsave(filename = "~/github_repos/POMS_manuscript/figures/Supp_ref.based_sims_ranking.pdf",
       plot = ref.based_ranks_and_genome_nums_plots,
       device = "pdf",
       width = 12,
       height = 6,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/figures/Supp_ref.based_sims_ranking.png",
       plot = ref.based_ranks_and_genome_nums_plots,
       device = "png",
       width = 12,
       height = 6,
       dpi = 300)
