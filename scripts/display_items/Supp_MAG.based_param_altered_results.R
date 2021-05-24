rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

library(cowplot)
library(ggplot2)

summary_df <- readRDS(file = "POMS_wilcoxon.musicc_out_summary.rds")

summary_df_clean <- summary_df[, c("rep", "MAGs", "pseudocount", "abun_increase", "num_focal_pos_mags",
                                   "POMS_num_sig_nodes", "POMS_sig_0.05", "POMS_rank_0.05", "POMS_rel_rank_0.05", "POMS_jaccard_0.05",
                                   "wilcoxon.musicc_sig_0.05", "wilcoxon.musicc_rank_0.05",
                                   "wilcoxon.musicc_rel_rank_0.05", "wilcoxon.musicc_jaccard_0.05")]

summary_df_clean_mean <- aggregate(. ~ MAGs * pseudocount * abun_increase,
                                   data = summary_df_clean, FUN = mean, na.rm = TRUE, na.action = NULL)

summary_df_clean_mean <- summary_df_clean_mean[order(summary_df_clean_mean$pseudocount, summary_df_clean_mean$abun_increase, decreasing = TRUE), ]


summary_df_clean_mean$sel_set <- paste("pseudo", as.character(summary_df_clean_mean$pseudocount), "; ",
                                       "abun*", as.character(summary_df_clean_mean$abun_increase), sep = "")

summary_df_clean_mean$sel_set <- factor(summary_df_clean_mean$sel_set,
                                        levels = rev(summary_df_clean_mean$sel_set[-which(duplicated(summary_df_clean_mean$sel_set))]))

summary_df_clean_mean$MAG_num_char <- as.character(summary_df_clean_mean$MAGs)
summary_df_clean_mean$MAG_num_char <- factor(summary_df_clean_mean$MAG_num_char,
                                             levels = c("50", "100", "250", "500", "750", "1000", "1250", "1595"))


summary_df_clean_median <- aggregate(. ~ MAGs * pseudocount * abun_increase,
                                   data = summary_df_clean, FUN = median, na.rm = TRUE, na.action = NULL)

summary_df_clean_median <- summary_df_clean_median[order(summary_df_clean_median$pseudocount, summary_df_clean_median$abun_increase, decreasing = TRUE), ]

summary_df_clean_median$sel_set <- paste("pseudo", as.character(summary_df_clean_median$pseudocount), "; ",
                                         "abun*", as.character(summary_df_clean_median$abun_increase), sep = "")

summary_df_clean_median$sel_set <- factor(summary_df_clean_median$sel_set,
                                        levels = rev(summary_df_clean_median$sel_set[-which(duplicated(summary_df_clean_median$sel_set))]))

summary_df_clean_median$MAG_num_char <- as.character(summary_df_clean_median$MAGs)
summary_df_clean_median$MAG_num_char <- factor(summary_df_clean_median$MAG_num_char,
                                             levels = c("50", "100", "250", "500", "750", "1000", "1250", "1595"))


num_sig_nodes_panel <- ggplot(data = summary_df_clean_mean, aes(y = sel_set, x = MAG_num_char, fill = POMS_num_sig_nodes)) +
                              geom_tile() +
                              geom_text(aes(label = round(POMS_num_sig_nodes, 2))) +
                              ylab("Selection setting") +
                              xlab("Number of MAGs") +
                              theme_bw() +
                              scale_fill_continuous(name = "Mean no.\nsig nodes",
                                                    low = "light green", high = "dark green",
                                                    limits = c(0, 32))

POMS_sig_prop_panel <- ggplot(data = summary_df_clean_mean, aes(y = sel_set, x = MAG_num_char, fill = POMS_sig_0.05)) +
                             geom_tile() +
                             geom_text(aes(label = round(POMS_sig_0.05, 3))) +
                             ylab("Selection setting") +
                             xlab("Number of MAGs") +
                             theme_bw() +
                             scale_fill_continuous(name = "Mean prop.\nsig KOs", low = "coral1", high = "coral4",
                                                   limits = c(0, 1)) +
                             ggtitle("POMS") +
                             theme(plot.title = element_text(hjust = 0.5, vjust = -1))

POMS_rank_panel <- ggplot(data = summary_df_clean_median, aes(y = sel_set, x = MAG_num_char, fill = POMS_rank_0.05)) +
                           geom_tile() +
                           geom_text(aes(label = POMS_rank_0.05)) +
                           ylab("") +
                           xlab("Number of MAGs") +
                           theme_bw() +
                           scale_fill_gradient(name = "Median focal\ngene rank",
                                               low = "dodgerblue", high = "dodgerblue4", limits = c(0, 1260)) +
                           ggtitle("POMS") +
                           theme(plot.title = element_text(hjust = 0.5, vjust = -1))

wilcoxon.musicc_sig_prop_panel <- ggplot(data = summary_df_clean_mean, aes(y = sel_set, x = MAG_num_char, fill = wilcoxon.musicc_sig_0.05)) +
                                          geom_tile() +
                                          geom_text(aes(label = round(wilcoxon.musicc_sig_0.05, 3))) +
                                          ylab("Selection setting") +
                                          xlab("") +
                                          theme_bw() +
                                          scale_fill_continuous(name = "Mean prop.\nsig KOs", low = "coral1", high = "coral4",
                                                                limits = c(0, 1)) +
                                          ggtitle("Wilcoxon test") +
                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1))

wilcoxon.musicc_rank_panel <- ggplot(data = summary_df_clean_median, aes(y = sel_set, x = MAG_num_char, fill = wilcoxon.musicc_rank_0.05)) +
                                    geom_tile() +
                                    geom_text(aes(label = wilcoxon.musicc_rank_0.05)) +
                                    ylab("") +
                                    xlab("") +
                                    theme_bw() +
                                    scale_fill_gradient(name = "Median focal\ngene rank",
                                                        low = "dodgerblue", high = "dodgerblue4", limits = c(0, 1260)) +
                                    ggtitle("Wilcoxon test") +
                                    theme(plot.title = element_text(hjust = 0.5, vjust = -1))


ranking_param_altered_plot <- plot_grid(wilcoxon.musicc_rank_panel,
                                        POMS_rank_panel,
                                        labels = c('a', 'b'),
                                        nrow = 2, ncol = 1)


sig_prop_param_altered_plot <- plot_grid(wilcoxon.musicc_sig_prop_panel,
                                         POMS_sig_prop_panel,
                                         labels = c('a', 'b'),
                                         nrow = 2, ncol = 1)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking_param_altered_heatmaps.pdf",
       plot = ranking_param_altered_plot,
       device = "pdf",
       width = 7,
       height = 10,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking_param_altered_heatmaps.png",
       plot = ranking_param_altered_plot,
       device = "png",
       width = 7,
       height = 10,
       dpi = 300)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_sig_prop_param_altered_heatmaps.pdf",
       plot = sig_prop_param_altered_plot,
       device = "pdf",
       width = 7,
       height = 10,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_sig_prop_param_altered_heatmaps.png",
       plot = sig_prop_param_altered_plot,
       device = "png",
       width = 7,
       height = 10,
       dpi = 300)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_param_altered_num_sig_nodes.pdf",
       plot = num_sig_nodes_panel,
       device = "pdf",
       width = 7,
       height = 6,
       dpi = 600)

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_param_altered_num_sig_nodes.png",
       plot = num_sig_nodes_panel,
       device = "png",
       width = 7,
       height = 6,
       dpi = 300)
