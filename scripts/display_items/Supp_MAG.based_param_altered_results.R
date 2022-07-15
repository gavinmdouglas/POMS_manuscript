rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

library(cowplot)
library(ggplot2)
library(reshape2)

func.based_summary <- readRDS(file = "summaries/MAG.based_func.based_param.altered_summary.rds")

func.based_summary <- func.based_summary[, which(colnames(func.based_summary) != "focal_func")]

func.based_summary_mean <- aggregate(. ~ MAGs * pseudocount * abun_increase,
                                     data = func.based_summary, FUN = mean, na.rm = TRUE, na.action = NULL)

func.based_summary_mean <- func.based_summary_mean[order(func.based_summary_mean$pseudocount, func.based_summary_mean$abun_increase, decreasing = TRUE), ]


func.based_summary_mean$sel_set <- paste("pseudo", as.character(func.based_summary_mean$pseudocount), "; ",
                                       "abun*", as.character(func.based_summary_mean$abun_increase), sep = "")

func.based_summary_mean$sel_set <- factor(func.based_summary_mean$sel_set,
                                        levels = rev(func.based_summary_mean$sel_set[-which(duplicated(func.based_summary_mean$sel_set))]))

func.based_summary_mean$MAG_num_char <- as.character(func.based_summary_mean$MAGs)
func.based_summary_mean$MAG_num_char <- factor(func.based_summary_mean$MAG_num_char,
                                             levels = c("100", "250", "500", "1000", "1595"))


func.based_summary_median <- aggregate(. ~ MAGs * pseudocount * abun_increase,
                                     data = func.based_summary, FUN = median, na.rm = TRUE, na.action = NULL)

func.based_summary_median <- func.based_summary_median[order(func.based_summary_median$pseudocount, func.based_summary_median$abun_increase, decreasing = TRUE), ]


func.based_summary_median$sel_set <- paste("pseudo", as.character(func.based_summary_median$pseudocount), "; ",
                                         "abun*", as.character(func.based_summary_median$abun_increase), sep = "")

func.based_summary_median$sel_set <- factor(func.based_summary_median$sel_set,
                                          levels = rev(func.based_summary_median$sel_set[-which(duplicated(func.based_summary_median$sel_set))]))

func.based_summary_median$MAG_num_char <- as.character(func.based_summary_median$MAGs)
func.based_summary_median$MAG_num_char <- factor(func.based_summary_median$MAG_num_char,
                                               levels = c("100", "250", "500", "1000", "1595"))


func.based_summary_mean_melt <- melt(data = func.based_summary_mean, id.vars = c("sel_set", "MAG_num_char"))

func.based_summary_mean_melt_sig_prop <- func.based_summary_mean_melt[grep("_sig_0.05$", func.based_summary_mean_melt$variable), ]

func.based_summary_mean_melt_sig_prop$variable <- as.character(func.based_summary_mean_melt_sig_prop$variable)
func.based_summary_mean_melt_sig_prop$Approach <- NA
func.based_summary_mean_melt_sig_prop$Approach[which(func.based_summary_mean_melt_sig_prop$variable == "POMS_sig_0.05")] <- "POMS"
func.based_summary_mean_melt_sig_prop$Approach[which(func.based_summary_mean_melt_sig_prop$variable == "regress_specificity_func_prop_sig_0.05")] <- "Phylo. regression (specificity)"
func.based_summary_mean_melt_sig_prop$Approach[which(func.based_summary_mean_melt_sig_prop$variable == "regress_sig_taxa_func_prop_sig_0.05")] <- "Phylo. regression (sig. taxa)"
func.based_summary_mean_melt_sig_prop$Approach[which(func.based_summary_mean_melt_sig_prop$variable == "wilcoxon.musicc_sig_0.05")] <- "Wilcoxon test (corrected)"


func.based_sig_prop_heatmaps <- ggplot(data = func.based_summary_mean_melt_sig_prop, aes(y = sel_set, x = MAG_num_char, fill = value)) +
                                geom_tile() +
                                geom_text(aes(label = sprintf("%0.2f", value)), size = 2.5) +
                                ylab("Selection setting") +
                                xlab("Number of MAGs") +
                                theme_bw() +
                                scale_fill_continuous(name = "Mean prop.\nsig. KOs", low = "coral1", high = "coral4",
                                                      limits = c(0, 1)) +
                                facet_wrap(. ~ Approach) +
                                theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                ggtitle("Focal function simulations")


func.based_summary_median_melt <- melt(data = func.based_summary_median, id.vars = c("sel_set", "MAG_num_char"))

func.based_summary_median_melt_rank <- func.based_summary_median_melt[grep("_rank_0.05$", func.based_summary_median_melt$variable), ]
func.based_summary_median_melt_rank <- func.based_summary_median_melt_rank[-grep("_rel_rank_0.05$", func.based_summary_median_melt_rank$variable), ]

func.based_summary_median_melt_rank$variable <- as.character(func.based_summary_median_melt_rank$variable)
func.based_summary_median_melt_rank$Approach <- NA
func.based_summary_median_melt_rank$Approach[which(func.based_summary_median_melt_rank$variable == "POMS_rank_0.05")] <- "POMS"
func.based_summary_median_melt_rank$Approach[which(func.based_summary_median_melt_rank$variable == "regress_specificity_rank_0.05")] <- "Phylogenetic regression (specificity)"
func.based_summary_median_melt_rank$Approach[which(func.based_summary_median_melt_rank$variable == "regress_sig_taxa_rank_0.05")] <- "Phylogenetic regression (sig. taxa)"
func.based_summary_median_melt_rank$Approach[which(func.based_summary_median_melt_rank$variable == "wilcoxon.musicc_rank_0.05")] <- "Wilcoxon test (corrected)"


func.based_rank_heatmaps <-  ggplot(data = func.based_summary_median_melt_rank, aes(y = sel_set, x = MAG_num_char, fill = value)) +
                                                geom_tile() +
                                                geom_text(aes(label = value)) +
                                                ylab("") +
                                                xlab("Number of metagenome-assembled genomes") +
                                                theme_bw() +
                                                scale_fill_gradient(name = "Median focal\ngene rank",
                                                                    low = "dodgerblue", high = "dodgerblue4", limits = c(0, 1050)) +
                                                theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                                facet_wrap(. ~ Approach) +
                                                ggtitle("Focal function simulations")



clade.based_summary <- readRDS(file = "summaries/MAG.based_clade.based_param.altered_summary.rds")

clade.based_summary <- clade.based_summary[, which(colnames(clade.based_summary) != "focal_func")]

clade.based_summary_mean <- aggregate(. ~ MAGs * pseudocount * abun_increase,
                                     data = clade.based_summary, FUN = mean, na.rm = TRUE, na.action = NULL)

clade.based_summary_mean <- clade.based_summary_mean[order(clade.based_summary_mean$pseudocount, clade.based_summary_mean$abun_increase, decreasing = TRUE), ]


clade.based_summary_mean$sel_set <- paste("pseudo", as.character(clade.based_summary_mean$pseudocount), "; ",
                                         "abun*", as.character(clade.based_summary_mean$abun_increase), sep = "")

clade.based_summary_mean$sel_set <- factor(clade.based_summary_mean$sel_set,
                                          levels = rev(clade.based_summary_mean$sel_set[-which(duplicated(clade.based_summary_mean$sel_set))]))

clade.based_summary_mean$MAG_num_char <- as.character(clade.based_summary_mean$MAGs)
clade.based_summary_mean$MAG_num_char <- factor(clade.based_summary_mean$MAG_num_char,
                                               levels = c("100", "250", "500", "1000", "1595"))

clade.based_summary_mean_melt <- melt(data = clade.based_summary_mean, id.vars = c("sel_set", "MAG_num_char"))

clade.based_summary_mean_melt_sig_prop <- clade.based_summary_mean_melt[grep("0.05$", clade.based_summary_mean_melt$variable), ]

clade.based_summary_mean_melt_sig_prop$variable <- as.character(clade.based_summary_mean_melt_sig_prop$variable)
clade.based_summary_mean_melt_sig_prop$Approach <- NA
clade.based_summary_mean_melt_sig_prop$Approach[which(clade.based_summary_mean_melt_sig_prop$variable == "POMS_sig_0.05")] <- "POMS"
clade.based_summary_mean_melt_sig_prop$Approach[which(clade.based_summary_mean_melt_sig_prop$variable == "regress_specificity_sig_0.05")] <- "Phylo. regression (specificity)"
clade.based_summary_mean_melt_sig_prop$Approach[which(clade.based_summary_mean_melt_sig_prop$variable == "regress_sig.taxa_sig_0.05")] <- "Phylo. regression (sig. taxa)"
clade.based_summary_mean_melt_sig_prop$Approach[which(clade.based_summary_mean_melt_sig_prop$variable == "wilcoxon.musicc_0.05")] <- "Wilcoxon test (corrected)"


clade.based_sig_prop_heatmaps <- ggplot(data = clade.based_summary_mean_melt_sig_prop, aes(y = sel_set, x = MAG_num_char, fill = value)) +
                                        geom_tile() +
                                        geom_text(aes(label = sprintf("%0.2f", value)), size = 2.5) +
                                        ylab("Selection setting") +
                                        xlab("Number of MAGs") +
                                        theme_bw() +
                                        scale_fill_continuous(name = "Mean prop.\nsig. KOs", low = "coral1", high = "coral4",
                                                              limits = c(0, 1)) +
                                        facet_wrap(. ~ Approach) +
                                        theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                        ggtitle("Clade-based simulations")



taxa.based_summary <- readRDS(file = "summaries/MAG.based_taxa.based_param.altered_summary.rds")

taxa.based_summary <- taxa.based_summary[, which(colnames(taxa.based_summary) != "focal_func")]

taxa.based_summary_mean <- aggregate(. ~ MAGs * pseudocount * abun_increase,
                                      data = taxa.based_summary, FUN = mean, na.rm = TRUE, na.action = NULL)

taxa.based_summary_mean <- taxa.based_summary_mean[order(taxa.based_summary_mean$pseudocount, taxa.based_summary_mean$abun_increase, decreasing = TRUE), ]


taxa.based_summary_mean$sel_set <- paste("pseudo", as.character(taxa.based_summary_mean$pseudocount), "; ",
                                          "abun*", as.character(taxa.based_summary_mean$abun_increase), sep = "")

taxa.based_summary_mean$sel_set <- factor(taxa.based_summary_mean$sel_set,
                                           levels = rev(taxa.based_summary_mean$sel_set[-which(duplicated(taxa.based_summary_mean$sel_set))]))

taxa.based_summary_mean$MAG_num_char <- as.character(taxa.based_summary_mean$MAGs)
taxa.based_summary_mean$MAG_num_char <- factor(taxa.based_summary_mean$MAG_num_char,
                                                levels = c("100", "250", "500", "1000", "1595"))


taxa.based_summary_mean_melt <- melt(data = taxa.based_summary_mean, id.vars = c("sel_set", "MAG_num_char"))

taxa.based_summary_mean_melt_sig_prop <- taxa.based_summary_mean_melt[grep("0.05$", taxa.based_summary_mean_melt$variable), ]

taxa.based_summary_mean_melt_sig_prop$variable <- as.character(taxa.based_summary_mean_melt_sig_prop$variable)
taxa.based_summary_mean_melt_sig_prop$Approach <- NA
taxa.based_summary_mean_melt_sig_prop$Approach[which(taxa.based_summary_mean_melt_sig_prop$variable == "POMS_sig_0.05")] <- "POMS"
taxa.based_summary_mean_melt_sig_prop$Approach[which(taxa.based_summary_mean_melt_sig_prop$variable == "regress_specificity_sig_0.05")] <- "Phylo. regression (specificity)"
taxa.based_summary_mean_melt_sig_prop$Approach[which(taxa.based_summary_mean_melt_sig_prop$variable == "regress_sig.taxa_sig_0.05")] <- "Phylo. regression (sig. taxa)"
taxa.based_summary_mean_melt_sig_prop$Approach[which(taxa.based_summary_mean_melt_sig_prop$variable == "wilcoxon.musicc_0.05")] <- "Wilcoxon test (corrected)"


taxa.based_sig_prop_heatmaps <- ggplot(data = taxa.based_summary_mean_melt_sig_prop, aes(y = sel_set, x = MAG_num_char, fill = value)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%0.2f", value)), size = 2.5) +
        ylab("Selection setting") +
        xlab("Number of MAGs") +
        theme_bw() +
        scale_fill_continuous(name = "Mean prop.\nsig. KOs", low = "coral1", high = "coral4",
                              limits = c(0, 1)) +
        facet_wrap(. ~ Approach) +
        theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
        ggtitle("Random taxa simulations")



taxa.based_BSNs_subset <- taxa.based_summary_mean[, c("sel_set", "MAG_num_char", "POMS_num_BSNs")]
taxa.based_BSNs_subset$Approach <- "Random taxa"

clade.based_BSNs_subset <- clade.based_summary_mean[, c("sel_set", "MAG_num_char", "POMS_num_BSNs")]
clade.based_BSNs_subset$Approach <- "Clade-based"

func.based_BSNs_subset <- func.based_summary_mean[, c("sel_set", "MAG_num_char", "POMS_num_BSNs")]
func.based_BSNs_subset$Approach <- "Focal function"

num_BSNs_combined_table <- do.call(rbind, list(taxa.based_BSNs_subset, clade.based_BSNs_subset, func.based_BSNs_subset))

num_BSNs_combined_table$Approach <- factor(num_BSNs_combined_table$Approach,
                                           levels = c("Random taxa", "Clade-based", "Focal function"))

num_BSNs_heatmaps <- ggplot(data = num_BSNs_combined_table, aes(y = sel_set, x = MAG_num_char, fill = POMS_num_BSNs)) +
                                geom_tile() +
                                geom_text(aes(label = sprintf("%0.1f", POMS_num_BSNs))) +
                                facet_wrap(. ~ Approach) +
                                ylab("Selection setting") +
                                xlab("Number of metagenome-assembled genomes") +
                                theme_bw() +
                                scale_fill_continuous(name = "Mean\nnumber of\nsignificant\nnodes",
                                                      low = "light green", high = "dark green",
                                                      limits = c(0, 70))


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_param_altered_num_sig_nodes.pdf",
       plot = num_BSNs_heatmaps,
       device = "pdf",
       width = 10,
       height = 4,
       dpi = 600)

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_param_altered_num_sig_nodes.png",
       plot = num_BSNs_heatmaps,
       device = "png",
       width = 10,
       height = 4,
       dpi = 300)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking_param_altered_heatmaps.pdf",
       plot = func.based_rank_heatmaps,
       device = "pdf",
       width = 8,
       height = 5,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking_param_altered_heatmaps.png",
       plot = func.based_rank_heatmaps,
       device = "png",
       width = 8,
       height = 5,
       dpi = 300)



prop_sig_heatmaps <- plot_grid(taxa.based_sig_prop_heatmaps, clade.based_sig_prop_heatmaps, func.based_sig_prop_heatmaps,
                               nrow = 3, labels = c('a', 'b', 'c'))

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_sig_prop_param_altered_heatmaps.pdf",
       plot = prop_sig_heatmaps,
       device = "pdf",
       width = 7,
       height = 12,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_sig_prop_param_altered_heatmaps.png",
       plot = prop_sig_heatmaps,
       device = "png",
       width = 7,
       height = 12,
       dpi = 300)

