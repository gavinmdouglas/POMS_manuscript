### Compare results identified when jittering an input table when
### run through basic Wilcoxon test on functions and with POMS.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/reference_genome_sim/")

source("/home/gavin/github_repos/POMS/balance_tree_functions.R")

library(cowplot)
library(ggExtra)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)


ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep="\t", header=TRUE, row.names=1)

func_sim_info_100rep <- readRDS("func_sim_info_100rep.rds")

wilcox_out <- readRDS("func_rand_output_wilcoxon_100reps_BEZI_cutoffs.rds")


BEZI_tables <- list()

BEZI_tables[["mu0.1_nu_0.5"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.5_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.9"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.9_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.99_sigma1.tsv",
                                             header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.01_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.01_nu_0.99_sigma1.tsv",
                                              header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

cutoff_func <- list()
wilcox_out_BY <- list()

for(cutoff in names(BEZI_tables)) {
  
  wilcox_out_BY[[cutoff]] <- lapply(wilcox_out[[cutoff]], function(x) { p.adjust(x, "BY") })
  
  cutoff_func[[cutoff]] <- ref_func[rownames(BEZI_tables[[cutoff]]), ]
  
  if(length(which(colSums(cutoff_func[[cutoff]]) == 0)) > 0) {
    cutoff_func[[cutoff]] <- cutoff_func[[cutoff]][, -which(colSums(cutoff_func[[cutoff]]) == 0)]
  }
  
}

wilcox_ranks <- list()

for(cutoff in names(BEZI_tables)) {

  wilcoxon_ranks <- c()
  func_ids <- c()
  num_sig <- c()
  
  for(rep in c(1:100)) {
    
    f <- func_sim_info_100rep[[cutoff]][[rep]]$func
    
    func_ids <- c(func_ids, f)
    
    f_i <- which(colnames(cutoff_func[[cutoff]]) == f)
    
    num_pos_mags <- length(which(cutoff_func[[cutoff]][, f] > 0))
    
    if(wilcox_out_BY[[cutoff]][[rep]][f_i] < 0.05) {
      focal_rank <- rank(wilcox_out_BY[[cutoff]][[rep]])[f_i]
    } else {
      focal_rank <- NA
    }
    
    num_sig <- c(num_sig, length(which(wilcox_out_BY[[cutoff]][[rep]] < 0.05)))
    
    wilcoxon_ranks <- c(wilcoxon_ranks, focal_rank)
     
  }
  
  wilcoxon_rel_rank <- (wilcoxon_ranks / num_sig) * 100

  wilcox_ranks[[cutoff]] <- list(func_ids=func_ids, wilcoxon_ranks=wilcoxon_ranks, wilcoxon_rel_rank=wilcoxon_rel_rank, num_sig=num_sig)
  
}

wilcox_ranks_df <- data.frame(mu0.1_nu0.5=wilcox_ranks$mu0.1_nu_0.5$wilcoxon_rel_rank,
                              mu0.1_nu0.9=wilcox_ranks$mu0.1_nu_0.9$wilcoxon_rel_rank,
                              mu0.1_nu0.99=wilcox_ranks$mu0.1_nu_0.99$wilcoxon_rel_rank,
                              mu0.01_nu0.99=wilcox_ranks$mu0.01_nu_0.99$wilcoxon_rel_rank)

wilcox_ranks_df_melt <- melt(wilcox_ranks_df)

wilcox_rank_plot <- ggplot(data=wilcox_ranks_df_melt, aes(x=variable, y=value)) +
                              geom_boxplot(outlier.shape = NA) +
                              geom_jitter(cex=1, width=0.2) +
                              ylab("Rel. ranking of focal gene in Wilcoxon test") +
                              xlab("Genome rel. abun. simulation setting") +
                              ylim(c(0, 100)) +
                              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                    legend.position = "none",
                                    panel.background = element_rect(fill = 'grey85'))


taxa_prevalence <- data.frame(mu0.1_nu0.5=c(as.numeric(rowSums(BEZI_tables$mu0.1_nu_0.5 > 0) / ncol(BEZI_tables$mu0.1_nu_0.5)), rep(NA, 3000 - nrow(BEZI_tables$mu0.1_nu_0.5))),
                              mu0.1_nu0.9=c(as.numeric(rowSums(BEZI_tables$mu0.1_nu_0.9 > 0) / ncol(BEZI_tables$mu0.1_nu_0.9)), rep(NA, 3000 - nrow(BEZI_tables$mu0.1_nu_0.9))),
                              mu0.1_nu0.99=c(as.numeric(rowSums(BEZI_tables$mu0.1_nu_0.99 > 0) / ncol(BEZI_tables$mu0.1_nu_0.99)), rep(NA, 3000 - nrow(BEZI_tables$mu0.1_nu_0.99))),
                              mu0.01_nu0.99=c(as.numeric(rowSums(BEZI_tables$mu0.01_nu_0.99 > 0) / ncol(BEZI_tables$mu0.01_nu_0.99)), rep(NA, 3000 - nrow(BEZI_tables$mu0.01_nu_0.99))))

taxa_prevalence <- taxa_prevalence * 100

taxa_prevalence_melt <- melt(taxa_prevalence)

prevalence_plot <- ggplot(data=taxa_prevalence_melt, aes(x=variable, y=value)) +
                            geom_boxplot(outlier.shape = NA) +
                            geom_jitter(cex=1, width=0.2) +
                            ylab("Genome prevalence across samples (%)") +
                            xlab("Genome rel. abun. simulation setting") +
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                  legend.position = "none",
                                  panel.background = element_rect(fill = 'grey85'))






taxa_richness <- data.frame(mu0.1_nu0.5=c(as.numeric(colSums(BEZI_tables$mu0.1_nu_0.5 > 0) / nrow(BEZI_tables$mu0.1_nu_0.5)), rep(NA, 1000 - ncol(BEZI_tables$mu0.1_nu_0.5))),
                              mu0.1_nu0.9=c(as.numeric(colSums(BEZI_tables$mu0.1_nu_0.9 > 0) / nrow(BEZI_tables$mu0.1_nu_0.9)), rep(NA, 1000 - ncol(BEZI_tables$mu0.1_nu_0.9))),
                              mu0.1_nu0.99=c(as.numeric(colSums(BEZI_tables$mu0.1_nu_0.99 > 0) / nrow(BEZI_tables$mu0.1_nu_0.99)), rep(NA, 1000 - ncol(BEZI_tables$mu0.1_nu_0.99))),
                              mu0.01_nu0.99=c(as.numeric(colSums(BEZI_tables$mu0.01_nu_0.99 > 0) / nrow(BEZI_tables$mu0.01_nu_0.99)), rep(NA, 1000 - ncol(BEZI_tables$mu0.01_nu_0.99))))

taxa_richness <- taxa_richness * 100

taxa_richness_melt <- melt(taxa_richness)

richness_plot <- ggplot(data=taxa_richness_melt, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(cex=1, width=0.2) +
  ylab("Genome richness across samples (%)") +
  xlab("Genome rel. abun. simulation setting") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        panel.background = element_rect(fill = 'grey85'))




plot_grid(wilcox_rank_plot, prevalence_plot)