rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

unperturbed_summary <- readRDS(file = "simulation_summaries/unperturbed_summary.rds")

unperturbed_summary_prop_sig <- unperturbed_summary[grep("_0.05$", colnames(unperturbed_summary))]

sapply(unperturbed_summary_prop_sig, mean, na.rm = TRUE)
sapply(unperturbed_summary_prop_sig, sd, na.rm = TRUE)


# The rest of this script is just sanity checks that the same values are computed by looking
# at raw files and using quick and dirty parsing approach:
POMS_raw_out <- readRDS("POMS_out/MAG.based_POMS_sim_unperturbed_1000reps.rds")
unperturbed_regress_raw_out <- readRDS("regression_out/MAG.based_regression_sim_unperturbed_1000reps.rds")
other_alt_raw_out <- readRDS("DA_tool_out/MAG.based_wilcoxon.musicc_wilcoxon.relab_limma.voom_unperturbed_1000reps.rds")
POMS_some_sig <- 0
regress_sig.taxa_rep_some_sig <- 0
regress_specificity_rep_some_sig <- 0
aldex2_some_sig <- 0
deseq2_some_sig <- 0
limmavoom_some_sig <- 0
wilcox.relab_some_sig <- 0
wilcox.musicc_some_sig <- 0

for (i in 1:1000) {

 if (length(which(unperturbed_regress_raw_out[[i]]$sig_taxa_regress$BH < 0.05)) > 0) {
         regress_sig.taxa_rep_some_sig <- regress_sig.taxa_rep_some_sig + 1
 }

  if (length(which(unperturbed_regress_raw_out[[i]]$specificity_regress$BH < 0.05)) > 0) {
          regress_specificity_rep_some_sig <- regress_specificity_rep_some_sig + 1
  }
        
  if (length(which(POMS_raw_out[[i]]$output$results$multinomial_corr < 0.05)) > 0) {
          POMS_some_sig <- POMS_some_sig + 1
  }
        
  if (length(which(other_alt_raw_out[[i]]$limma.voom$BH_corr_p < 0.05)) > 0) {
          limmavoom_some_sig <- limmavoom_some_sig + 1
  }
        
        if (length(which(other_alt_raw_out[[i]]$wilcoxon.relab$BH_corr_p < 0.05)) > 0) {
                wilcox.relab_some_sig <- wilcox.relab_some_sig + 1
        }
        
        if (length(which(other_alt_raw_out[[i]]$wilcoxon.musicc$BH_corr_p < 0.05)) > 0) {
                wilcox.musicc_some_sig <- wilcox.musicc_some_sig + 1
        }
        
        aldex2_deseq2_unperturbed_rep <- readRDS(paste("DA_tool_out/aldex2_deseq2/MAG.based_unperturbed_aldex2_deseq2_rds/rep", as.character(i), ".rds", sep = ""))
        
        if (length(which(aldex2_deseq2_unperturbed_rep$aldex2$BH_corr_p < 0.05)) > 0) {
                aldex2_some_sig <- aldex2_some_sig + 1
        }
        
        if (length(which(aldex2_deseq2_unperturbed_rep$deseq2$BH_corr_p < 0.05)) > 0) {
                deseq2_some_sig <- deseq2_some_sig + 1
        }
}

POMS_some_sig
regress_sig.taxa_rep_some_sig
regress_specificity_rep_some_sig
limmavoom_some_sig
wilcox.relab_some_sig
wilcox.musicc_some_sig
aldex2_some_sig
deseq2_some_sig
