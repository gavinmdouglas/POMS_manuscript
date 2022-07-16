rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

func.based_summary <- readRDS(file = "simulation_summaries/func.based_summary.rds")
func.based_summary_prop.sig <- func.based_summary[grep("sig_0.05", colnames(func.based_summary))]

sapply(func.based_summary_prop.sig, mean)
sapply(func.based_summary_prop.sig, sd)


taxa.based_summary <- readRDS(file = "simulation_summaries/taxa.based_summary.rds")
taxa.based_summary_prop.sig <- taxa.based_summary[grep("sig_0.05", colnames(taxa.based_summary))]

sapply(taxa.based_summary_prop.sig, mean)
sapply(taxa.based_summary_prop.sig, sd)


clade.based_summary <- readRDS(file = "simulation_summaries/clade.based_summary.rds")
clade.based_summary_prop.sig <- clade.based_summary[grep("_0.05$", colnames(clade.based_summary))]

sapply(clade.based_summary_prop.sig, mean)
sapply(clade.based_summary_prop.sig, sd)


# Func.based percent increase:
((sapply(func.based_summary_prop.sig, mean) - sapply(taxa.based_summary_prop.sig, mean)) / sapply(taxa.based_summary_prop.sig, mean)) * 100

((mean(func.based_summary_prop.sig$POMS_sig_0.05) - mean(clade.based_summary_prop.sig$POMS_sig_0.05)) / mean(clade.based_summary_prop.sig$POMS_sig_0.05)) * 100

