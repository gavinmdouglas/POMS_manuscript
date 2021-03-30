rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

# Read in needed prepped files.
almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")

func_sim_info <- readRDS("MAG.based_prepped_func_sim_info.rds")

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1

ptm <- proc.time()

func_sim_alt.tools <- mclapply(X = 1:100, FUN = function(rep_i) {
  
  rep_func_abun <- calc_func_abun(in_abun = func_sim_info[[rep_i]]$taxa_perturb_abun, in_func = almeida_func_subset, ncores = 1)

  wilcoxon.relab.out <- wilcoxon_2group_pvalues(intable = rep_func_abun,
                                                group1_samples = func_sim_info[[rep_i]]$group1,
                                                group2_samples = func_sim_info[[rep_i]]$group2,
                                                convert_relab = TRUE)
  
  musicc_uscgs_rep <- musicc_uscgs[which(musicc_uscgs %in% rownames(rep_func_abun))]
  func_abun_table_musicc <- data.frame(sweep(rep_func_abun, 2, colMedians(as.matrix(rep_func_abun[musicc_uscgs_rep, ])), `/`))
  
  wilcoxon.musicc.out <- wilcoxon_2group_pvalues(intable = func_abun_table_musicc,
                                                 group1_samples = func_sim_info[[rep_i]]$group1,
                                                 group2_samples = func_sim_info[[rep_i]]$group2,
                                                 convert_relab = FALSE)
  
  return(list(func = func_sim_info[[rep_i]]$func,
              wilcoxon.relab.out = wilcoxon.relab.out,
              wilcoxon.musicc.out = wilcoxon.musicc.out))
}
, mc.cores = 30)


print(proc.time() - ptm)

saveRDS(object = func_sim_alt.tools, file = "MAG.based_Wilcoxon_sim_rand_func_100reps.rds")

rm(func_sim_info)
rm(func_sim_alt.tools)

taxa_sim_info <- readRDS("MAG.based_prepped_rand.taxa_sim_info.rds")

ptm <- proc.time()

taxa_sim_alt.tools <- mclapply(X = 1:100, FUN = function(rep_i) {
  
  rep_func_abun <- calc_func_abun(in_abun = taxa_sim_info[[rep_i]]$taxa_perturb_abun,
                                  in_func = almeida_func_subset, ncores = 1)
  
  wilcoxon.relab.out <- wilcoxon_2group_pvalues(intable = rep_func_abun,
                                                group1_samples = taxa_sim_info[[rep_i]]$group1,
                                                group2_samples = taxa_sim_info[[rep_i]]$group2,
                                                convert_relab = TRUE)
  
  musicc_uscgs_rep <- musicc_uscgs[which(musicc_uscgs %in% rownames(rep_func_abun))]
  func_abun_table_musicc <- data.frame(sweep(rep_func_abun, 2, colMedians(as.matrix(rep_func_abun[musicc_uscgs_rep, ])), `/`))
  
  wilcoxon.musicc.out <- wilcoxon_2group_pvalues(intable = func_abun_table_musicc,
                                                 group1_samples = taxa_sim_info[[rep_i]]$group1,
                                                 group2_samples = taxa_sim_info[[rep_i]]$group2,
                                                 convert_relab = FALSE)
  
  return(list(func = taxa_sim_info[[rep_i]]$func,
              wilcoxon.relab.out = wilcoxon.relab.out,
              wilcoxon.musicc.out = wilcoxon.musicc.out))
}
, mc.cores = 30)

print(proc.time() - ptm)

saveRDS(object = taxa_sim_alt.tools, file = "MAG.based_Wilcoxon_sim_rand_taxa_100reps.rds")

rm(taxa_sim_alt.tools)


# Finally, run unperturbed simulations:
almeida_abun_subset <- readRDS("MAG.based_prepped_abun.rds")

ptm <- proc.time()

unperturbed_sim_alt.tools <- mclapply(X = 1:100, FUN = function(rep_i) {
  
  rep_func_abun <- calc_func_abun(in_abun = almeida_abun_subset,
                                  in_func = almeida_func_subset, ncores = 1)
  
  wilcoxon.relab.out <- wilcoxon_2group_pvalues(intable = rep_func_abun,
                                                group1_samples = taxa_sim_info[[rep_i]]$group1,
                                                group2_samples = taxa_sim_info[[rep_i]]$group2,
                                                convert_relab = TRUE)
  
  musicc_uscgs_rep <- musicc_uscgs[which(musicc_uscgs %in% rownames(rep_func_abun))]
  func_abun_table_musicc <- data.frame(sweep(rep_func_abun, 2, colMedians(as.matrix(rep_func_abun[musicc_uscgs_rep, ])), `/`))
  
  wilcoxon.musicc.out <- wilcoxon_2group_pvalues(intable = func_abun_table_musicc,
                                                 group1_samples = taxa_sim_info[[rep_i]]$group1,
                                                 group2_samples = taxa_sim_info[[rep_i]]$group2,
                                                 convert_relab = FALSE)
  
  return(list(wilcoxon.relab.out = wilcoxon.relab.out,
              wilcoxon.musicc.out = wilcoxon.musicc.out))
}
, mc.cores = 30)
  
print(proc.time() - ptm)

saveRDS(object = unperturbed_sim_alt.tools, file = "MAG.based_Wilcoxon_sim_unperturbed_100reps.rds")
