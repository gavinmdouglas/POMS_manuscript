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


run_alt.tools_exceptALDEX2 <- function(func_abun_table, group1_samples, group2_samples, USCGs) {
  
  func_abun_table_ceil <- ceiling(func_abun_table)
  
  DA_out <- list()
  
  # DA_out[["aldex2"]] <- run_2group_ALDEx2(in_table = func_abun_table_ceil,
  #                                         group1_samples = group1_samples,
  #                                         group2_samples = group2_samples)
  
  DA_out[["deseq2"]] <- deseq2_default_two_groups(table = func_abun_table_ceil,
                                                  group1 = group1_samples,
                                                  group2 = group2_samples)
  
  DA_out[["limma.voom"]] <- limma_voom_two_group_TMM(table = func_abun_table,
                                                     group1 = group1_samples,
                                                     group2 = group2_samples)
  
  
  DA_out[["wilcoxon.relab"]] <- wilcoxon_2group_pvalues(intable = func_abun_table,
                                                        group1_samples = group1_samples,
                                                        group2_samples = group2_samples,
                                                        convert_relab = TRUE)
  
  dataset_uscg_set <- USCGs[which(USCGs %in% rownames(func_abun_table))]
  func_abun_table_musicc <- data.frame(sweep(func_abun_table, 2, colMedians(as.matrix(func_abun_table[dataset_uscg_set, ])), `/`))
  
  DA_out[["wilcoxon.musicc"]] <- wilcoxon_2group_pvalues(intable = func_abun_table_musicc,
                                                         group1_samples = group1_samples,
                                                         group2_samples = group2_samples,
                                                         convert_relab = FALSE)
  
  return(DA_out)
  
}

ptm <- proc.time()

func_sim_alt.tools <- mclapply(X = 1:100, FUN = function(rep_i) {
  
  rep_func_abun <- calc_func_abun(in_abun = func_sim_info[[rep_i]]$taxa_perturb_abun, in_func = almeida_func_subset, ncores = 1)

  output <- run_alt.tools_exceptALDEX2(func_abun_table = rep_func_abun,
                                        group1_samples = func_sim_info[[rep_i]]$group1,
                                        group2_samples = func_sim_info[[rep_i]]$group2,
                                        USCGs = musicc_uscgs)
  return(list(func = func_sim_info[[rep_i]]$func,
              output = output))
}
, mc.cores = 30)


print(proc.time() - ptm)

saveRDS(object = func_sim_alt.tools, file = "MAG.based_alt.tools_sim_rand_func_100reps.rds")

rm(func_sim_info)
rm(func_sim_alt.tools)

taxa_sim_info <- readRDS("MAG.based_prepped_rand.taxa_sim_info.rds")

ptm <- proc.time()

taxa_sim_alt.tools <- mclapply(X = 1:100, FUN = function(rep_i) {
  
  rep_func_abun <- calc_func_abun(in_abun = taxa_sim_info[[rep_i]]$taxa_perturb_abun,
                                  in_func = almeida_func_subset, ncores = 1)
  
  output <- run_alt.tools_exceptALDEX2(func_abun_table = rep_func_abun,
                          group1_samples = taxa_sim_info[[rep_i]]$group1,
                          group2_samples = taxa_sim_info[[rep_i]]$group2,
                          USCGs = musicc_uscgs)
  
  return(list(func = taxa_sim_info[[rep_i]]$func,
              output = output))
}
, mc.cores = 30)

print(proc.time() - ptm)

saveRDS(object = taxa_sim_alt.tools, file = "MAG.based_alt.tools_sim_rand_taxa_100reps.rds")

rm(taxa_sim_alt.tools)


# Finally, run unperturbed simulations:
almeida_abun_subset <- readRDS("MAG.based_prepped_abun.rds")

ptm <- proc.time()

unperturbed_sim_alt.tools <- mclapply(X = 1:100, FUN = function(rep_i) {
  
  rep_func_abun <- calc_func_abun(in_abun = almeida_abun_subset,
                                  in_func = almeida_func_subset, ncores = 1)
  
  output <- run_alt.tools_exceptALDEX2(func_abun_table = rep_func_abun,
                          group1_samples = taxa_sim_info[[rep_i]]$group1,
                          group2_samples = taxa_sim_info[[rep_i]]$group2,
                          USCGs = musicc_uscgs)
  
  return(list(output = output))
}
, mc.cores = 30)
  
print(proc.time() - ptm)

saveRDS(object = unperturbed_sim_alt.tools, file = "MAG.based_alt.tools_sim_unperturbed_100reps.rds")
