rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

# Read in needed prepped files.
almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")

func_sim_info <- readRDS("MAG.based_prepped_func_sim_info.rds")

ptm <- proc.time()

func_sim_POMS <- mclapply(X = 1:1000, FUN = function(rep_i) {

  output <- POMS_pipeline(abun = func_sim_info[[rep_i]]$taxa_perturb_abun,
                          func = almeida_func_subset,
                          phylogeny = almeida_tree_subset,
                          group1_samples = func_sim_info[[rep_i]]$group1,
                          group2_samples = func_sim_info[[rep_i]]$group2,
                          ncores = 1,
                          balance_p_cutoff = 0.05,
                          balance_correction = "none",
                          function_p_cutoff = 0.05,
                          function_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          run_multinomial_test = TRUE,
                          multinomial_correction = "BH",
                          calc_node_dist = FALSE,
                          detailed_output = FALSE,
                          verbose = FALSE)

    return(list(func = func_sim_info[[rep_i]]$func,
                output = output))
  }
  , mc.cores = 20)


print(proc.time() - ptm)

saveRDS(object = func_sim_POMS, file = "MAG.based_POMS_sim_rand_func_1000reps.rds")

# Remove large files no longer needed.
rm(func_sim_info)
rm(func_sim_POMS)


taxa_sim_info <- readRDS("MAG.based_prepped_rand.taxa_sim_info.rds")

ptm <- proc.time()

taxa_sim_POMS <- mclapply(X = 1:1000, FUN = function(rep_i) {

  output <- POMS_pipeline(abun = taxa_sim_info[[rep_i]]$taxa_perturb_abun,
                          func = almeida_func_subset,
                          phylogeny = almeida_tree_subset,
                          group1_samples = taxa_sim_info[[rep_i]]$group1,
                          group2_samples = taxa_sim_info[[rep_i]]$group2,
                          ncores = 1,
                          balance_p_cutoff = 0.05,
                          balance_correction = "none",
                          function_p_cutoff = 0.05,
                          function_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          run_multinomial_test = TRUE,
                          multinomial_correction = "BH",
                          calc_node_dist = FALSE,
                          detailed_output = FALSE,
                          verbose = FALSE)

  return(list(func = taxa_sim_info[[rep_i]]$func,
              output = output))
}
, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = taxa_sim_POMS, file = "MAG.based_POMS_sim_rand_taxa_1000reps.rds")

# Remove POMS output object.
rm(taxa_sim_POMS)


# Finally, run unperturbed simulations:
almeida_abun_subset <- readRDS("MAG.based_prepped_abun.rds")

ptm <- proc.time()

unperturbed_sim_POMS <- mclapply(X = 1:1000, FUN = function(rep_i) {
  
  output <- POMS_pipeline(abun = almeida_abun_subset,
                          func = almeida_func_subset,
                          phylogeny = almeida_tree_subset,
                          group1_samples = taxa_sim_info[[rep_i]]$group1,
                          group2_samples = taxa_sim_info[[rep_i]]$group2,
                          ncores = 1,
                          balance_p_cutoff = 0.05,
                          balance_correction = "none",
                          function_p_cutoff = 0.05,
                          function_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          run_multinomial_test = TRUE,
                          multinomial_correction = "BH",
                          calc_node_dist = FALSE,
                          detailed_output = FALSE,
                          verbose = FALSE)
  
  return(list(output = output))
}
, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = unperturbed_sim_POMS, file = "MAG.based_POMS_sim_unperturbed_1000reps.rds")
