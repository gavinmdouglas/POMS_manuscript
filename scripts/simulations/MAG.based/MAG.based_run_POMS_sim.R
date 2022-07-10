rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

# Read in needed prepped files.
almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")


ptm <- proc.time()

clade.based_sim_POMS <-  mclapply(X = 1:693, FUN = function(rep_i) {
 
  clade.based_info <- readRDS(paste("MAG.based_prepped_clade.based_sim_info_sel1.5/", "clade.based_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  output <- POMS_pipeline(abun = clade.based_info$taxa_perturb_abun,
                          func = almeida_func_subset,
                          tree = almeida_tree_subset,
                          group1_samples = clade.based_info$group1,
                          group2_samples = clade.based_info$group2,
                          ncores = 1,
                          BSN_p_cutoff = 0.05,
                          BSN_correction = "none",
                          FSN_p_cutoff = 0.05,
                          FSN_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          multinomial_correction = "BH",
                          detailed_output = FALSE,
                          verbose = FALSE)
  
  return(list(output = output))
   
  }, mc.cores = 50)

print(proc.time() - ptm)

saveRDS(object = clade.based_sim_POMS, file = "POMS_out/MAG.based_POMS_sim_rand_clade.based_693reps.rds")
rm(clade.based_sim_POMS)


ptm <- proc.time()

func_sim_POMS <- mclapply(X = 1:1000, FUN = function(rep_i) {

  func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  output <- POMS_pipeline(abun = func_sim_info$taxa_perturb_abun,
                          func = almeida_func_subset,
                          tree = almeida_tree_subset,
                          group1_samples = func_sim_info$group1,
                          group2_samples = func_sim_info$group2,
                          ncores = 1,
                          BSN_p_cutoff = 0.05,
                          BSN_correction = "none",
                          FSN_p_cutoff = 0.05,
                          FSN_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          multinomial_correction = "BH",
                          detailed_output = FALSE,
                          verbose = FALSE)

    return(list(func = func_sim_info$func,
                output = output))
  }
  , mc.cores = 20)


print(proc.time() - ptm)

saveRDS(object = func_sim_POMS, file = "POMS_out/MAG.based_POMS_sim_rand_func_1000reps.rds")

# Remove large files no longer needed.
rm(func_sim_POMS)



ptm <- proc.time()

taxa_sim_POMS <- mclapply(X = 1:1000, FUN = function(rep_i) {

  taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  output <- POMS_pipeline(abun = taxa_sim_info$taxa_perturb_abun,
                          func = almeida_func_subset,
                          tree = almeida_tree_subset,
                          group1_samples = taxa_sim_info$group1,
                          group2_samples = taxa_sim_info$group2,
                          ncores = 1,
                          BSN_p_cutoff = 0.05,
                          BSN_correction = "none",
                          FSN_p_cutoff = 0.05,
                          FSN_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          multinomial_correction = "BH",
                          detailed_output = FALSE,
                          verbose = FALSE)

  return(list(func = taxa_sim_info$func,
              output = output))
}
, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = taxa_sim_POMS, file = "POMS_out/MAG.based_POMS_sim_rand_taxa_1000reps.rds")

# Remove POMS output object.
rm(taxa_sim_POMS)


# Finally, run unperturbed simulations:
almeida_abun_subset <- readRDS("MAG.based_prepped_abun.rds")

ptm <- proc.time()

unperturbed_sim_POMS <- mclapply(X = 1:1000, FUN = function(rep_i) {
  
  taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  output <- POMS_pipeline(abun = almeida_abun_subset,
                          func = almeida_func_subset,
                          tree = almeida_tree_subset,
                          group1_samples = taxa_sim_info$group1,
                          group2_samples = taxa_sim_info$group2,
                          ncores = 1,
                          BSN_p_cutoff = 0.05,
                          BSN_correction = "none",
                          FSN_p_cutoff = 0.05,
                          FSN_correction = "none",
                          min_func_instances = 0,
                          min_func_prop = 0,
                          multinomial_correction = "BH",
                          detailed_output = FALSE,
                          verbose = FALSE)
  
  return(list(output = output))
}
, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = unperturbed_sim_POMS, file = "POMS_out/MAG.based_POMS_sim_unperturbed_1000reps.rds")
