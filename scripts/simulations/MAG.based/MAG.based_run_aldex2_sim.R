rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(parallel)


almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")
random_groups <- readRDS("almeida_healthy_random_groups.rds")

almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- POMS::subset_by_col_and_filt(in_tab = almeida_abun,
                                             col2keep = c(random_groups$group1[[1]],
                                                          random_groups$group2[[1]]))

almeida_func_abun <- calc_func_abun_crossprod(in_abun = almeida_abun,
                                              in_func = almeida_func_subset)

unperturbed_alt.tools <- mclapply(1:1000, function(rep_i) {
  
  alt_tools_out <- run_alt.tools(func_abun_table = almeida_func_abun,
                                 group1_samples = random_groups$group1[[rep_i]],
                                 group2_samples = random_groups$group2[[rep_i]],
                                 tools_to_run = c("aldex2"))
  
  alt_tools_out[["func"]] <- "unperturbed"
  
  saveRDS(object = alt_tools_out,
          file = paste("DA_tool_out/aldex2_output/unperturbed/MAG.based_unperturbed_aldex2_", as.character(rep_i), "reps.rds", sep = ""))
  
  return("Success")
  
}, mc.cores = 30)

print(proc.time() - ptm)



ptm <- proc.time()

func_sim_alt.tools <- mclapply(X = 1:1000, FUN = function(rep_i) {

  prepped_func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))

  prepped_func_abun <- readRDS(paste("func_abun_tables_func_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))

  output <- run_alt.tools(func_abun_table = prepped_func_abun,
                                       group1_samples = prepped_func_sim_info$group1,
                                       group2_samples = prepped_func_sim_info$group2,
                                       tools_to_run = "aldex2")

  output$func <- prepped_func_sim_info$func
  
  saveRDS(object = output,
          file = paste("DA_tool_out/aldex2_output/func.based/MAG.based_func.based_aldex2_", as.character(rep_i), "reps.rds", sep = ""))
  
  return("Success")
}
, mc.cores = 30)

print(proc.time() - ptm)




ptm <- proc.time()

taxa_sim_alt.tools <- mclapply(X = 1:1000, FUN = function(rep_i) {
  
  prepped_taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  prepped_func_abun <- readRDS(paste("func_abun_tables_taxa_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
  output <- run_alt.tools(func_abun_table = prepped_func_abun,
                                       group1_samples = prepped_taxa_sim_info$group1,
                                       group2_samples = prepped_taxa_sim_info$group2,
                                       tools_to_run = "aldex2")
  
  output$func <- prepped_taxa_sim_info$func
  
  saveRDS(object = output,
          file = paste("DA_tool_out/aldex2_output/taxa.based/MAG.based_taxa.based_aldex2_", as.character(rep_i), "reps.rds", sep = ""))
  
  return("Success")
}
, mc.cores = 30)

print(proc.time() - ptm)



ptm <- proc.time()

clade.based_alt.tools <- mclapply(X = 1:693, FUN = function(rep_i) {
  
  prepped_taxa_sim_info <- readRDS(paste("MAG.based_prepped_clade.based_sim_info_sel1.5/clade.based_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  prepped_func_abun <- readRDS(paste("func_abun_tables_clade.based_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
  output <- run_alt.tools(func_abun_table = prepped_func_abun,
                          group1_samples = prepped_taxa_sim_info$group1,
                          group2_samples = prepped_taxa_sim_info$group2,
                          tools_to_run = "aldex2")
  
  output$func <- prepped_taxa_sim_info$func
  
  saveRDS(object = output,
          file = paste("DA_tool_out/aldex2_output/clade.based/MAG.based_clade.based_aldex2_", as.character(rep_i), "reps.rds", sep = ""))
  
  return("Success")
}
, mc.cores = 30)

print(proc.time() - ptm)
