rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

ptm <- proc.time()

func_sim_alt.tools <- list()

for (rep_i in 1:10) {
  
  print(rep_i)
  
  prepped_func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))

  prepped_func_abun <- readRDS(paste("func_abun_tables_func_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))

  alt_tools_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                 group1_samples = prepped_func_sim_info$group1,
                                 group2_samples = prepped_func_sim_info$group2,
                                 tools_to_run = c("aldex2", "deseq2"))

  alt_tools_out[["func"]] <- prepped_func_sim_info$func
  
  func_sim_alt.tools[[rep_i]] <- alt_tools_out
}

saveRDS(object = func_sim_alt.tools, file = "MAG.based_aldex2_deseq2_sim_rand_func_10reps.rds")

print(proc.time() - ptm)




ptm <- proc.time()

taxa_sim_alt.tools <- list()

for (rep_i in 1:10) {
  
  print(rep_i)
  
  prepped_taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  prepped_func_abun <- readRDS(paste("func_abun_tables_taxa_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
  alt_tools_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                 group1_samples = prepped_taxa_sim_info$group1,
                                 group2_samples = prepped_taxa_sim_info$group2,
                                 tools_to_run = c("aldex2", "deseq2"))
  
  alt_tools_out[["func"]] <- prepped_taxa_sim_info$func
  
  taxa_sim_alt.tools[[rep_i]] <- alt_tools_out
}

saveRDS(object = taxa_sim_alt.tools, file = "MAG.based_aldex2_deseq2_sim_rand_taxa_10reps.rds")

print(proc.time() - ptm)
