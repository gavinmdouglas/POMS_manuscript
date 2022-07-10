rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

ptm <- proc.time()

func_alt.tools <- mclapply(1:693, function(rep_i) {
  
  prepped_clade.based_sim_info <- readRDS(paste("MAG.based_prepped_clade.based_sim_info_sel1.5/clade.based_sim_info_rep", as.character(rep_i), ".rds", sep = ""))

  prepped_func_abun <- readRDS(paste("func_abun_tables_clade.based_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))

  alt_tools_out_func <- run_alt.tools(func_abun_table = prepped_func_abun,
                                      group1_samples = prepped_clade.based_sim_info$group1,
                                      group2_samples = prepped_clade.based_sim_info$group2,
                                      tools_to_run = c("aldex2", "deseq2"))
  
  out_rds_func <- paste("MAG.based_clade.based_sim_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = "")
  
  saveRDS(object = alt_tools_out_func, file = out_rds_func)
  
  return("Success")
  
}, mc.cores = 10)

print(proc.time() - ptm)

