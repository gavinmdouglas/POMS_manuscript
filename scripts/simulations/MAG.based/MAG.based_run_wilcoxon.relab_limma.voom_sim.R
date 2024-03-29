rm(list = ls(all.names = TRUE))

# Run wilcoxon.relab and limma.voom for 1000 replicates (these are run together because they run very quickly)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

ptm <- proc.time()

func_sim_alt.tools <- mclapply(1:1000, function(rep_i) {

    prepped_func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
    prepped_func_abun <- readRDS(paste("func_abun_tables_func_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
    alt_tools_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                   group1_samples = prepped_func_sim_info$group1,
                                   group2_samples = prepped_func_sim_info$group2,
                                   tools_to_run = c("limma.voom", "wilcoxon.relab"))
  
    alt_tools_out[["func"]] <- prepped_func_sim_info$func
    
    return(alt_tools_out)
  }, mc.cores = 30)


saveRDS(object = func_sim_alt.tools, file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_func_1000reps.rds")

print(proc.time() - ptm)




ptm <- proc.time()

taxa_sim_alt.tools <- mclapply(1:1000, function(rep_i) {
  
  prepped_taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  prepped_func_abun <- readRDS(paste("func_abun_tables_taxa_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
  alt_tools_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                 group1_samples = prepped_taxa_sim_info$group1,
                                 group2_samples = prepped_taxa_sim_info$group2,
                                 tools_to_run = c("limma.voom", "wilcoxon.relab"))
  
  alt_tools_out[["func"]] <- prepped_taxa_sim_info$func
  
  return(alt_tools_out)
}, mc.cores = 30)


saveRDS(object = taxa_sim_alt.tools, file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_taxa_1000reps.rds")

print(proc.time() - ptm)


ptm <- proc.time()

unperturbed_alt.tools <- mclapply(1:1000, function(rep_i) {
  
  prepped_unperturbed_info <- readRDS(paste("MAG.based_prepped_unperturbed_info_sel1.5/unperturbed_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  prepped_func_abun <- readRDS(paste("func_abun_tables_unperturbed_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
  alt_tools_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                 group1_samples = prepped_unperturbed_info$group1,
                                 group2_samples = prepped_unperturbed_info$group2,
                                 tools_to_run = c("limma.voom", "wilcoxon.relab"))
  
  alt_tools_out[["func"]] <- prepped_unperturbed_info$func
  
  return(alt_tools_out)
}, mc.cores = 30)


saveRDS(object = unperturbed_alt.tools, file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_taxa_1000reps.rds")

print(proc.time() - ptm)