rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1

ptm <- proc.time()

clade.based_sim_alt.tools <- mclapply(1:693, function(rep_i) {

    clade.based_info <- readRDS(paste("MAG.based_prepped_clade.based_sim_info_sel1.5/", "clade.based_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
    prepped_func_abun <- readRDS(paste("func_abun_tables_clade.based_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
    alt_tools_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                   group1_samples = clade.based_info$group1,
                                   group2_samples = clade.based_info$group2,
                                   USCGs = musicc_uscgs,
                                   tools_to_run = c("limma.voom", "wilcoxon.relab", "wilcoxon.musicc"))
    
    return(alt_tools_out)
  }, mc.cores = 30)


saveRDS(object = clade.based_sim_alt.tools, file = "MAG.based_wilcoxon.relab_wilcoxon.musicc_limma.voom_sim_rand_clade.based_693reps.rds")

print(proc.time() - ptm)


