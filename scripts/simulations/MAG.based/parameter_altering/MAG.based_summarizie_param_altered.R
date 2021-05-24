rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(pheatmap)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")


# Create dummy df for results based on single replicate.
test_POMS_rds <- list(output = readRDS("POMS_out/POMS_func_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_POMS_rds$func <- test_POMS_rds$output$func
test_POMS_rds_list <- list(test_POMS_rds)
test_wilcoxon_rds <- readRDS("wilcoxon_out/wilcoxon_func_out_rep10_MAGs1000_pseudo0_increase1.05.rds")
test_wilcoxon_rds_list <- list(test_wilcoxon_rds)
ex_func <- readRDS("prepped_func_tables/subset_1000MAGs_func_rep10.rds")
summary_df <- simulation_summaries(POMS_sims = test_POMS_rds_list,
                                   alt_tool_sims = test_wilcoxon_rds_list,
                                   focal_func_present = TRUE,
                                   func_table = ex_func,
                                   num_cores = 20)
summary_df[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- NA
summary_df <- summary_df[-1, ]
        
parameter_settings <- list()

MAG_nums <- c(1595, 1250, 1000, 750, 500, 250, 100, 50)

pseudocount_settings <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

abun_increase_settings <- c(1.5, 1.3, 1.1, 1.05)


for (rep_i in 1:25) {
  for (MAG_num in MAG_nums) {

    POMS_output <- list()
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    for (pseudocount_set in pseudocount_settings) {
      for (abun_increase_set in abun_increase_settings) {
        
        
        POMS_rds_input <- list(output = readRDS(paste("POMS_out/POMS_func_out_",
                                                      "rep", as.character(rep_i),
                                                      "_MAGs", as.character(MAG_num),
                                                      "_pseudo", as.character(pseudocount_set),
                                                      "_increase", as.character(abun_increase_set),
                                                      ".rds", sep = "")))
        
        wilcoxon_rds_input <- list(readRDS(paste("wilcoxon_out/wilcoxon_func_out_",
                                                 "rep", as.character(rep_i),
                                                 "_MAGs", as.character(MAG_num),
                                                 "_pseudo", as.character(pseudocount_set),
                                                 "_increase", as.character(abun_increase_set),
                                                 ".rds", sep = "")))
        
        POMS_output <- list(POMS_rds_input)
        
        POMS_output[[1]]$func <- POMS_rds_input$output$func
  
        parameter_rep_summary <- simulation_summaries(POMS_sims = POMS_output,
                                                      alt_tool_sims = wilcoxon_rds_input,
                                                      focal_func_present = TRUE,
                                                      func_table = MAG_rep_func,
                                                      num_cores = 1)
        
        parameter_rep_summary[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- c(rep_i, MAG_num, pseudocount_set, abun_increase_set)
        
        summary_df <- rbind(summary_df, parameter_rep_summary)
        
      }
    }
  }
}

saveRDS(object = summary_df, file = "POMS_wilcoxon.musicc_out_summary.rds")
