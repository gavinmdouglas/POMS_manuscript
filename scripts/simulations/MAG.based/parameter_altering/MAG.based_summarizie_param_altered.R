rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(pheatmap)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")


# func.based
# Create dummy df for results based on single replicate.
test_func.based_POMS_rds <- list(output = readRDS("POMS_out/func.based/POMS_func.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_func.based_POMS_rds$func <- test_func.based_POMS_rds$output$func
test_func.based_POMS_rds_list <- list(test_func.based_POMS_rds)

test_func.based_regress_rds <- list(output = readRDS("regress_out/func.based/regress_func.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_func.based_regress_rds_list <- list(test_func.based_regress_rds$output)
test_func.based_regress_rds_list[[1]]$func <- test_func.based_POMS_rds$output$func


test_func.based_wilcoxon_rds <- readRDS("wilcoxon_out/func.based/wilcoxon_func.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds")
test_func.based_wilcoxon_rds_list <- list(test_func.based_wilcoxon_rds)

ex_func <- readRDS("prepped_func_tables/subset_1000MAGs_func_rep10.rds")

func.based_summary_df <- simulation_summaries(POMS_sims = test_func.based_POMS_rds_list,
                                              regress_sims = test_func.based_regress_rds_list,
                                              alt_tool_sims = test_func.based_wilcoxon_rds_list,
                                              focal_func_present = TRUE,
                                              func_table = ex_func,
                                              num_cores = 20,
                                              sig_cutoffs = c(0.05))

func.based_summary_df[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- NA
func.based_summary_df <- func.based_summary_df[-1, ]
        
parameter_settings <- list()

MAG_nums <- c(1595, 1000, 500, 250, 100)

pseudocount_settings <- c(0, 0.3, 0.7, 1)

abun_increase_settings <- c(1.5, 1.25, 1.05)

num_reps <- 10

for (rep_i in 1:num_reps) {

  for (MAG_num in MAG_nums) {

    POMS_output <- list()
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    for (pseudocount_set in pseudocount_settings) {
      
      for (abun_increase_set in abun_increase_settings) {
        
        
        POMS_rds_input <- list(output = readRDS(paste("POMS_out/func.based/POMS_func.based_out_",
                                                      "rep", as.character(rep_i),
                                                      "_MAGs", as.character(MAG_num),
                                                      "_pseudo", as.character(pseudocount_set),
                                                      "_increase", as.character(abun_increase_set),
                                                      ".rds", sep = "")))
        
        regress_rds_input <- list(output = readRDS(paste("regress_out/func.based/regress_func.based_out_",
                                                      "rep", as.character(rep_i),
                                                      "_MAGs", as.character(MAG_num),
                                                      "_pseudo", as.character(pseudocount_set),
                                                      "_increase", as.character(abun_increase_set),
                                                      ".rds", sep = "")))
        
        wilcoxon_rds_input <- list(readRDS(paste("wilcoxon_out/func.based/wilcoxon_func.based_out_",
                                                 "rep", as.character(rep_i),
                                                 "_MAGs", as.character(MAG_num),
                                                 "_pseudo", as.character(pseudocount_set),
                                                 "_increase", as.character(abun_increase_set),
                                                 ".rds", sep = "")))
        
        POMS_output <- list(POMS_rds_input)
        POMS_output[[1]]$func <- POMS_rds_input$output$func
  
        regress_output <- list(regress_rds_input$output)
        regress_output[[1]]$func <- POMS_output[[1]]$func
         
        parameter_rep_summary <- simulation_summaries(POMS_sims = POMS_output,
                                                      regress_sims = regress_output,
                                                      alt_tool_sims = wilcoxon_rds_input,
                                                      focal_func_present = TRUE,
                                                      func_table = MAG_rep_func,
                                                      num_cores = 1, sig_cutoffs = c(0.05))
        
        parameter_rep_summary[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- c(rep_i, MAG_num, pseudocount_set, abun_increase_set)
        
        func.based_summary_df <- rbind(func.based_summary_df, parameter_rep_summary)
        
      }
    }
  }
}

saveRDS(object = func.based_summary_df, file = "summaries/MAG.based_func.based_param.altered_summary.rds")





# clade.based

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(pheatmap)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")



# Create dummy df for results based on single replicate.
test_clade.based_POMS_rds <- list(output = readRDS("POMS_out/clade.based/POMS_func_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_clade.based_POMS_rds_list <- list(test_clade.based_POMS_rds)

test_clade.based_regress_rds <- list(output = readRDS("regress_out/clade.based/regress_clade.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_clade.based_regress_rds_list <- list(test_clade.based_regress_rds$output)

test_clade.based_wilcoxon_rds <- readRDS("wilcoxon_out/clade.based/wilcoxon_clade.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds")
test_clade.based_wilcoxon_rds_list <- list(test_clade.based_wilcoxon_rds)

ex_func <- readRDS("prepped_func_tables/subset_1000MAGs_func_rep10.rds")

clade.based_summary_df <- simulation_summaries(POMS_sims = test_clade.based_POMS_rds_list,
                                              regress_sims = test_clade.based_regress_rds_list,
                                              alt_tool_sims = test_clade.based_wilcoxon_rds_list,
                                              focal_func_present = FALSE,
                                              func_table = ex_func,
                                              num_cores = 20,
                                              sig_cutoffs = c(0.05))

clade.based_summary_df[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- NA
clade.based_summary_df <- clade.based_summary_df[-1, ]

parameter_settings <- list()

MAG_nums <- c(1595, 1000, 500, 250, 100)

pseudocount_settings <- c(0, 0.3, 0.7, 1)

abun_increase_settings <- c(1.5, 1.25, 1.05)

num_reps <- 10

for (rep_i in 1:num_reps) {
  
  for (MAG_num in MAG_nums) {
    
    POMS_output <- list()
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    for (pseudocount_set in pseudocount_settings) {
      
      for (abun_increase_set in abun_increase_settings) {
        
        
        POMS_rds_input <- list(output = readRDS(paste("POMS_out/clade.based/POMS_func_out_",
                                                      "rep", as.character(rep_i),
                                                      "_MAGs", as.character(MAG_num),
                                                      "_pseudo", as.character(pseudocount_set),
                                                      "_increase", as.character(abun_increase_set),
                                                      ".rds", sep = "")))
        
        regress_rds_input <- list(output = readRDS(paste("regress_out/clade.based/regress_clade.based_out_",
                                                         "rep", as.character(rep_i),
                                                         "_MAGs", as.character(MAG_num),
                                                         "_pseudo", as.character(pseudocount_set),
                                                         "_increase", as.character(abun_increase_set),
                                                         ".rds", sep = "")))
        
        wilcoxon_rds_input <- list(readRDS(paste("wilcoxon_out/clade.based/wilcoxon_clade.based_out_",
                                                 "rep", as.character(rep_i),
                                                 "_MAGs", as.character(MAG_num),
                                                 "_pseudo", as.character(pseudocount_set),
                                                 "_increase", as.character(abun_increase_set),
                                                 ".rds", sep = "")))
        
        POMS_output <- list(POMS_rds_input)
        
        regress_output <- list(regress_rds_input$output)
        
        parameter_rep_summary <- simulation_summaries(POMS_sims = POMS_output,
                                                      regress_sims = regress_output,
                                                      alt_tool_sims = wilcoxon_rds_input,
                                                      focal_func_present = FALSE,
                                                      func_table = MAG_rep_func,
                                                      num_cores = 1, sig_cutoffs = c(0.05))
        
        parameter_rep_summary[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- c(rep_i, MAG_num, pseudocount_set, abun_increase_set)
        
        clade.based_summary_df <- rbind(clade.based_summary_df, parameter_rep_summary)
        
      }
    }
  }
}

saveRDS(object = clade.based_summary_df, file = "summaries/MAG.based_clade.based_param.altered_summary.rds")





# taxa.based
rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(pheatmap)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")


# Create dummy df for results based on single replicate.
test_taxa.based_POMS_rds <- list(output = readRDS("POMS_out/taxa.based/POMS_func_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_taxa.based_POMS_rds_list <- list(test_taxa.based_POMS_rds)

test_taxa.based_regress_rds <- list(output = readRDS("regress_out/taxa.based/regress_taxa.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds"))
test_taxa.based_regress_rds_list <- list(test_taxa.based_regress_rds$output)

test_taxa.based_wilcoxon_rds <- readRDS("wilcoxon_out/taxa.based/wilcoxon_taxa.based_out_rep10_MAGs1000_pseudo0_increase1.05.rds")
test_taxa.based_wilcoxon_rds_list <- list(test_taxa.based_wilcoxon_rds)

ex_func <- readRDS("prepped_func_tables/subset_1000MAGs_func_rep10.rds")

taxa.based_summary_df <- simulation_summaries(POMS_sims = test_taxa.based_POMS_rds_list,
                                               regress_sims = test_taxa.based_regress_rds_list,
                                               alt_tool_sims = test_taxa.based_wilcoxon_rds_list,
                                               focal_func_present = FALSE,
                                               func_table = ex_func,
                                               num_cores = 20,
                                               sig_cutoffs = c(0.05))

taxa.based_summary_df[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- NA
taxa.based_summary_df <- taxa.based_summary_df[-1, ]

parameter_settings <- list()

MAG_nums <- c(1595, 1000, 500, 250, 100)

pseudocount_settings <- c(0, 0.3, 0.7, 1)

abun_increase_settings <- c(1.5, 1.25, 1.05)

num_reps <- 10

for (rep_i in 1:num_reps) {
  
  for (MAG_num in MAG_nums) {
    
    POMS_output <- list()
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    for (pseudocount_set in pseudocount_settings) {
      
      for (abun_increase_set in abun_increase_settings) {
        
        
        POMS_rds_input <- list(output = readRDS(paste("POMS_out/taxa.based/POMS_func_out_",
                                                      "rep", as.character(rep_i),
                                                      "_MAGs", as.character(MAG_num),
                                                      "_pseudo", as.character(pseudocount_set),
                                                      "_increase", as.character(abun_increase_set),
                                                      ".rds", sep = "")))
        
        regress_rds_input <- list(output = readRDS(paste("regress_out/taxa.based/regress_taxa.based_out_",
                                                         "rep", as.character(rep_i),
                                                         "_MAGs", as.character(MAG_num),
                                                         "_pseudo", as.character(pseudocount_set),
                                                         "_increase", as.character(abun_increase_set),
                                                         ".rds", sep = "")))
        
        wilcoxon_rds_input <- list(readRDS(paste("wilcoxon_out/taxa.based/wilcoxon_taxa.based_out_",
                                                 "rep", as.character(rep_i),
                                                 "_MAGs", as.character(MAG_num),
                                                 "_pseudo", as.character(pseudocount_set),
                                                 "_increase", as.character(abun_increase_set),
                                                 ".rds", sep = "")))
        
        POMS_output <- list(POMS_rds_input)
        
        regress_output <- list(regress_rds_input$output)
        
        parameter_rep_summary <- simulation_summaries(POMS_sims = POMS_output,
                                                      regress_sims = regress_output,
                                                      alt_tool_sims = wilcoxon_rds_input,
                                                      focal_func_present = FALSE,
                                                      func_table = MAG_rep_func,
                                                      num_cores = 1, sig_cutoffs = c(0.05))
        
        parameter_rep_summary[, c("rep", "MAGs", "pseudocount", "abun_increase")] <- c(rep_i, MAG_num, pseudocount_set, abun_increase_set)
        
        taxa.based_summary_df <- rbind(taxa.based_summary_df, parameter_rep_summary)
        
      }
    }
  }
}

saveRDS(object = taxa.based_summary_df, file = "summaries/MAG.based_taxa.based_param.altered_summary.rds")
