rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

MAG_nums <- c(1595, 1000, 500, 250, 100)

mean_num <- c()
for (MAG_num in MAG_nums) {
    
  num_func <- c()
  
  for (rep_i in 1:10) {
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
 
    num_func <- c(num_func, ncol(MAG_rep_func))
       
  }
  
  mean_num <- c(mean_num, mean(num_func))
}
    
    