# Create simulated tables based on altering the number of MAGs and the "selection pressure".
# In this version of the simulation, add pseudocount of 1 to ALL taxa (in both sample groups) prior to multiplying by abundance factor.

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

# Read in pre-determined random sample groupings.
random_groups <- readRDS("almeida_healthy_random_groups.rds")

almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- subset_abun_table(in_abun = almeida_abun,
                                  col2keep = c(random_groups$group1[[1]], random_groups$group2[[1]]))

MAG_nums <- c(1595, 1250, 1000, 750, 500, 250, 100, 50)

abun_increase_settings <- c(1.5, 1.3, 1.1, 1.05)

for (rep_i in 1:25) {
  
  for (MAG_num in MAG_nums) {
    
    MAG_subset_func <- readRDS(paste("parameter_altered_files/prepped_func_tables/subset_",
                                     as.character(MAG_num),
                                     "MAGs_func_rep",
                                     as.character(rep_i),
                                     ".rds", sep = ""))
  
    representative_sim_info <- readRDS(paste("parameter_altered_files/sim_info/func_sim_info_rep",
                                             as.character(rep_i),
                                             "_MAGs",
                                             as.character(MAG_num),
                                             "_pseudo0_increase1.05.rds",
                                             sep = ""))
  
    focal_gene <- representative_sim_info$func
    group1_samples <- representative_sim_info$group1
    group2_samples <- representative_sim_info$group2
    
    for (abun_increase_set in abun_increase_settings) {
      
      orig_contributors <- rownames(MAG_subset_func)[which(MAG_subset_func[, focal_gene] > 0)]
      
      rep_abun_MAG_subset <- almeida_abun[rownames(MAG_subset_func), ]
      
      if (length(which(colSums(rep_abun_MAG_subset) == 0)) > 0) {
        rep_abun_MAG_subset <- rep_abun_MAG_subset[, -which(colSums(rep_abun_MAG_subset) == 0)]
      }
      
      # Key part of these simulations is to add a pseudocount to all taxa across all samples prior to multiplying by abundance factor.
      rep_abun_MAG_subset <- rep_abun_MAG_subset + 1
      
      rep_abun_MAG_subset[orig_contributors, group1_samples] <- rep_abun_MAG_subset[orig_contributors, group1_samples] * abun_increase_set
      
      rep_output <- list(taxa_perturb_abun = rep_abun_MAG_subset,
                         orig_contributors = orig_contributors,
                         func = focal_gene,
                         group1 = group1_samples,
                         group2 = group2_samples,
                         pseudocount_set = "prepseudo",
                         abun_increase = abun_increase_set,
                         selected_group = "group1")
      
      outfile <- paste("parameter_altered_files/sim_info_prepseudo/func_sim_info_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_prepseudo",
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = "")
      
      saveRDS(file = outfile, object = rep_output)
              
    }
  }
}
