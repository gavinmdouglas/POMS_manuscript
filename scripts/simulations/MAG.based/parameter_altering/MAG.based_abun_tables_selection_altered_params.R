# Create simulated tables based on altering the number of MAGs and the "selection pressure".

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

# Read in pre-determined random sample groupings.
random_groups <- readRDS("almeida_healthy_random_groups.rds")

# Read in Almeida 2019 MAG fiels
almeida_func <- read.table(file = "../../key_inputs/Almeida2019_dataset/functional_analyses/kegg_summary.csv.gz",
                           sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
almeida_func <- data.frame(t(almeida_func), check.names = FALSE)

almeida_tree <- ape::read.tree(file = "../../key_inputs/Almeida2019_dataset/phylogenies/raxml_hgr-umgs_phylogeny.nwk")
almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- POMS::subset_by_col_and_filt(in_tab = almeida_abun,
                                             col2keep = c(random_groups$group1[[1]], random_groups$group2[[1]]))

almeida_func <- almeida_func[rownames(almeida_abun), ]
almeida_func <- almeida_func[, which(colSums(almeida_func) != 0)]

# Remove missing tips from tree.
tips2drop <- almeida_tree$tip.label[which(!almeida_tree$tip.label %in% rownames(almeida_abun))]
almeida_tree <- ape::multi2di(ape::drop.tip(phy = almeida_tree, tip = tips2drop, trim.internal = TRUE))

# Filter rare functions from table.
almeida_func <- filter_rare_table_cols(in_tab = almeida_func, min_nonzero_count = 5, min_nonzero_prop = 0.001, verbose = TRUE)

MAG_nums <- c(1595, 1000, 500, 250, 100)

pseudocount_settings <- c(0, 0.3, 0.7, 1)

abun_increase_settings <- c(1.5, 1.25, 1.05)

num_reps <- 10

# Choose 10 random functions that will be used for each replicate
# (e.g., replicate 1 will always be the same KO, to keep things more comparable) 
# Restrict to functions that are found in at least 100 genomes.

# Commenting out since only needed to be run once.
# almeida_func_atleast100 <- almeida_func[, which(colSums(almeida_func > 0) >= 100)]
# rep_sampled_genes <- sample(colnames(almeida_func_atleast100), num_reps)

# boxplot(colSums(almeida_func_atleast100 > 0), colSums(almeida_func_atleast100[, rep_sampled_genes] > 0))

# Based on above boxplot - chose set of sampled genes with similar distribution of prevalence as the background of all genes.
# Played around with re-sampling until the distribution seemed representative.

# Then wrote to file so that this step could be skipped if this script needs to be re-run.
#write.table(x = rep_sampled_genes, file = "parameter_altered_files/sampled_10_functions_reps.txt",
#            quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE).

rep_sampled_genes <- read.table(file = "parameter_altered_files/sampled_10_functions_reps.txt", stringsAsFactors = FALSE)$V1

for (rep_i in 1:num_reps) {
  
   sampled_func <- rep_sampled_genes[rep_i]
   
   original_prev <- sum(almeida_func[, sampled_func] > 0) / nrow(almeida_func)
  
   original_encoders <- rownames(almeida_func)[which(almeida_func[, sampled_func] > 0)]
   
   original_nonencoder <- rownames(almeida_func)[which(almeida_func[, sampled_func] == 0)]
   
   group1_samples <- random_groups$group1[[rep_i]]
   group2_samples <- random_groups$group2[[rep_i]]
   
  for (MAG_num in MAG_nums) {
  
    if (MAG_num == 1595) {
      abun_MAG_subset <- almeida_abun
    } else {
      
      # Keep the prevalence of genomes with the function close to the original prevalence.
      # While also making sure that there are a min of 5 genomes that encode the function.
  
      num_required_encoding_genomes <- round(original_prev * MAG_num)
      
      if (num_required_encoding_genomes < 5) {
        num_required_encoding_genomes <- 5
      }
      
      sampled_encoders <- sample(x = original_encoders, size = num_required_encoding_genomes, replace = FALSE)
      sampled_nonencoders <- sample(x = original_nonencoder, size = MAG_num - num_required_encoding_genomes, replace = FALSE)
      
      abun_MAG_subset <- almeida_abun[c(sampled_encoders, sampled_nonencoders), ]
      abun_MAG_subset <- abun_MAG_subset[, which(colSums(abun_MAG_subset) > 0)]
      
    }
    
    group1_samples <- group1_samples[which(group1_samples %in% colnames(abun_MAG_subset))]
    group2_samples <- group2_samples[which(group2_samples %in% colnames(abun_MAG_subset))]
  
    MAG_subset_tree <- POMS::prep_tree(phy = almeida_tree, tips2keep = rownames(abun_MAG_subset))
    
    MAG_subset_func <- almeida_func[rownames(abun_MAG_subset), ]

    # Filter rare functions from table.
    MAG_subset_func <- filter_rare_table_cols(in_tab = MAG_subset_func, min_nonzero_count = 5, min_nonzero_prop = 0.001, verbose = TRUE)
    
    # Save intermediate files that will be used later.
    saveRDS(object = MAG_subset_tree,
            file = paste("parameter_altered_files/prepped_trees/subset_", as.character(MAG_num), "MAGs_tree_rep", as.character(rep_i), ".rds", sep = ""))
    
    saveRDS(object = MAG_subset_func,
            file = paste("parameter_altered_files/prepped_func_tables/subset_", as.character(MAG_num), "MAGs_func_rep", as.character(rep_i), ".rds", sep = ""))

    # Loop through different selection settings.
    for (pseudocount_set in pseudocount_settings) {
     
       for (abun_increase_set in abun_increase_settings) {
    
         contributors <- rownames(MAG_subset_func)[which(MAG_subset_func[, sampled_func] > 0)]
         
         rep_abun_MAG_subset <- abun_MAG_subset
         
         num_to_bump <- floor(length(contributors) * pseudocount_set)
         if (num_to_bump > 0) {
           for (g1_samp in group1_samples) {
             
             increased_contributors <- sample(contributors, size = num_to_bump, replace = FALSE)
             
             rep_abun_MAG_subset[increased_contributors, g1_samp] <- rep_abun_MAG_subset[increased_contributors, g1_samp] + 1
             
           }
         }
         
         rep_abun_MAG_subset[contributors, group1_samples] <- rep_abun_MAG_subset[contributors, group1_samples] * abun_increase_set

         rep_output <- list(taxa_perturb_abun = rep_abun_MAG_subset,
                            orig_contributors = contributors,
                            func = sampled_func,
                            group1 = group1_samples,
                            group2 = group2_samples,
                            pseudocount_set = pseudocount_set,
                            abun_increase = abun_increase_set,
                            selected_group = "group1")
         
         outfile <- paste("parameter_altered_files/sim_info/func.based/func.based_sim_info_",
                          "rep", as.character(rep_i),
                          "_MAGs", as.character(MAG_num),
                          "_pseudo", as.character(pseudocount_set),
                          "_increase", as.character(abun_increase_set),
                          ".rds", sep = "")
         
         saveRDS(file = outfile, object = rep_output)
              
       }
    }
  }
}
