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

almeida_abun <- subset_abun_table(in_abun = almeida_abun,
                                  col2keep = c(random_groups$group1[[1]], random_groups$group2[[1]]))

almeida_func <- almeida_func[rownames(almeida_abun), ]
almeida_func <- almeida_func[, -which(colSums(almeida_func) == 0)]

# Remove missing tips from tree.
tips2drop <- almeida_tree$tip.label[which(!almeida_tree$tip.label %in% rownames(almeida_abun))]
almeida_tree <- ape::multi2di(ape::drop.tip(phy = almeida_tree, tip = tips2drop, trim.internal = TRUE))

# Filter rare functions from table.
almeida_func <- filter_rare_table_cols(in_tab = almeida_func, min_nonzero_count = 5, min_nonzero_prop = 0.001, verbose = TRUE)

MAG_nums <- c(1595, 1250, 1000, 750, 500, 250, 100, 50)

pseudocount_settings <- c(0, 1)

abun_increase_settings <- c(1.5, 1.3, 1.1, 1.05)

sample_count <- data.frame(matrix(NA, nrow = 25, ncol = 8))

focal_genes <- list()
for (MAG_num in MAG_nums) { focal_genes[[as.character(MAG_num)]] <- as.character() }

for (rep_i in 1:25) {
  
  for (MAG_num in MAG_nums) {
  
    almeida_abun_MAG_subset <- almeida_abun[sample(x = rownames(almeida_abun), size = MAG_num, replace = FALSE), ]
    
    if (length(which(colSums(almeida_abun_MAG_subset) == 0)) > 0) {
      almeida_abun_MAG_subset <- almeida_abun_MAG_subset[, -which(colSums(almeida_abun_MAG_subset) == 0)]
    }
    
    group1_samples <- random_groups$group1[[rep_i]]
    group2_samples <- random_groups$group2[[rep_i]]
    
    group1_samples <- group1_samples[which(group1_samples %in% colnames(almeida_abun_MAG_subset))]
    group2_samples <- group2_samples[which(group2_samples %in% colnames(almeida_abun_MAG_subset))]
    
    sample_count[rep_i, which(MAG_nums == MAG_num)] <- ncol(almeida_abun_MAG_subset)
   
    MAG_subset_tips2drop <- almeida_tree$tip.label[which(!almeida_tree$tip.label %in% rownames(almeida_abun_MAG_subset))]
    MAG_subset_tree <- ape::multi2di(ape::drop.tip(phy = almeida_tree, tip = MAG_subset_tips2drop, trim.internal = TRUE))
    
    MAG_subset_func <- almeida_func[rownames(almeida_abun_MAG_subset), ]
    if (length(which(colSums(MAG_subset_func) == 0)) > 0) {
      MAG_subset_func <- MAG_subset_func[, -which(colSums(MAG_subset_func) == 0)]
    }
    
    # Filter rare functions from table.
    MAG_subset_func <- filter_rare_table_cols(in_tab = MAG_subset_func, min_nonzero_count = 5, min_nonzero_prop = 0.001, verbose = TRUE)
    
    # Save intermediate files that will be used later.
    saveRDS(object = MAG_subset_tree,
            file = paste("parameter_altered_files/prepped_trees/subset_", as.character(MAG_num), "MAGs_tree_rep", as.character(rep_i), ".rds", sep = ""))
    
    saveRDS(object = MAG_subset_func,
            file = paste("parameter_altered_files/prepped_func_tables/subset_", as.character(MAG_num), "MAGs_func_rep", as.character(rep_i), ".rds", sep = ""))
    
    focal_gene <- sample(x = colnames(MAG_subset_func), size = 1)
    
    # Make sure this focal func isn't used for a different replicate.
    while (focal_gene %in% focal_genes[[as.character(MAG_num)]]) {
      focal_gene <- sample(x = colnames(MAG_subset_func), size = 1)
    }
    
    focal_genes[[as.character(MAG_num)]] <- c(focal_genes[[as.character(MAG_num)]], focal_gene)

    # Loop through different selection settings.
    for (pseudocount_set in pseudocount_settings) {
     
       for (abun_increase_set in abun_increase_settings) {
    
         orig_contributors <- rownames(MAG_subset_func)[which(MAG_subset_func[, focal_gene] > 0)]
         
         rep_abun_MAG_subset <- almeida_abun_MAG_subset
         
         rep_abun_MAG_subset[orig_contributors, group1_samples] <- (rep_abun_MAG_subset[orig_contributors, group1_samples] + pseudocount_set) * abun_increase_set
         
         rep_output <- list(taxa_perturb_abun = rep_abun_MAG_subset,
                            orig_contributors = orig_contributors,
                            func = focal_gene,
                            group1 = group1_samples,
                            group2 = group2_samples,
                            pseudocount_set = pseudocount_set,
                            abun_increase = abun_increase_set,
                            selected_group = "group1")
         
         outfile <- paste("parameter_altered_files/sim_info/func_sim_info_",
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
