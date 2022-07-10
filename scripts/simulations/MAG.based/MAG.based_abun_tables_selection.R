
rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")
source("../../../scripts/POMS_manuscript_functions.R")
source("../../../scripts/alt_tool_functions.R")

library(parallel)

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

# Sample 1000 funcs to be under selection (done and then commented out in case commands are re-run).
# set.seed(141515)
# random_funcs <- sample(colnames(almeida_func), size = 1000)
# write.table(x = random_funcs, file = "1000_random_KOs.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

random_funcs <- read.table("1000_random_KOs.txt", stringsAsFactors = FALSE)$V1

# Save prepped infiles.
saveRDS(object = almeida_tree, file = "MAG.based_prepped_tree.rds")
saveRDS(object = almeida_func, file = "MAG.based_prepped_func.rds")
saveRDS(object = almeida_abun, file = "MAG.based_prepped_abun.rds")


# Get prepped tables for two cases: (1) selection on taxa encoding a function
# and (2) selection on random taxa, but same number as in case #2 per each replicate.
func_sim_info <- mclapply(1:1000,
                                 
                          function(rep_i) {
                            
                            group1_samples <- random_groups$group1[[rep_i]]
                            group2_samples <- random_groups$group2[[rep_i]]
                            
                            func <- random_funcs[rep_i]
                            
                            orig_contributors <- rownames(almeida_func)[which(almeida_func[, func] > 0)]
                            
                            tmp_abun <- almeida_abun
                            
                            # Sanity check that no samples are missing.
                            missing_samples <- which(!c(group1_samples, group2_samples) %in% colnames(tmp_abun))
                            if (length(missing_samples) > 0) { stop("Error - there are samples missing.") }

                            tmp_abun[orig_contributors, group1_samples] <- (tmp_abun[orig_contributors, group1_samples] + 1) * 1.5
                            
                            rep_output <- list(taxa_perturb_abun = tmp_abun,
                                               orig_contributors = orig_contributors,
                                               func = func,
                                               group1 = group1_samples,
                                               group2 = group2_samples,
                                               selected_group = "group1")
                            
                            outfile <- paste("MAG.based_prepped_func_sim_info_sel1.5/", "func_sim_info_rep", as.character(rep_i), ".rds", sep = "")
                            
                            saveRDS(file = outfile, object = rep_output)
                            
                            return("Success")

                          }, mc.cores = 30)


taxa_sim_info <- mclapply(1:1000,
                              
                              function(rep_i) {

                                group1_samples <- random_groups$group1[[rep_i]]
                                group2_samples <- random_groups$group2[[rep_i]]
                                
                                func <- random_funcs[rep_i]
                                
                                orig_contributors <- rownames(almeida_func)[which(almeida_func[, func] > 0)]
        
                                ran_contributors <- sample(rownames(almeida_func), length(orig_contributors))
                                
                                tmp_abun <- almeida_abun
                                
                                # Sanity check that no samples are missing.
                                missing_samples <- which(!c(group1_samples, group2_samples) %in% colnames(tmp_abun))
                                if (length(missing_samples) > 0) { stop("Error - there are samples missing.") }
                                
                                tmp_abun[ran_contributors, group1_samples] <- (tmp_abun[ran_contributors, group1_samples] + 1) * 1.5
                                
                                rep_output <- list(taxa_perturb_abun = tmp_abun,
                                                   orig_contributors = orig_contributors,
                                                   ran_contributors = ran_contributors,
                                                   func = func,
                                                   group1 = group1_samples,
                                                   group2 = group2_samples,
                                                   selected_group = "group1")
                                
                                outfile <- paste("MAG.based_prepped_taxa_sim_info_sel1.5/", "taxa_sim_info_rep", as.character(rep_i), ".rds", sep = "")
                                
                                saveRDS(file = outfile, object = rep_output)
                                
                                return("Success")

                              }, mc.cores = 30)
