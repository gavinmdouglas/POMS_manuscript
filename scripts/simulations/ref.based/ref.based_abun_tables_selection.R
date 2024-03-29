
rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(parallel)

# Read in pre-determined random sample groupings.
random_groups <- readRDS("ref.based_random_groups.rds")

BEZI_tables <- list()

BEZI_tables[["nu0.50"]] <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.50_sigma1.tsv.gz",
                                            header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

BEZI_tables[["nu0.65"]] <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.65_sigma1.tsv.gz",
                                            header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

BEZI_tables[["nu0.80"]] <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.80_sigma1.tsv.gz",
                                             header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

BEZI_tables[["nu0.95"]] <- read.table(file = "relabun_tables/BEZI_table_mu0.1_nu0.95_sigma1.tsv.gz",
                                             header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

in_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                      sep = "\t", header = TRUE, row.names = 1)

in_tree <- read.tree(file = "GToTree_output_subset.tre")



cutoff_func <- list()
cutoff_tree <- list()

for (cutoff in names(BEZI_tables)) {
  
  cutoff_func[[cutoff]] <- in_func[rownames(BEZI_tables[[cutoff]]), ]
  
  if (length(which(colSums(cutoff_func[[cutoff]]) == 0)) > 0) {
    cutoff_func[[cutoff]] <- cutoff_func[[cutoff]][, -which(colSums(cutoff_func[[cutoff]]) == 0)]
  }

  tips2drop <- in_tree$tip.label[which(!in_tree$tip.label %in% rownames(BEZI_tables[[cutoff]]))]

  if (length(tips2drop) > 0) {
    cutoff_tree[[cutoff]] <- multi2di(drop.tip(phy = in_tree, tip = tips2drop, trim.internal = TRUE))
  } else {
    cutoff_tree[[cutoff]] <- multi2di(in_tree)
  }
  
}

saveRDS(object = cutoff_tree, "ref.based_prepped_trees.rds")
saveRDS(object = cutoff_func, "ref.based_prepped_func_tables.rds")


# strictest_cutoff <- "nu0.95"
# set.seed(141515)
# in_func_no_rare <- filter_rare_table_cols(in_tab = cutoff_func[[strictest_cutoff]], min_nonzero_count = 5, min_nonzero_prop = 0.001, verbose = TRUE)
# random_funcs <- sample(colnames(in_func_no_rare), size = 1000)
# write.table(x = random_funcs, file = "1000_random_KOs.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

random_funcs <- read.table(file = "1000_random_KOs.txt", stringsAsFactors = FALSE)$V1

taxa_sim_info <- list()
func_sim_info <- list()

for (cutoff in names(BEZI_tables)) {
  
  taxa_sim_info[[cutoff]] <- mclapply(1:1000,
  
                                                            function(rep_i) {
  
                                                              group1_samples <- random_groups$group1[[rep_i]]
                                                              group2_samples <- random_groups$group2[[rep_i]]
  
                                                              func <- random_funcs[rep_i]
  
                                                              orig_contributors <- rownames(cutoff_func[[cutoff]])[which(cutoff_func[[cutoff]][, func] > 0)]
  
                                                              ran_contributors <- sample(rownames(cutoff_func[[cutoff]]), length(orig_contributors))
  
                                                              tmp_abun <- BEZI_tables[[cutoff]]
  
                                                              missing_group1 <- which(!group1_samples %in% colnames(tmp_abun))
                                                              missing_group2 <- which(!group2_samples %in% colnames(tmp_abun))
  
                                                              if (length(missing_group1) > 0) { group1_samples <- group1_samples[-missing_group1] }
                                                              if (length(missing_group2) > 0) { group2_samples <- group2_samples[-missing_group2] }
  
                                                              tmp_abun[ran_contributors, group1_samples] <- (tmp_abun[ran_contributors, group1_samples] + 1) * 1.5
                                                              
                                                              sim_info <- list(taxa_perturb_abun = tmp_abun,
                                                                               orig_contributors = orig_contributors,
                                                                               ran_contributors = ran_contributors,
                                                                               func = func,
                                                                               group1 = group1_samples,
                                                                               group2 = group2_samples,
                                                                               selected_group = "group1")
                                                              
                                                              sim_info_outfile <- paste("taxa_sim_info_sel1.5/taxa_sim_info_cutoff_", cutoff, "_rep", rep_i, ".rds", sep = "")
                                                              
                                                              saveRDS(object = sim_info, file = sim_info_outfile)
                                                              
                                                              return("Success")
                                                            }, mc.cores = 1)


    func_sim_info[[cutoff]] <- mclapply(1:1000,
                                 
                                                   function(rep_i) {
                                                     
                                                     group1_samples <- random_groups$group1[[rep_i]]
                                                     group2_samples <- random_groups$group2[[rep_i]]
                                                     
                                                     func <- random_funcs[rep_i]
                                                     
                                                     orig_contributors <- rownames(cutoff_func[[cutoff]])[which(cutoff_func[[cutoff]][, func] > 0)]
                                                     
                                                     tmp_abun <- BEZI_tables[[cutoff]]
                                                     
                                                     missing_group1 <- which(!group1_samples %in% colnames(tmp_abun))
                                                     missing_group2 <- which(!group2_samples %in% colnames(tmp_abun))
                                                     
                                                     if (length(missing_group1) > 0) { group1_samples <- group1_samples[-missing_group1] }
                                                     if (length(missing_group2) > 0) { group2_samples <- group2_samples[-missing_group2] }
                                                     
                                                     tmp_abun[orig_contributors, group1_samples] <- (tmp_abun[orig_contributors, group1_samples] + 1) * 1.5
                                                     
                                                     sim_info <- list(taxa_perturb_abun = tmp_abun,
                                                                      orig_contributors = orig_contributors,
                                                                      func = func,
                                                                      group1 = group1_samples,
                                                                      group2 = group2_samples,
                                                                      selected_group = "group1")
                                                     
                                                     sim_info_outfile <- paste("func_sim_info_sel1.5/func_sim_info_cutoff_", cutoff, "_rep", rep_i, ".rds", sep = "")
                                                     
                                                     saveRDS(object = sim_info, file = sim_info_outfile)
                                                     
                                                     return("Success")
                                                   }, mc.cores = 1)

}
