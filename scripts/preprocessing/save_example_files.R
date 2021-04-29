### Save example input files for running POMS.

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

test_tree <- readRDS("MAG.based_prepped_tree.rds")

tip_subset <- test_tree$tip.label[661:720]

test_tree <- drop.tip(phy = test_tree,
                      tip = test_tree$tip.label[which(!test_tree$tip.label %in% tip_subset)],
                      trim.internal = TRUE)

test_func <- readRDS("MAG.based_prepped_func.rds")

test_func <- test_func[, c("K07106", "K02036")]

func_sim_info <- readRDS("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep76.rds")

test_taxa <- func_sim_info$taxa_perturb_abun[test_tree$tip.label, ]

test_taxa <- test_taxa[, -which(colSums(test_taxa) == 0)]

test_taxa[test_taxa == 1.5] <- rnorm(n = length(which(test_taxa == 1.5)), mean = 1.5, sd = 0.5)

group1_samples <- func_sim_info$group1[which(func_sim_info$group1 %in% colnames(test_taxa))]
group2_samples <- func_sim_info$group2[which(func_sim_info$group2 %in% colnames(test_taxa))]

test_taxa <- test_taxa[, c(group1_samples, group2_samples)]

write.tree(phy = test_tree, file = "/home/gavin/github_repos/POMS/example_files/ex_tree.newick")

# I manually gzipped these two files:
write.table(x = group1_samples, file = "/home/gavin/github_repos/POMS/example_files/ex_group1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(x = group2_samples, file = "/home/gavin/github_repos/POMS/example_files/ex_group2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Note that I gzipped and added "genome" to be the first column manually in the below two files
write.table(x = test_func, file = "/home/gavin/github_repos/POMS/example_files/ex_func.tsv", sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
write.table(x = test_taxa, file = "/home/gavin/github_repos/POMS/example_files/ex_taxa_abun.tsv", sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
