# Run POMS on TARA oceans dataset.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

library(ape)
library(parallel)


# Read in input files.
TARA_ko <- read.table(file = "Table_S11_KO_abun.txt",
                         sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)

TARA_ko <- data.frame(t(TARA_ko), check.names = FALSE)

TARA_pathways <- read.table(file = "kegg_pathways/path_abun_unstrat.tsv.gz",
                               sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
TARA_pathways <- data.frame(t(TARA_pathways), check.names = FALSE)

TARA_modules <- read.table(file = "kegg_modules/path_abun_unstrat.tsv.gz",
                              sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
TARA_modules <- data.frame(t(TARA_modules), check.names = FALSE)

TARA_tree <- read.tree(file = "GToTree_output/GToTree_output_modified.tre")


TARA_abun <- read.table(file = "NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/modified/TARA_abundance_min_mean_coverage1.tsv",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

TARA_sample_info <- read.table("Table_S1_sample_info.txt",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)


# Identify significant nodes outside of main POMS pipeline, based on Spearman correlations with sample metadata.

TARA_tree <- ape::makeNodeLabel(TARA_tree, method = "number", prefix = 'n')

TARA_node_balances <- compute_tree_node_balances(phylogeny = TARA_tree,
                                                 abun = TARA_abun,
                                                 ncores = 40,
                                                 pseudocount = 1,
                                                 min_num_tips = 10)

TARA_POMS_out <- list()

var2compare <- c("Chlorophyll.Sensor.s", "Mean_Temperature", "Mean_Salinity", "Mean_Oxygen", "Mean_Nitrates", "NO2", "PO4", "NO2NO3", "SI")

for (sample_var in var2compare) {

  TARA_POMS_out[[sample_var]] <- list()
  
  TARA_node_cor <- node_balances_spearman_cor(node_balances = TARA_node_balances,
                                               sample_info = TARA_sample_info,
                                               sample_var = sample_var)

  var_sig_nodes <- rownames(TARA_node_cor)[which(TARA_node_cor$p < 0.05)]

  var_nodes_dir <- rep(NA, nrow(TARA_node_cor))
  names(var_nodes_dir) <- rownames(TARA_node_cor)
  var_nodes_dir[which(TARA_node_cor[, "rho"] >= 0)] <- "group1"
  var_nodes_dir[which(TARA_node_cor[, "rho"] < 0)] <- "group2"

  TARA_abun_subset <- TARA_abun[, rownames(TARA_sample_info[which(!is.na(TARA_sample_info[, sample_var])), ])]
  if (min(rowSums(TARA_abun_subset)) == 0) { TARA_abun_subset <- TARA_abun_subset[-which(rowSums(TARA_abun_subset) == 0), ]}
  if (min(colSums(TARA_abun_subset)) == 0) { TARA_abun_subset <- TARA_abun_subset[, -which(colSums(TARA_abun_subset) == 0)]}
  
  TARA_pathways_subset <- TARA_pathways[rownames(TARA_abun_subset), ]
  
  if (length(which(colSums(TARA_pathways_subset > 0) > 0.85 * nrow(TARA_pathways_subset))) > 0) {
    TARA_pathways_subset <- TARA_pathways_subset[, -which(colSums(TARA_pathways_subset > 0) > 0.85 * nrow(TARA_pathways_subset))]
  }
  
  if (length(which(colSums(TARA_pathways_subset > 0) < 0.15 * nrow(TARA_pathways_subset))) > 0) {
    TARA_pathways_subset <- TARA_pathways_subset[, -which(colSums(TARA_pathways_subset > 0) < 0.15 * nrow(TARA_pathways_subset))]
  }
  
  TARA_POMS_out[[sample_var]][["pathway"]] <- POMS_pipeline_continuous(abun = TARA_abun,
                                                                       func = TARA_pathways_subset,
                                                                       phylogeny = TARA_tree,
                                                                       group1_samples = c("placeholder"),
                                                                       group2_samples = c("placeholder2"),
                                                                       ncores = 10,
                                                                       significant_nodes = var_sig_nodes,
                                                                       nodes_dir = var_nodes_dir,
                                                                       tested_balances = TARA_node_balances$balances,
                                                                       balance_p_cutoff = 0.05,
                                                                       balance_correction = "none",
                                                                       function_p_cutoff = 0.05,
                                                                       function_correction = "none",
                                                                       min_num_tips = 10,
                                                                       min_func_instances = 10,
                                                                       min_func_prop = 0.001,
                                                                       run_multinomial_test = TRUE,
                                                                       multinomial_correction = "BH",
                                                                       detailed_output = TRUE,
                                                                       verbose = TRUE,
                                                                       func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz")
                        
  TARA_modules_subset <- TARA_modules[rownames(TARA_abun_subset), ]
  
  if (length(which(colSums(TARA_modules_subset > 0) > 0.85 * nrow(TARA_modules_subset))) > 0) {
    TARA_modules_subset <- TARA_modules_subset[, -which(colSums(TARA_modules_subset > 0) > 0.85 * nrow(TARA_modules_subset))]
  }
  
  if (length(which(colSums(TARA_modules_subset > 0) < 0.15 * nrow(TARA_modules_subset))) > 0) {
    TARA_modules_subset <- TARA_modules_subset[, -which(colSums(TARA_modules_subset > 0) < 0.15 * nrow(TARA_modules_subset))]
  }
  
  TARA_POMS_out[[sample_var]][["module"]] <- POMS_pipeline_continuous(abun = TARA_abun,
                                                                      func = TARA_modules_subset,
                                                                      phylogeny = TARA_tree,
                                                                      group1_samples = c("placeholder"),
                                                                      group2_samples = c("placeholder2"),
                                                                      ncores = 10,
                                                                      significant_nodes = var_sig_nodes,
                                                                      nodes_dir = var_nodes_dir,
                                                                      tested_balances = TARA_node_balances$balances,
                                                                      balance_p_cutoff = 0.05,
                                                                      balance_correction = "none",
                                                                      function_p_cutoff = 0.05,
                                                                      function_correction = "none",
                                                                      min_num_tips = 10,
                                                                      min_func_instances = 10,
                                                                      min_func_prop = 0.001,
                                                                      run_multinomial_test = TRUE,
                                                                      multinomial_correction = "BH",
                                                                      detailed_output = TRUE,
                                                                      verbose = TRUE,
                                                                      func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz")
  
  
  TARA_ko_subset <- TARA_ko[rownames(TARA_abun_subset), ]
  
  if (length(which(colSums(TARA_ko_subset > 0) > 0.85 * nrow(TARA_ko_subset))) > 0) {
    TARA_ko_subset <- TARA_ko_subset[, -which(colSums(TARA_ko_subset > 0) > 0.85 * nrow(TARA_ko_subset))]
  }
  
  if (length(which(colSums(TARA_ko_subset > 0) < 0.15 * nrow(TARA_ko_subset))) > 0) {
    TARA_ko_subset <- TARA_ko_subset[, -which(colSums(TARA_ko_subset > 0) < 0.15 * nrow(TARA_ko_subset))]
  }
  
  TARA_POMS_out[[sample_var]][["ko"]] <- POMS_pipeline_continuous(abun = TARA_abun,
                                                                      func = TARA_ko_subset,
                                                                      phylogeny = TARA_tree,
                                                                      group1_samples = c("placeholder"),
                                                                      group2_samples = c("placeholder2"),
                                                                      ncores = 10,
                                                                      significant_nodes = var_sig_nodes,
                                                                      nodes_dir = var_nodes_dir,
                                                                      tested_balances = TARA_node_balances$balances,
                                                                      balance_p_cutoff = 0.05,
                                                                      balance_correction = "none",
                                                                      function_p_cutoff = 0.05,
                                                                      function_correction = "none",
                                                                      min_num_tips = 10,
                                                                      min_func_instances = 10,
                                                                      min_func_prop = 0.001,
                                                                      run_multinomial_test = TRUE,
                                                                      multinomial_correction = "BH",
                                                                      detailed_output = TRUE,
                                                                      verbose = TRUE,
                                                                      func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz")
  

}


for (sample_var in names(TARA_POMS_out)) {

  print(sample_var)

  for (func_type in names(TARA_POMS_out[[sample_var]])) {
   
    print(func_type)
    print(TARA_POMS_out[[sample_var]][[func_type]]$df[which(TARA_POMS_out[[sample_var]][[func_type]]$df$multinomial_corr < 0.25), ])
    
  }
  
  print("")
  print("")
   
}

saveRDS(object = TARA_POMS_out, file = "../../results/TARA_POMS_out.rds")
