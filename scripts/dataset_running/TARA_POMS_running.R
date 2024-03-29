# Run POMS on TARA oceans dataset.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

node_balances_spearman_cor <- function(node_balances, sample_info, sample_var) {

  if (length(which(is.na(sample_info[, sample_var]))) > 0) {
    sample_info <- sample_info[-which(is.na(sample_info[, sample_var])), , drop = FALSE]
  }
  
  node_spearman_cor <- data.frame(matrix(NA, nrow = length(node_balances$balances), ncol = 2))
  colnames(node_spearman_cor) <- c("rho", "p")
  rownames(node_spearman_cor) <- names(node_balances$balances)
  
  for (node in names(node_balances$balances)) {
    
    cor_out <- cor.test(node_balances$balances[[node]][rownames(sample_info)],
                        sample_info[ , sample_var],
                        method = "spearman", exact = FALSE)
    
    node_spearman_cor[node, ] <- c(cor_out$estimate, cor_out$p.value) 
  }
  
  return(node_spearman_cor)
}

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


descrip_tables <- list()
descrip_tables[["ko"]] <- "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz"
descrip_tables[["pathway"]] <- "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz"
descrip_tables[["module"]] <- "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz"


TARA_ko <- POMS::filter_rare_table_cols(in_tab = TARA_ko, min_nonzero_count = 5, min_nonzero_prop = 0.001)
TARA_pathways <- POMS::filter_rare_table_cols(in_tab = TARA_pathways, min_nonzero_count = 5, min_nonzero_prop = 0.001)
TARA_modules <- POMS::filter_rare_table_cols(in_tab = TARA_modules, min_nonzero_count = 5, min_nonzero_prop = 0.001)

TARA_func <- list()
TARA_func[["ko"]] <- TARA_ko
TARA_func[["pathway"]] <- TARA_pathways
TARA_func[["module"]] <- TARA_modules


TARA_abun <- read.table(file = "NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/modified/TARA_abundance_min_mean_coverage1.tsv",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

TARA_sample_info <- read.table("Table_S1_sample_info.txt",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)


# Identify significant nodes outside of main POMS pipeline, based on Spearman correlations with sample metadata.
TARA_tree <- POMS::prep_tree(phy = TARA_tree, tips2keep = TARA_tree$tip.label)

TARA_node_balances <- POMS::compute_node_balances(tree = TARA_tree, 
                                                  abun_table = TARA_abun, 
                                                  ncores = 40,
                                                  pseudocount = 1,
                                                  min_num_tips = 10)

TARA_POMS_out <- list()

var2compare <- c("Chlorophyll.Sensor.s", "Mean_Temperature", 
                 "Mean_Salinity", "Mean_Oxygen", "Mean_Nitrates",
                 "NO2", "PO4", "NO2NO3", "SI")

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
  TARA_abun_subset <- TARA_abun_subset[which(rowSums(TARA_abun_subset) > 0), which(colSums(TARA_abun_subset) > 0)]
  
  for (func_level in names(TARA_func)) {
    
    TARA_func_table <- TARA_func[[func_level]][rownames(TARA_abun_subset), ]
    TARA_func_table <- TARA_func_table[, which(colSums(TARA_func_table) > 0)]
    
    TARA_POMS_out[[sample_var]][[func_level]] <- POMS_pipeline(abun = TARA_abun_subset, func = TARA_func_table, tree = TARA_tree, ncores = 10, pseudocount = 1,
                                                               manual_BSNs = var_sig_nodes, manual_balances = TARA_node_balances$balances, manual_BSN_dir = var_nodes_dir[var_sig_nodes],
                                                               min_func_instances = 5, min_func_prop = 0.001, func_descrip_infile = descrip_tables[[func_level]])
  }

}


# Quick dig into the results:
for (sample_var in names(TARA_POMS_out)) {
  for (func_type in names(TARA_POMS_out[[sample_var]])) {
   
    sig_subset <- TARA_POMS_out[[sample_var]][[func_type]]$results[which(TARA_POMS_out[[sample_var]][[func_type]]$results$multinomial_corr < 0.25), ]
    
    if (nrow(sig_subset) > 0) {
      print(func_type)
      print(sample_var)
      print(sig_subset)
      print("")
      print("")
      
    }
    
    
  }
  

}

saveRDS(object = TARA_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/TARA_POMS_out.rds")
