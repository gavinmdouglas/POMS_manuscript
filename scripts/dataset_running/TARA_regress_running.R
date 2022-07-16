# Run regress on TARA oceans dataset.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/regress_manuscript/data/key_inputs/TARA/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

abun_spearman_cor <- function(abun_table, sample_info, sample_var) {
  
  abun_table = TARA_abun_subset
  sample_info <- TARA_sample_info
  sample_var <- sample_var
  
  if (length(which(is.na(sample_info[, sample_var]))) > 0) {
    sample_info <- sample_info[-which(is.na(sample_info[, sample_var])), , drop = FALSE]
  }
  
  abun_spearman_cor <- data.frame(matrix(NA, nrow = nrow(abun_table), ncol = 2))
  colnames(abun_spearman_cor) <- c("rho", "p")
  rownames(abun_spearman_cor) <- rownames(abun_table)
  
  for (taxon in rownames(abun_table)) {
    
    cor_out <- cor.test(as.numeric(abun_table[taxon, rownames(sample_info)]),
                        sample_info[ , sample_var],
                        method = "spearman", exact = FALSE)
    
    abun_spearman_cor[taxon, ] <- c(cor_out$estimate, cor_out$p.value) 
  }
  
  return(abun_spearman_cor)
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
descrip_tables[["ko"]] <- read.table("/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz",
                                     header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")

descrip_tables[["pathway"]] <- read.table("/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz",
                                          header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")

descrip_tables[["module"]] <- read.table("/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz",
                                         header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")

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

TARA_tree <- POMS::prep_tree(phy = TARA_tree, tips2keep = TARA_tree$tip.label).

TARA_regress_out <- list()

var2compare <- c("Chlorophyll.Sensor.s", "Mean_Temperature", 
                 "Mean_Salinity", "Mean_Oxygen", "Mean_Nitrates",
                 "NO2", "PO4", "NO2NO3", "SI")

for (sample_var in var2compare) {

  TARA_regress_out[[sample_var]] <- list()
  
  TARA_abun_subset <- TARA_abun[, rownames(TARA_sample_info[which(!is.na(TARA_sample_info[, sample_var])), ])]
  TARA_abun_subset <- TARA_abun_subset[which(rowSums(TARA_abun_subset) > 0), which(colSums(TARA_abun_subset) > 0)]
  
  correlated_taxa <- abun_spearman_cor(abun_table = TARA_abun_subset, sample_info = TARA_sample_info, sample_var = sample_var)
  
  for (func_level in names(TARA_func)) {
    
    TARA_func_table <- TARA_func[[func_level]][rownames(TARA_abun_subset), ]
    TARA_func_table <- TARA_func_table[, which(colSums(TARA_func_table) > 0)]

    TARA_func_table <- POMS::filter_rare_table_cols(in_tab = TARA_func_table, min_nonzero_count = 5, min_nonzero_prop = 0.001)
    
    sig_taxa <- rep(0, nrow(correlated_taxa))
    sig_taxa[which(correlated_taxa$p < 0.05)] <- 1
    
    TARA_tree_subset <- POMS::prep_tree(phy = TARA_tree, tips2keep = rownames(TARA_func_table))
    
    var_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                    func =  TARA_func_table,
                                                    in_tree = TARA_tree_subset,
                                                    ncores = 10,
                                                    model_type = "BM")

    var_regress_out$BH <- p.adjust(var_regress_out$p, "BH")
    
    var_regress_out$Description <- descrip_tables[[func_level]][rownames(var_regress_out), "V2"]
    
    TARA_regress_out[[sample_var]][[func_level]] <- var_regress_out
    
  }

}


# Quick dig into the results:
for (sample_var in names(TARA_regress_out)) {
  for (func_type in names(TARA_regress_out[[sample_var]])) {
   
    sig_subset <- TARA_regress_out[[sample_var]][[func_type]][which(TARA_regress_out[[sample_var]][[func_type]]$BH < 0.25), ]
    
    if (nrow(sig_subset) > 0) {
      print(func_type)
      print(sample_var)
      print(sig_subset)
      print("")
      print("")
      
    }
    
    
  }
  

}

saveRDS(object = TARA_regress_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/TARA_regress_out.rds")
