# Run spearman correlations against function abundances as a comparison on TARA oceans dataset.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

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

TARA_abun <- read.table(file = "NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/modified/TARA_abundance_min_mean_coverage1.tsv",
                        header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

TARA_sample_info <- read.table("Table_S1_sample_info.txt",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

# Generate function abundance tables.
TARA_func_abun <- list()

TARA_func_abun[["ko"]] <- calc_func_abun(in_abun = TARA_abun, in_func = TARA_ko, ncores = 40)

TARA_func_abun[["pathway"]] <- calc_func_abun(in_abun = TARA_abun, in_func = TARA_pathways, ncores = 40)

TARA_func_abun[["module"]] <- calc_func_abun(in_abun = TARA_abun, in_func = TARA_modules, ncores = 40)


TARA_spearman_out <- list()

var2compare <- c("Chlorophyll.Sensor.s", "Mean_Temperature", "Mean_Salinity", "Mean_Oxygen", "Mean_Nitrates", "NO2", "PO4", "NO2NO3", "SI")

for (sample_var in var2compare) {

  print(sample_var)
  
  if (length(which(is.na(TARA_sample_info[, sample_var]))) > 0) {
    TARA_sample_info_subset <- TARA_sample_info[-which(is.na(TARA_sample_info[, sample_var])), ]
  }
  
  TARA_spearman_out[[sample_var]] <- list()

  TARA_abun_subset <- TARA_abun[, rownames(TARA_sample_info_subset)]
  if (min(rowSums(TARA_abun_subset)) == 0) { TARA_abun_subset <- TARA_abun_subset[-which(rowSums(TARA_abun_subset) == 0), ]}
  if (min(colSums(TARA_abun_subset)) == 0) { TARA_abun_subset <- TARA_abun_subset[, -which(colSums(TARA_abun_subset) == 0)]}

  for (func in c("ko", "pathway", "module")) {
  
    TARA_func_abun_subset <- TARA_func_abun[[func]][, colnames(TARA_abun_subset)]
  
    if (length(which(rowSums(TARA_func_abun_subset > 0) < 0.15 * ncol(TARA_func_abun_subset))) > 0) {
      TARA_func_abun_subset <- TARA_func_abun_subset[-which(rowSums(TARA_func_abun_subset > 0) < 0.15 * ncol(TARA_func_abun_subset)), ]
    }
    
    sample_var_vec <- as.numeric(TARA_sample_info_subset[, sample_var])
    names(sample_var_vec) <- rownames(TARA_sample_info_subset)
    
    TARA_spearman_out[[sample_var]][[func]] <- spearman_cor_df_vs_vector(in_df = TARA_func_abun_subset, in_vec = sample_var_vec, p_corr = 'BH')
    
  }
  
}


for (sample_var in names(TARA_spearman_out)) {
  print(sample_var)
  
  for (func_type in names(TARA_spearman_out[[sample_var]])) {
    
    print(func_type)
    print(TARA_spearman_out[[sample_var]][[func_type]][which(TARA_spearman_out[[sample_var]][[func_type]]$p_corr < 0.05), ])
    
  }
  
  print("")
  print("")
  
}

saveRDS(object = TARA_spearman_out, file = "../../results/TARA_spearman_out.rds")

