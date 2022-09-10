# Run spearman correlations against function abundances as a comparison on TARA oceans dataset.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")
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

TARA_abun <- read.table(file = "NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/modified/TARA_abundance_min_mean_coverage1.tsv",
                        header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

TARA_sample_info <- read.table("Table_S1_sample_info.txt",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

# Generate function abundance tables.
TARA_func_abun <- list()

# Filter out rare functions from each table.
TARA_ko <- POMS::filter_rare_table_cols(in_tab = TARA_ko, min_nonzero_count = 5, min_nonzero_prop = 0.001)
TARA_pathways <- POMS::filter_rare_table_cols(in_tab = TARA_pathways, min_nonzero_count = 5, min_nonzero_prop = 0.001)
TARA_modules <- POMS::filter_rare_table_cols(in_tab = TARA_modules, min_nonzero_count = 5, min_nonzero_prop = 0.001)

TARA_func_abun[["ko"]] <- calc_func_abun_crossprod(in_abun = TARA_abun, in_func = TARA_ko)
TARA_func_abun[["pathway"]] <- calc_func_abun_crossprod(in_abun = TARA_abun, in_func = TARA_pathways)
TARA_func_abun[["module"]] <- calc_func_abun_crossprod(in_abun = TARA_abun, in_func = TARA_modules)

TARA_func_abun[["ko"]] <- TARA_func_abun[["ko"]][which(rowSums(TARA_func_abun[["ko"]]) > 0), ]
TARA_func_abun[["pathway"]] <- TARA_func_abun[["pathway"]][which(rowSums(TARA_func_abun[["pathway"]]) > 0), ]
TARA_func_abun[["module"]] <- TARA_func_abun[["module"]][which(rowSums(TARA_func_abun[["module"]]) > 0), ]

TARA_spearman_out <- list()

var2compare <- c("Chlorophyll.Sensor.s", "Mean_Temperature", "Mean_Salinity", "Mean_Oxygen", "Mean_Nitrates", "NO2", "PO4", "NO2NO3", "SI")

for (sample_var in var2compare) {

  print(sample_var)
  
  if (length(which(is.na(TARA_sample_info[, sample_var]))) > 0) {
    TARA_sample_info_subset <- TARA_sample_info[-which(is.na(TARA_sample_info[, sample_var])), ]
  }
  
  TARA_spearman_out[[sample_var]] <- list()

  TARA_abun_subset <- TARA_abun[, rownames(TARA_sample_info_subset)]
  TARA_abun_subset <- TARA_abun_subset[which(rowSums(TARA_abun_subset) > 0), which(colSums(TARA_abun_subset) > 0)]

  for (func in c("ko", "pathway", "module")) {
  
    TARA_func_abun_subset <- TARA_func_abun[[func]][, colnames(TARA_abun_subset)]
    
    TARA_func_abun_subset <- TARA_func_abun_subset[which(rowSums(TARA_func_abun_subset) > 0), ]
    
    sample_var_vec <- as.numeric(TARA_sample_info_subset[, sample_var])
    names(sample_var_vec) <- rownames(TARA_sample_info_subset)
    
    TARA_spearman_out[[sample_var]][[func]] <- spearman_cor_df_vs_vector(in_df = TARA_func_abun_subset, in_vec = sample_var_vec, p_corr = 'BH')
    
  }
  
}

### Quick look at results:
# for (sample_var in names(TARA_spearman_out)) {
#   print(sample_var)
#   
#   for (func_type in names(TARA_spearman_out[[sample_var]])) {
#     
#     print(func_type)
#     print(TARA_spearman_out[[sample_var]][[func_type]][which(TARA_spearman_out[[sample_var]][[func_type]]$p_corr < 0.05), ])
#     
#   }
#   
#   print("")
#   print("")
#   
# }

saveRDS(object = TARA_spearman_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/TARA_spearman_out.rds")


# Info for main-text:
rm(list = ls(all.names = TRUE))
setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

TARA_spearman_out <- readRDS("/home/gavin/github_repos/POMS_manuscript/data/results/TARA_spearman_out.rds")

length(which(TARA_spearman_out$Mean_Salinity$ko$p_corr < 0.25))
length(which(TARA_spearman_out$PO4$ko$p_corr < 0.25))

combined_TARA_sig_df <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/display_items/Maintext_TARA_sig_func_RAW.tsv",
                                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)

POMS_salinity_pathway <- combined_TARA_sig_df[which(combined_TARA_sig_df$func_type == "pathway" & combined_TARA_sig_df$variable == "Mean_Salinity"), "func"]
POMS_salinity_module <- combined_TARA_sig_df[which(combined_TARA_sig_df$func_type == "module" & combined_TARA_sig_df$variable == "Mean_Salinity"), "func"]

POMS_phosphate_pathway <- combined_TARA_sig_df[which(combined_TARA_sig_df$func_type == "pathway" & combined_TARA_sig_df$variable == "PO4"), "func"]
POMS_phosphate_ko <- combined_TARA_sig_df[which(combined_TARA_sig_df$func_type == "ko" & combined_TARA_sig_df$variable == "PO4"), "func"]

length(which(TARA_spearman_out$Mean_Salinity$pathway[POMS_salinity_pathway, "p_corr"] < 0.25))
length(which(TARA_spearman_out$Mean_Salinity$module[POMS_salinity_module, "p_corr"] < 0.25))

length(which(TARA_spearman_out$PO4$pathway[POMS_phosphate_pathway, "p_corr"] < 0.25))

length(which(TARA_spearman_out$PO4$ko[POMS_phosphate_ko, "p_corr"] < 0.25))

