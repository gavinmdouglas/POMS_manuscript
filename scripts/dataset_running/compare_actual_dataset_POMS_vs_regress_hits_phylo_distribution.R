rm(list = ls(all.names = TRUE))

# Compare the phylogenetic distribution of significant hits identified by POMS vs phylogenetic regression.

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(picante)
library(parallel)

# First compute Faith's phylogenetic diversity for all tested functions (based on the subset of MAGs used).
almeida_ko <- read.table(file = "functional_analyses/kegg_summary.csv.gz",
                         sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)

almeida_ko <- data.frame(t(almeida_ko), check.names = FALSE)

almeida_pathways <- read.table(file = "functional_analyses/modified/kegg_pathways/path_abun_unstrat.tsv.gz",
                               sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
almeida_pathways <- data.frame(t(almeida_pathways), check.names = FALSE)

almeida_modules <- read.table(file = "functional_analyses/modified/kegg_modules/path_abun_unstrat.tsv.gz",
                              sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
almeida_modules <- data.frame(t(almeida_modules), check.names = FALSE)


almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")
almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt.gz",
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")


studies <- c("ERP002061", "ERP012177", "ERP003612")

func_tables <- list()
func_tables[["kos"]] <- almeida_ko
func_tables[["pathways"]] <- almeida_pathways
func_tables[["modules"]] <- almeida_modules


pd_out <- list()
pd_out_edges1 <- list()

for (study in studies) {
  
  POMS_out[[study]] <- list()
  
  sample_info <- almeida_sample_info[which(almeida_sample_info$Study == study), ]
  
  abun_table <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = sample_info$Run)
  
  sample_info <- sample_info[which(sample_info$Run %in% colnames(abun_table)), ]
  group1_samples <- sample_info[which(sample_info$Health.state == "Diseased"), "Run"]
  group2_samples <- sample_info[which(sample_info$Health.state == "Healthy"), "Run"]
  
  for (f in names(func_tables)) {
    
    filt_func_table <-  POMS::filter_rare_table_cols(in_tab = func_tables[[f]][rownames(abun_table), ],
                                                     min_nonzero_count = 5,
                                                     min_nonzero_prop = 0.001)
    
    
    

almeida_regress_out <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_regress_specificity_output/combined_output.rds")

names(almeida_regress_out)[which(names(almeida_regress_out) == "ERP002061")] <- "ERP002061 | Obesity 1"
names(almeida_regress_out)[which(names(almeida_regress_out) == "ERP002061")] <- "ERP003612 | Obesity 2"
names(almeida_regress_out)[which(names(almeida_regress_out) == "ERP012177")] <- "ERP012177 | CRC"

