rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

# Read in input files.
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

descrip_tables <- list()
descrip_tables[["kos"]] <- "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz"
descrip_tables[["pathways"]] <- "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz"
descrip_tables[["modules"]] <- "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz"

POMS_out <- list()

for (study in studies) {
 
  POMS_out[[study]] <- list()
  
  sample_info <- almeida_sample_info[which(almeida_sample_info$Study == study), ]
  
  abun_table <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = sample_info$Run)
  
  sample_info <- sample_info[which(sample_info$Run %in% colnames(abun_table)), ]
  group1_samples <- sample_info[which(sample_info$Health.state == "Diseased"), "Run"]
  group2_samples <- sample_info[which(sample_info$Health.state == "Healthy"), "Run"]
  
  for (f in names(func_tables)) {
  
    ERP002061_POMS_out[[f]] <- POMS_pipeline(abun = abun_table,
                                                func = func_tables[[f]][rownames(abun_table), ],
                                                tree = almeida_tree,
                                                group1_samples = group1_samples,
                                                group2_samples = group2_samples,
                                                ncores = 40,
                                                BSN_p_cutoff = 0.05,
                                                BSN_correction = "none",
                                                FSN_p_cutoff = 0.05,
                                                FSN_correction = "none",
                                                min_num_tips = 10,
                                                min_func_instances = 5,
                                                min_func_prop = 0.001,
                                                multinomial_correction = "BH",
                                                detailed_output = TRUE,
                                                verbose = TRUE,
                                                func_descrip_infile = descrip_tables[[f]])
    
  }
}



for (func_type in c("ko", "pathways", "modules")) {
  print("ERP002061")
  print(func_type)
  print(rownames(ERP002061_POMS_out[[func_type]]$results[which(ERP002061_POMS_out[[func_type]]$results$multinomial_corr < 0.1), ]))
  
  print("ERP012177")
  print(func_type)
  print(rownames(ERP012177_POMS_out[[func_type]]$results[which(ERP012177_POMS_out[[func_type]]$results$multinomial_corr < 0.3), ]))
  
  print("ERP003612")
  print(func_type)
  print(rownames(ERP003612_POMS_out[[func_type]]$results[which(ERP003612_POMS_out[[func_type]]$results$multinomial_corr < 0.3), ]))
}

saveRDS(object = ERP002061_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
saveRDS(object = ERP012177_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/ERP012177_POMS_out.rds")
saveRDS(object = ERP003612_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/ERP003612_POMS_out.rds")
