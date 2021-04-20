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

# Obesity #1 
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]


ERP002061_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)

ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP002061_almeida_ko <- almeida_ko[rownames(ERP002061_almeida_abun), ]
ERP002061_almeida_pathways <- almeida_pathways[rownames(ERP002061_almeida_abun), ]
ERP002061_almeida_modules <- almeida_modules[rownames(ERP002061_almeida_abun), ]

ERP002061_POMS_out <- list()

if (length(which(colSums(ERP002061_almeida_ko > 0) > 0.85 * nrow(ERP002061_almeida_ko))) > 0) {
  ERP002061_almeida_ko_subset <- ERP002061_almeida_ko[, -which(colSums(ERP002061_almeida_ko > 0) > 0.85 * nrow(ERP002061_almeida_ko))]
} else {
  ERP002061_almeida_ko_subset <- ERP002061_almeida_ko 
}

if (length(which(colSums(ERP002061_almeida_ko > 0) < 0.15 * nrow(ERP002061_almeida_ko))) > 0) {
  ERP002061_almeida_ko_subset <- ERP002061_almeida_ko[, -which(colSums(ERP002061_almeida_ko > 0) < 0.15 * nrow(ERP002061_almeida_ko))]
}

ERP002061_POMS_out[["ko"]] <- POMS_pipeline(abun = ERP002061_almeida_abun,
                                            func = ERP002061_almeida_ko,
                                            phylogeny = almeida_tree,
                                            group1_samples = ERP002061_group1_samples,
                                            group2_samples = ERP002061_group2_samples,
                                            ncores = 10,
                                            balance_p_cutoff = 0.05,
                                            balance_correction = "none",
                                            function_p_cutoff = 0.05,
                                            function_correction = "none",
                                            min_num_tips = 10,
                                            min_func_instances = 1,
                                            min_func_prop = 0.001,
                                            run_multinomial_test = TRUE,
                                            multinomial_correction = "BH",
                                            detailed_output = TRUE,
                                            verbose = TRUE,
                                            func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz")


if (length(which(colSums(ERP002061_almeida_pathways > 0) > 0.85 * nrow(ERP002061_almeida_pathways))) > 0) {
  ERP002061_almeida_pathways_subset <- ERP002061_almeida_pathways[, -which(colSums(ERP002061_almeida_pathways > 0) > 0.85 * nrow(ERP002061_almeida_pathways))]
} else {
  ERP002061_almeida_pathways_subset <- ERP002061_almeida_pathways 
}

if (length(which(colSums(ERP002061_almeida_pathways > 0) < 0.15 * nrow(ERP002061_almeida_pathways))) > 0) {
  ERP002061_almeida_pathways_subset <- ERP002061_almeida_pathways[, -which(colSums(ERP002061_almeida_pathways > 0) < 0.15 * nrow(ERP002061_almeida_pathways))]
}

ERP002061_POMS_out[["pathways"]] <- POMS_pipeline(abun = ERP002061_almeida_abun,
                                                 func = ERP002061_almeida_pathways_subset,
                                                 phylogeny = almeida_tree,
                                                 group1_samples = ERP002061_group1_samples,
                                                 group2_samples = ERP002061_group2_samples,
                                                 ncores = 10,
                                                 balance_p_cutoff = 0.05,
                                                 balance_correction = "none",
                                                 function_p_cutoff = 0.05,
                                                 function_correction = "none",
                                                 min_num_tips = 10,
                                                 min_func_instances = 1,
                                                 min_func_prop = 0.001,
                                                 run_multinomial_test = TRUE,
                                                 multinomial_correction = "BH",
                                                 detailed_output = TRUE,
                                                 verbose = TRUE,
                                                 func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz")


if (length(which(colSums(ERP002061_almeida_modules > 0) > 0.85 * nrow(ERP002061_almeida_modules))) > 0) {
  ERP002061_almeida_modules_subset <- ERP002061_almeida_modules[, -which(colSums(ERP002061_almeida_modules > 0) > 0.85 * nrow(ERP002061_almeida_modules))]
} else {
  ERP002061_almeida_modules_subset <- ERP002061_almeida_modules 
}

if (length(which(colSums(ERP002061_almeida_modules > 0) < 0.15 * nrow(ERP002061_almeida_modules))) > 0) {
  ERP002061_almeida_modules_subset <- ERP002061_almeida_modules[, -which(colSums(ERP002061_almeida_modules > 0) < 0.15 * nrow(ERP002061_almeida_modules))]
}

ERP002061_POMS_out[["modules"]] <- POMS_pipeline(abun = ERP002061_almeida_abun,
                                                func = ERP002061_almeida_modules_subset,
                                                phylogeny = almeida_tree,
                                                group1_samples = ERP002061_group1_samples,
                                                group2_samples = ERP002061_group2_samples,
                                                ncores = 10,
                                                balance_p_cutoff = 0.05,
                                                balance_correction = "none",
                                                function_p_cutoff = 0.05,
                                                function_correction = "none",
                                                min_num_tips = 10,
                                                min_func_instances = 100,
                                                min_func_prop = 0.001,
                                                run_multinomial_test = TRUE,
                                                multinomial_correction = "BH",
                                                detailed_output = TRUE,
                                                verbose = TRUE,
                                                func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz")


# CRC
ERP012177_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP012177"), ]


ERP012177_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP012177_almeida_sample_info$Run)

ERP012177_almeida_sample_info <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Run %in% colnames(ERP012177_almeida_abun)), ]
ERP012177_group1_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP012177_group2_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP012177_almeida_ko <- almeida_ko[rownames(ERP012177_almeida_abun), ]
ERP012177_almeida_pathways <- almeida_pathways[rownames(ERP012177_almeida_abun), ]
ERP012177_almeida_modules <- almeida_modules[rownames(ERP012177_almeida_abun), ]

ERP012177_POMS_out <- list()

if (length(which(colSums(ERP012177_almeida_ko > 0) > 0.85 * nrow(ERP012177_almeida_ko))) > 0) {
  ERP012177_almeida_ko_subset <- ERP012177_almeida_ko[, -which(colSums(ERP012177_almeida_ko > 0) > 0.85 * nrow(ERP012177_almeida_ko))]
} else {
  ERP012177_almeida_ko_subset <- ERP012177_almeida_ko 
}

if (length(which(colSums(ERP012177_almeida_ko > 0) < 0.15 * nrow(ERP012177_almeida_ko))) > 0) {
  ERP012177_almeida_ko_subset <- ERP012177_almeida_ko[, -which(colSums(ERP012177_almeida_ko > 0) < 0.15 * nrow(ERP012177_almeida_ko))]
}

ERP012177_POMS_out[["ko"]] <- POMS_pipeline(abun = ERP012177_almeida_abun,
                                            func = ERP012177_almeida_ko,
                                            phylogeny = almeida_tree,
                                            group1_samples = ERP012177_group1_samples,
                                            group2_samples = ERP012177_group2_samples,
                                            ncores = 10,
                                            balance_p_cutoff = 0.05,
                                            balance_correction = "none",
                                            function_p_cutoff = 0.05,
                                            function_correction = "none",
                                            min_num_tips = 10,
                                            min_func_instances = 1,
                                            min_func_prop = 0.001,
                                            run_multinomial_test = TRUE,
                                            multinomial_correction = "BH",
                                            detailed_output = TRUE,
                                            verbose = TRUE,
                                            func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz")


if (length(which(colSums(ERP012177_almeida_pathways > 0) > 0.85 * nrow(ERP012177_almeida_pathways))) > 0) {
  ERP012177_almeida_pathways_subset <- ERP012177_almeida_pathways[, -which(colSums(ERP012177_almeida_pathways > 0) > 0.85 * nrow(ERP012177_almeida_pathways))]
} else {
  ERP012177_almeida_pathways_subset <- ERP012177_almeida_pathways 
}

if (length(which(colSums(ERP012177_almeida_pathways > 0) < 0.15 * nrow(ERP012177_almeida_pathways))) > 0) {
  ERP012177_almeida_pathways_subset <- ERP012177_almeida_pathways[, -which(colSums(ERP012177_almeida_pathways > 0) < 0.15 * nrow(ERP012177_almeida_pathways))]
}

ERP012177_POMS_out[["pathways"]] <- POMS_pipeline(abun = ERP012177_almeida_abun,
                                                  func = ERP012177_almeida_pathways_subset,
                                                  phylogeny = almeida_tree,
                                                  group1_samples = ERP012177_group1_samples,
                                                  group2_samples = ERP012177_group2_samples,
                                                  ncores = 10,
                                                  balance_p_cutoff = 0.05,
                                                  balance_correction = "none",
                                                  function_p_cutoff = 0.05,
                                                  function_correction = "none",
                                                  min_num_tips = 10,
                                                  min_func_instances = 1,
                                                  min_func_prop = 0.001,
                                                  run_multinomial_test = TRUE,
                                                  multinomial_correction = "BH",
                                                  detailed_output = TRUE,
                                                  verbose = TRUE,
                                                  func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz")


if (length(which(colSums(ERP012177_almeida_modules > 0) > 0.85 * nrow(ERP012177_almeida_modules))) > 0) {
  ERP012177_almeida_modules_subset <- ERP012177_almeida_modules[, -which(colSums(ERP012177_almeida_modules > 0) > 0.85 * nrow(ERP012177_almeida_modules))]
} else {
  ERP012177_almeida_modules_subset <- ERP012177_almeida_modules 
}

if (length(which(colSums(ERP012177_almeida_modules > 0) < 0.15 * nrow(ERP012177_almeida_modules))) > 0) {
  ERP012177_almeida_modules_subset <- ERP012177_almeida_modules[, -which(colSums(ERP012177_almeida_modules > 0) < 0.15 * nrow(ERP012177_almeida_modules))]
}

ERP012177_POMS_out[["modules"]] <- POMS_pipeline(abun = ERP012177_almeida_abun,
                                                 func = ERP012177_almeida_modules_subset,
                                                 phylogeny = almeida_tree,
                                                 group1_samples = ERP012177_group1_samples,
                                                 group2_samples = ERP012177_group2_samples,
                                                 ncores = 10,
                                                 balance_p_cutoff = 0.05,
                                                 balance_correction = "none",
                                                 function_p_cutoff = 0.05,
                                                 function_correction = "none",
                                                 min_num_tips = 10,
                                                 min_func_instances = 100,
                                                 min_func_prop = 0.001,
                                                 run_multinomial_test = TRUE,
                                                 multinomial_correction = "BH",
                                                 detailed_output = TRUE,
                                                 verbose = TRUE,
                                                 func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz")

# Obesity number 2 - ERP003612
ERP003612_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP003612"), ]


ERP003612_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP003612_almeida_sample_info$Run)

ERP003612_almeida_sample_info <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Run %in% colnames(ERP003612_almeida_abun)), ]
ERP003612_group1_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP003612_group2_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP003612_almeida_ko <- almeida_ko[rownames(ERP003612_almeida_abun), ]
ERP003612_almeida_pathways <- almeida_pathways[rownames(ERP003612_almeida_abun), ]
ERP003612_almeida_modules <- almeida_modules[rownames(ERP003612_almeida_abun), ]

ERP003612_POMS_out <- list()

if (length(which(colSums(ERP003612_almeida_ko > 0) > 0.85 * nrow(ERP003612_almeida_ko))) > 0) {
  ERP003612_almeida_ko_subset <- ERP003612_almeida_ko[, -which(colSums(ERP003612_almeida_ko > 0) > 0.85 * nrow(ERP003612_almeida_ko))]
} else {
  ERP003612_almeida_ko_subset <- ERP003612_almeida_ko 
}

if (length(which(colSums(ERP003612_almeida_ko > 0) < 0.15 * nrow(ERP003612_almeida_ko))) > 0) {
  ERP003612_almeida_ko_subset <- ERP003612_almeida_ko[, -which(colSums(ERP003612_almeida_ko > 0) < 0.15 * nrow(ERP003612_almeida_ko))]
}

ERP003612_POMS_out[["ko"]] <- POMS_pipeline(abun = ERP003612_almeida_abun,
                                            func = ERP003612_almeida_ko,
                                            phylogeny = almeida_tree,
                                            group1_samples = ERP003612_group1_samples,
                                            group2_samples = ERP003612_group2_samples,
                                            ncores = 10,
                                            balance_p_cutoff = 0.05,
                                            balance_correction = "none",
                                            function_p_cutoff = 0.05,
                                            function_correction = "none",
                                            min_num_tips = 10,
                                            min_func_instances = 1,
                                            min_func_prop = 0.001,
                                            run_multinomial_test = TRUE,
                                            multinomial_correction = "BH",
                                            detailed_output = TRUE,
                                            verbose = TRUE,
                                            func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_ko_descrip.tsv.gz")


if (length(which(colSums(ERP003612_almeida_pathways > 0) > 0.85 * nrow(ERP003612_almeida_pathways))) > 0) {
  ERP003612_almeida_pathways_subset <- ERP003612_almeida_pathways[, -which(colSums(ERP003612_almeida_pathways > 0) > 0.85 * nrow(ERP003612_almeida_pathways))]
} else {
  ERP003612_almeida_pathways_subset <- ERP003612_almeida_pathways 
}

if (length(which(colSums(ERP003612_almeida_pathways > 0) < 0.15 * nrow(ERP003612_almeida_pathways))) > 0) {
  ERP003612_almeida_pathways_subset <- ERP003612_almeida_pathways[, -which(colSums(ERP003612_almeida_pathways > 0) < 0.15 * nrow(ERP003612_almeida_pathways))]
}

ERP003612_POMS_out[["pathways"]] <- POMS_pipeline(abun = ERP003612_almeida_abun,
                                                  func = ERP003612_almeida_pathways_subset,
                                                  phylogeny = almeida_tree,
                                                  group1_samples = ERP003612_group1_samples,
                                                  group2_samples = ERP003612_group2_samples,
                                                  ncores = 10,
                                                  balance_p_cutoff = 0.05,
                                                  balance_correction = "none",
                                                  function_p_cutoff = 0.05,
                                                  function_correction = "none",
                                                  min_num_tips = 10,
                                                  min_func_instances = 1,
                                                  min_func_prop = 0.001,
                                                  run_multinomial_test = TRUE,
                                                  multinomial_correction = "BH",
                                                  detailed_output = TRUE,
                                                  verbose = TRUE,
                                                  func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_pathway_descrip.tsv.gz")


if (length(which(colSums(ERP003612_almeida_modules > 0) > 0.85 * nrow(ERP003612_almeida_modules))) > 0) {
  ERP003612_almeida_modules_subset <- ERP003612_almeida_modules[, -which(colSums(ERP003612_almeida_modules > 0) > 0.85 * nrow(ERP003612_almeida_modules))]
} else {
  ERP003612_almeida_modules_subset <- ERP003612_almeida_modules 
}

if (length(which(colSums(ERP003612_almeida_modules > 0) < 0.15 * nrow(ERP003612_almeida_modules))) > 0) {
  ERP003612_almeida_modules_subset <- ERP003612_almeida_modules[, -which(colSums(ERP003612_almeida_modules > 0) < 0.15 * nrow(ERP003612_almeida_modules))]
}

ERP003612_POMS_out[["modules"]] <- POMS_pipeline(abun = ERP003612_almeida_abun,
                                                 func = ERP003612_almeida_modules_subset,
                                                 phylogeny = almeida_tree,
                                                 group1_samples = ERP003612_group1_samples,
                                                 group2_samples = ERP003612_group2_samples,
                                                 ncores = 10,
                                                 balance_p_cutoff = 0.05,
                                                 balance_correction = "none",
                                                 function_p_cutoff = 0.05,
                                                 function_correction = "none",
                                                 min_num_tips = 10,
                                                 min_func_instances = 100,
                                                 min_func_prop = 0.001,
                                                 run_multinomial_test = TRUE,
                                                 multinomial_correction = "BH",
                                                 detailed_output = TRUE,
                                                 verbose = TRUE,
                                                 func_descrip_infile = "/home/gavin/github_repos/POMS_manuscript/data/KEGG_mappings/prepped/2021_04_12_KEGG_module_descrip.tsv.gz")


for (func_type in c("ko", "pathways", "modules")) {
  print(func_type)
  print(ERP002061_POMS_out[[func_type]]$df[which(ERP002061_POMS_out[[func_type]]$df$multinomial_corr < 0.25), ])
}

saveRDS(object = ERP002061_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
saveRDS(object = ERP012177_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/ERP012177_POMS_out.rds")
saveRDS(object = ERP003612_POMS_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/ERP003612_POMS_out.rds")
