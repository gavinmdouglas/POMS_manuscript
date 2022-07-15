rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

almeida_DA_out <- list()

# Read in input files.
almeida_func <- list()

almeida_func[["ko"]] <- read.table(file = "functional_analyses/kegg_summary.csv.gz",
                                    sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
almeida_func[["ko"]] <- data.frame(t(almeida_func[["ko"]]), check.names = FALSE)

almeida_func[["pathways"]] <- read.table(file = "functional_analyses/modified/kegg_pathways/path_abun_unstrat.tsv.gz",
                                          sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
almeida_func[["pathways"]] <- data.frame(t(almeida_func[["pathways"]]), check.names = FALSE)

almeida_func[["modules"]] <- read.table(file = "functional_analyses/modified/kegg_modules/path_abun_unstrat.tsv.gz",
                                         sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
almeida_func[["modules"]] <- data.frame(t(almeida_func[["modules"]]), check.names = FALSE)

almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")

almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt.gz",
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

almeida_taxa <- readRDS("taxonomy/taxa_table.rds")

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1


# ERP002061
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]
ERP002061_almeida_abun <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)
ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

almeida_DA_out[["ERP002061"]] <- list()

for (func_type in names(almeida_func)) {
  
  ERP002061_almeida_func <- almeida_func[[func_type]][rownames(ERP002061_almeida_abun), ]

  ERP002061_almeida_func <- POMS::filter_rare_table_cols(in_tab = ERP002061_almeida_func,
                                                         min_nonzero_count = 5,
                                                         min_nonzero_prop = 0.001)
  
  ERP002061_almeida_func_abun <- calc_func_abun_crossprod(in_abun = ERP002061_almeida_abun,
                                                          in_func = ERP002061_almeida_func)

  # Subset to samples
  ERP002061_group1_samples_subset <- ERP002061_group1_samples[which(ERP002061_group1_samples %in% colnames(ERP002061_almeida_func_abun))]
  ERP002061_group2_samples_subset <- ERP002061_group2_samples[which(ERP002061_group2_samples %in% colnames(ERP002061_almeida_func_abun))]
  ERP002061_almeida_func_abun <- ERP002061_almeida_func_abun[, c(ERP002061_group1_samples_subset, ERP002061_group2_samples_subset)]
  
  
  if (func_type == "ko") { 

    almeida_DA_out[["ERP002061"]][[func_type]] <- run_alt.tools(func_abun_table = ERP002061_almeida_func_abun,
                                                                group1_samples = ERP002061_group1_samples_subset,
                                                                group2_samples = ERP002061_group2_samples_subset,
                                                                USCGs = musicc_uscgs,
                                                                tools_to_run = c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc"))
  } else {

    almeida_DA_out[["ERP002061"]][[func_type]] <- run_alt.tools(func_abun_table = ERP002061_almeida_func_abun,
                                                                group1_samples = ERP002061_group1_samples_subset,
                                                                group2_samples = ERP002061_group2_samples_subset,
                                                                tools_to_run = c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab"))
  }

}


# ERP012177
ERP012177_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP012177"), ]
ERP012177_almeida_abun <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = ERP012177_almeida_sample_info$Run)
ERP012177_almeida_sample_info <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Run %in% colnames(ERP012177_almeida_abun)), ]
ERP012177_group1_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP012177_group2_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Healthy"), "Run"]

almeida_DA_out[["ERP012177"]] <- list()

for (func_type in names(almeida_func)) {
  
  ERP012177_almeida_func <- almeida_func[[func_type]][rownames(ERP012177_almeida_abun), ]
  
  ERP012177_almeida_func <- POMS::filter_rare_table_cols(in_tab = ERP012177_almeida_func,
                                                         min_nonzero_count = 5,
                                                         min_nonzero_prop = 0.001)
  
  ERP012177_almeida_func_abun <- calc_func_abun_crossprod(in_abun = ERP012177_almeida_abun,
                                                          in_func = ERP012177_almeida_func)
  
  # Subset to samples
  ERP012177_group1_samples_subset <- ERP012177_group1_samples[which(ERP012177_group1_samples %in% colnames(ERP012177_almeida_func_abun))]
  ERP012177_group2_samples_subset <- ERP012177_group2_samples[which(ERP012177_group2_samples %in% colnames(ERP012177_almeida_func_abun))]
  ERP012177_almeida_func_abun <- ERP012177_almeida_func_abun[, c(ERP012177_group1_samples_subset, ERP012177_group2_samples_subset)]
  
  if (func_type == "ko") { 
    
    almeida_DA_out[["ERP012177"]][[func_type]] <- run_alt.tools(func_abun_table = ERP012177_almeida_func_abun,
                                                                group1_samples = ERP012177_group1_samples_subset,
                                                                group2_samples = ERP012177_group2_samples_subset,
                                                                USCGs = musicc_uscgs,
                                                                tools_to_run = c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc"))
  } else {
    
    almeida_DA_out[["ERP012177"]][[func_type]] <- run_alt.tools(func_abun_table = ERP012177_almeida_func_abun,
                                                                group1_samples = ERP012177_group1_samples_subset,
                                                                group2_samples = ERP012177_group2_samples_subset,
                                                                tools_to_run = c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab"))
  }
  
}

# ERP003612
ERP003612_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP003612"), ]
ERP003612_almeida_abun <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = ERP003612_almeida_sample_info$Run)
ERP003612_almeida_sample_info <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Run %in% colnames(ERP003612_almeida_abun)), ]
ERP003612_group1_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP003612_group2_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Healthy"), "Run"]



almeida_DA_out[["ERP003612"]] <- list()

for (func_type in names(almeida_func)) {
  
  ERP003612_almeida_func <- almeida_func[[func_type]][rownames(ERP003612_almeida_abun), ]
  
  ERP003612_almeida_func <- POMS::filter_rare_table_cols(in_tab = ERP003612_almeida_func,
                                                         min_nonzero_count = 5,
                                                         min_nonzero_prop = 0.001)
  
  ERP003612_almeida_func_abun <- calc_func_abun_crossprod(in_abun = ERP003612_almeida_abun,
                                                          in_func = ERP003612_almeida_func)
  
  # Subset to samples
  ERP003612_group1_samples_subset <- ERP003612_group1_samples[which(ERP003612_group1_samples %in% colnames(ERP003612_almeida_func_abun))]
  ERP003612_group2_samples_subset <- ERP003612_group2_samples[which(ERP003612_group2_samples %in% colnames(ERP003612_almeida_func_abun))]
  ERP003612_almeida_func_abun <- ERP003612_almeida_func_abun[, c(ERP003612_group1_samples_subset, ERP003612_group2_samples_subset)]
  
  if (func_type == "ko") { 
    
    almeida_DA_out[["ERP003612"]][[func_type]] <- run_alt.tools(func_abun_table = ERP003612_almeida_func_abun,
                                                                group1_samples = ERP003612_group1_samples_subset,
                                                                group2_samples = ERP003612_group2_samples_subset,
                                                                USCGs = musicc_uscgs,
                                                                tools_to_run = c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc"))
  } else {
    
    almeida_DA_out[["ERP003612"]][[func_type]] <- run_alt.tools(func_abun_table = ERP003612_almeida_func_abun,
                                                                group1_samples = ERP003612_group1_samples_subset,
                                                                group2_samples = ERP003612_group2_samples_subset,
                                                                tools_to_run = c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab"))
  }
  
}

# Write out all DA output results for these three datasets
saveRDS(object = almeida_DA_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_DA_tool_output.rds")
