rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/MAGs/Almeida2019/")
source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")
devtools::load_all(path = "/home/gavin/github_repos/POMS/")

almeida_DA_out <- list()

# Read in input files.
almeida_func <- read.table(file = "functional_analyses/kegg_summary.csv", sep=",", stringsAsFactors = FALSE, quote="", comment.char = "",
                           header=TRUE, check.names = FALSE, row.names=1)

almeida_func <- data.frame(t(almeida_func), check.names = FALSE)

almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")
almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv", header=TRUE, sep="\t", check.names=FALSE,
                           row.names=1, quote="", comment.char="")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")

almeida_taxa_umgs <- read.table("taxonomy/taxonomy_umgs.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_umgs) <- almeida_taxa_umgs$MAG_ID
almeida_taxa_umgs <- almeida_taxa_umgs[, -which(colnames(almeida_taxa_umgs) %in% c("UMGS_ID", "MAG_ID"))]

almeida_taxa_hgr <- read.table("taxonomy/taxonomy_hgr.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_hgr) <- almeida_taxa_hgr$Genome
almeida_taxa_hgr <- almeida_taxa_hgr[, -which(colnames(almeida_taxa_hgr) == "Genome")]

almeida_taxa <- rbind(almeida_taxa_hgr, almeida_taxa_umgs)

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/MUSiCC_KEGG_single_copy_genes.txt", stringsAsFactors = FALSE)$V1



# ERP002061
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]

ERP002061_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)


ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP002061_almeida_func <- almeida_func[rownames(ERP002061_almeida_abun), ]

ERP002061_almeida_func_abun <- calc_func_abun(in_abun=ERP002061_almeida_abun, in_func=ERP002061_almeida_func, ncores = 70)

# Subset to samples
ERP002061_group1_samples_subset <- ERP002061_group1_samples[which(ERP002061_group1_samples %in% colnames(ERP002061_almeida_func_abun))]
ERP002061_group2_samples_subset <- ERP002061_group2_samples[which(ERP002061_group2_samples %in% colnames(ERP002061_almeida_func_abun))]
ERP002061_almeida_func_abun <- ERP002061_almeida_func_abun[, c(ERP002061_group1_samples_subset, ERP002061_group2_samples_subset)]

# Remove any rows or cols that are all missing after flooring data.
ERP002061_almeida_func_abun_floor <- floor(ERP002061_almeida_func_abun)
near_empty_rows <- which(rowSums(ERP002061_almeida_func_abun_floor) == 0)
near_empty_cols <- which(colSums(ERP002061_almeida_func_abun_floor) == 0)
if(length(near_empty_rows) > 0) { ERP002061_almeida_func_abun <- ERP002061_almeida_func_abun[-near_empty_rows, ] }
if(length(near_empty_cols) > 0) { ERP002061_almeida_func_abun <- ERP002061_almeida_func_abun[, -near_empty_cols] }

almeida_DA_out[["ERP002061"]] <- list()

almeida_DA_out[["ERP002061"]][["aldex2"]]<- run_2group_ALDEx2(in_table = ERP002061_almeida_func_abun,
                                                              group1_samples = ERP002061_group1_samples,
                                                              group2_samples = ERP002061_group2_samples)

almeida_DA_out[["ERP002061"]][["deseq2"]]<- deseq2_default_two_groups(table = ERP002061_almeida_func_abun,
                                                                   group1 = ERP002061_group1_samples,
                                                                   group2 = ERP002061_group2_samples,
                                                                   dataset_name = "ERP002061")

almeida_DA_out[["ERP002061"]][["limma.voom"]] <- limma_voom_two_group_TMM(table = ERP002061_almeida_func_abun,
                                                                       group1 = ERP002061_group1_samples,
                                                                       group2 = ERP002061_group2_samples)


almeida_DA_out[["ERP002061"]][["wilcoxon.relab"]] <- wilcoxon_2group_pvalues(intable = ERP002061_almeida_func_abun,
                                                                          group1 = ERP002061_group1_samples,
                                                                          group2 = ERP002061_group2_samples,
                                                                          convert_relab = TRUE)

ERP002061_uscg_set <- musicc_uscgs[which(musicc_uscgs %in% rownames(ERP002061_almeida_func_abun))]
ERP002061_almeida_func_abun_musicc <- data.frame(sweep(ERP002061_almeida_func_abun, 2, colMedians(as.matrix(ERP002061_almeida_func_abun[ERP002061_uscg_set, ])), `/`))
almeida_DA_out[["ERP002061"]][["wilcoxon.musicc"]] <- wilcoxon_2group_pvalues(intable = ERP002061_almeida_func_abun_musicc,
                                                                           group1 = ERP002061_group1_samples,
                                                                           group2 = ERP002061_group2_samples,
                                                                           convert_relab = FALSE)



# ERP012177
ERP012177_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP012177"), ]

ERP012177_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP012177_almeida_sample_info$Run)


ERP012177_almeida_sample_info <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Run %in% colnames(ERP012177_almeida_abun)), ]
ERP012177_group1_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP012177_group2_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP012177_almeida_func <- almeida_func[rownames(ERP012177_almeida_abun), ]

ERP012177_almeida_func_abun <- calc_func_abun(in_abun=ERP012177_almeida_abun, in_func=ERP012177_almeida_func, ncores = 70)

# Subset to samples
ERP012177_group1_samples_subset <- ERP012177_group1_samples[which(ERP012177_group1_samples %in% colnames(ERP012177_almeida_func_abun))]
ERP012177_group2_samples_subset <- ERP012177_group2_samples[which(ERP012177_group2_samples %in% colnames(ERP012177_almeida_func_abun))]
ERP012177_almeida_func_abun <- ERP012177_almeida_func_abun[, c(ERP012177_group1_samples_subset, ERP012177_group2_samples_subset)]

# Remove any rows or cols that are all missing after flooring data.
ERP012177_almeida_func_abun_floor <- floor(ERP012177_almeida_func_abun)
near_empty_rows <- which(rowSums(ERP012177_almeida_func_abun_floor) == 0)
near_empty_cols <- which(colSums(ERP012177_almeida_func_abun_floor) == 0)
if(length(near_empty_rows) > 0) { ERP012177_almeida_func_abun <- ERP012177_almeida_func_abun[-near_empty_rows, ] }
if(length(near_empty_cols) > 0) { ERP012177_almeida_func_abun <- ERP012177_almeida_func_abun[, -near_empty_cols] }

almeida_DA_out[["ERP012177"]] <- list()

almeida_DA_out[["ERP012177"]][["aldex2"]]<- run_2group_ALDEx2(in_table = ERP012177_almeida_func_abun,
                                                              group1_samples = ERP012177_group1_samples,
                                                              group2_samples = ERP012177_group2_samples)

almeida_DA_out[["ERP012177"]][["deseq2"]] <- deseq2_default_two_groups(table = ERP012177_almeida_func_abun,
                                                                   group1 = ERP012177_group1_samples,
                                                                   group2 = ERP012177_group2_samples,
                                                                   dataset_name = "ERP012177")

almeida_DA_out[["ERP012177"]][["limma.voom"]] <- limma_voom_two_group_TMM(table = ERP012177_almeida_func_abun,
                                                                      group1 = ERP012177_group1_samples,
                                                                      group2 = ERP012177_group2_samples)


almeida_DA_out[["ERP012177"]][["wilcoxon.relab"]] <- wilcoxon_2group_pvalues(intable = ERP012177_almeida_func_abun,
                                                                         group1 = ERP012177_group1_samples,
                                                                         group2 = ERP012177_group2_samples,
                                                                         convert_relab = TRUE)

ERP012177_uscg_set <- musicc_uscgs[which(musicc_uscgs %in% rownames(ERP012177_almeida_func_abun))]
ERP012177_almeida_func_abun_musicc <- data.frame(sweep(ERP012177_almeida_func_abun, 2, colMedians(as.matrix(ERP012177_almeida_func_abun[ERP012177_uscg_set, ])), `/`))
almeida_DA_out["ERP012177"]["wilcoxon.musicc"] <- wilcoxon_2group_pvalues(intable = ERP012177_almeida_func_abun_musicc,
                                                                          group1 = ERP012177_group1_samples,
                                                                          group2 = ERP012177_group2_samples,
                                                                          convert_relab = FALSE)



# ERP003612
ERP003612_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP003612"), ]

ERP003612_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP003612_almeida_sample_info$Run)


ERP003612_almeida_sample_info <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Run %in% colnames(ERP003612_almeida_abun)), ]
ERP003612_group1_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP003612_group2_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP003612_almeida_func <- almeida_func[rownames(ERP003612_almeida_abun), ]

ERP003612_almeida_func_abun <- calc_func_abun(in_abun=ERP003612_almeida_abun, in_func=ERP003612_almeida_func, ncores = 70)

# Subset to samples
ERP003612_group1_samples_subset <- ERP003612_group1_samples[which(ERP003612_group1_samples %in% colnames(ERP003612_almeida_func_abun))]
ERP003612_group2_samples_subset <- ERP003612_group2_samples[which(ERP003612_group2_samples %in% colnames(ERP003612_almeida_func_abun))]
ERP003612_almeida_func_abun <- ERP003612_almeida_func_abun[, c(ERP003612_group1_samples_subset, ERP003612_group2_samples_subset)]

# Remove any rows or cols that are all missing after flooring data.
ERP003612_almeida_func_abun_floor <- floor(ERP003612_almeida_func_abun)
near_empty_rows <- which(rowSums(ERP003612_almeida_func_abun_floor) == 0)
near_empty_cols <- which(colSums(ERP003612_almeida_func_abun_floor) == 0)
if(length(near_empty_rows) > 0) { ERP003612_almeida_func_abun <- ERP003612_almeida_func_abun[-near_empty_rows, ] }
if(length(near_empty_cols) > 0) { ERP003612_almeida_func_abun <- ERP003612_almeida_func_abun[, -near_empty_cols] }

almeida_DA_out[["ERP003612"]] <- list()

almeida_DA_out[["ERP003612"]][["aldex2"]]<- run_2group_ALDEx2(in_table = ERP003612_almeida_func_abun,
                                                              group1_samples = ERP003612_group1_samples,
                                                              group2_samples = ERP003612_group2_samples)

almeida_DA_out[["ERP003612"]][["deseq2"]] <- deseq2_default_two_groups(table = ERP003612_almeida_func_abun,
                                                                   group1 = ERP003612_group1_samples,
                                                                   group2 = ERP003612_group2_samples,
                                                                   dataset_name = "ERP003612")

almeida_DA_out[["ERP003612"]][["limma.voom"]] <- limma_voom_two_group_TMM(table = ERP003612_almeida_func_abun,
                                                                      group1 = ERP003612_group1_samples,
                                                                      group2 = ERP003612_group2_samples)


almeida_DA_out[["ERP003612"]][["wilcoxon.relab"]] <- wilcoxon_2group_pvalues(intable = ERP003612_almeida_func_abun,
                                                                         group1 = ERP003612_group1_samples,
                                                                         group2 = ERP003612_group2_samples,
                                                                         convert_relab = TRUE)

ERP003612_uscg_set <- musicc_uscgs[which(musicc_uscgs %in% rownames(ERP003612_almeida_func_abun))]
ERP003612_almeida_func_abun_musicc <- data.frame(sweep(ERP003612_almeida_func_abun, 2, colMedians(as.matrix(ERP003612_almeida_func_abun[ERP003612_uscg_set, ])), `/`))
almeida_DA_out[["ERP003612"]][["wilcoxon.musicc"]] <- wilcoxon_2group_pvalues(intable = ERP003612_almeida_func_abun_musicc,
                                                                          group1 = ERP003612_group1_samples,
                                                                          group2 = ERP003612_group2_samples,
                                                                          convert_relab = FALSE)



# Write out all DA output results for these three datasets
saveRDS(object = almeida_DA_out, file = "/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_DA_tool_output.rds")
