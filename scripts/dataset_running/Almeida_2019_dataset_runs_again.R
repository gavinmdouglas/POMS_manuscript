### Explore what the distribution of pos / neg enrichments looks like when the same number of balances are set to be arbitrarily significant

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/MAGs/Almeida2019/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

### First re-generate results for large obesity dataset, as before

# Read in input files.
almeida_func <- read.table(file = "functional_analyses/kegg_summary.csv", sep=",", stringsAsFactors = FALSE, quote="", comment.char = "",
                           header=TRUE, check.names = FALSE, row.names=1)

almeida_func <- data.frame(t(almeida_func), check.names = FALSE)

almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")
almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv", header=TRUE, sep="\t", check.names=FALSE,
                           row.names=1, quote="", comment.char="")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")

# Obesity #1 
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]


ERP002061_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)

ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP002061_almeida_func <- almeida_func[rownames(ERP002061_almeida_abun), ]

ERP002061_almeida_out <- two_group_balance_tree_pipeline(abun=ERP002061_almeida_abun,
                                                         func=ERP002061_almeida_func,
                                                         phylogeny=almeida_tree,
                                                         group1_samples = ERP002061_group1_samples,
                                                         group2_samples = ERP002061_group2_samples,
                                                         ncores=50,
                                                         detailed_output = TRUE)


# CRC
ERP012177_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP012177"), ]


ERP012177_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP012177_almeida_sample_info$Run)

ERP012177_almeida_sample_info <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Run %in% colnames(ERP012177_almeida_abun)), ]
ERP012177_group1_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP012177_group2_samples <- ERP012177_almeida_sample_info[which(ERP012177_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP012177_almeida_func <- almeida_func[rownames(ERP012177_almeida_abun), ]

ERP012177_almeida_out <- two_group_balance_tree_pipeline(abun=ERP012177_almeida_abun,
                                                         func=ERP012177_almeida_func,
                                                         phylogeny=almeida_tree,
                                                         group1_samples = ERP012177_group1_samples,
                                                         group2_samples = ERP012177_group2_samples,
                                                         ncores=50,
                                                         detailed_output = TRUE)



# Obesity number 2 - ERP003612
ERP003612_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP003612"), ]


ERP003612_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP003612_almeida_sample_info$Run)

ERP003612_almeida_sample_info <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Run %in% colnames(ERP003612_almeida_abun)), ]
ERP003612_group1_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP003612_group2_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP003612_almeida_func <- almeida_func[rownames(ERP003612_almeida_abun), ]

ERP003612_almeida_out <- two_group_balance_tree_pipeline(abun=ERP003612_almeida_abun,
                                                         func=ERP003612_almeida_func,
                                                         phylogeny=almeida_tree,
                                                         group1_samples = ERP003612_group1_samples,
                                                         group2_samples = ERP003612_group2_samples,
                                                         ncores=50,
                                                         detailed_output = TRUE)


saveRDS(object = ERP002061_almeida_out, file="/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_POMS_output/ERP002061_POMS_out_rerun.rds")
saveRDS(object = ERP012177_almeida_out, file="/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_POMS_output/ERP012177_POMS_out_rerun.rds")
saveRDS(object = ERP003612_almeida_out, file="/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_POMS_output/ERP003612_POMS_out_rerun.rds")
