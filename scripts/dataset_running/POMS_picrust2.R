rm(list=ls(all.names=TRUE))

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin//projects/ml_picrust2_application")

library(ape)
library(parallel)
library(phangorn)
library(stringr)

ob_datasets <- c("ob_ross", "ob_turnbaugh", "ob_zhu", "ob_goodrich")
group1_label <- "H"
group2_label <- "OB"
comparison_column <- "DiseaseState"

for(dataset in ob_datasets) {
  
  dataset_tree_path <- paste("MicrobiomeHD_picrust2_out/", dataset, ".picrust2_out/out.tre", sep="")
  dataset_KO_path <- paste("MicrobiomeHD_picrust2_out/", dataset, ".picrust2_out/KO_predicted.tsv.gz", sep="")
  dataset_abun_path <- paste("MicrobiomeHD_github_repo/microbiomeHD/data/clean_tables_csv/", dataset, ".out_table.clean.biom.tsv", sep="")
  dataset_meta_path <- paste("MicrobiomeHD_github_repo/microbiomeHD/data/clean_tables_csv/", dataset, ".metadata.clean.tsv", sep="")
  dataset_taxa_path <- paste("MicrobiomeHD_github_repo/microbiomeHD/data/clean_tables_csv/", dataset, ".tax_table.tsv", sep="")
  
  dataset_tree <- read.tree(dataset_tree_path)
  dataset_KO <- read.table(dataset_KO_path, sep="\t", row.names=1, header=TRUE, check.names=FALSE)
  dataset_abun <- read.table(dataset_abun_path, sep="\t", row.names=1, header=TRUE, comment.char="", check.names=FALSE)
  dataset_meta <- read.table(dataset_meta_path, sep="\t", row.names=1, header=TRUE, comment.char="", check.names=FALSE, stringsAsFactors=FALSE, quote="")
  dataset_taxa <- read.table(dataset_taxa_path, sep="\t", header=TRUE, comment.char="", check.names=FALSE, stringsAsFactors=FALSE, quote="")
  
  dataset_meta <- dataset_meta[which(rownames(dataset_meta) %in% colnames(dataset_abun)), ]
  
  group1_samples <- rownames(dataset_meta)[which(dataset_meta[, comparison_column] == group1_label)]
  group2_samples <- rownames(dataset_meta)[which(dataset_meta[, comparison_column] == group2_label)]
  
  dataset_abun <- dataset_abun[, c(group1_samples, group2_samples)]
  if(length(which(rowSums(dataset_abun) == 0)) > 0) {
    dataset_abun <- dataset_abun[-which(rowSums(dataset_abun) == 0), ]
  }
  
  dataset_taxa <- dataset_taxa[rownames(dataset_abun), ]
  
  dataset_KO <- dataset_KO[rownames(dataset_abun), ]
  if(length(which(rowSums(dataset_KO) == 0)) > 0) {
    dataset_KO <- dataset_KO[-which(rowSums(dataset_KO) == 0), ]
  }
  
  dataset_POMS_out <- two_group_balance_tree_pipeline(abun = dataset_abun,
                                                      func = dataset_KO,
                                                      phylogeny = dataset_tree,
                                                      taxa = dataset_taxa,
                                                      group1_samples = group1_samples,
                                                      group2_samples = group2_samples,
                                                      ncores=20,
                                                      skip_node_dist = TRUE,
                                                      verbose = TRUE,
                                                      balance_correction = "none")
  
  
  dataset_POMS_out$group1_label <- group1_label
  dataset_POMS_out$group2_label <- group2_label
  
  saveRDS(object = dataset_POMS_out, file = paste("/home/gavin/projects/POMS/dataset_POMS_out/picrust2_out/", dataset, ".rds", sep=""))

}

