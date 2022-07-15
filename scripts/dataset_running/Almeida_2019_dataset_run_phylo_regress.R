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

regress_out <- list()

for (study in studies) {
 
  regress_out[[study]] <- list()
  
  sample_info <- almeida_sample_info[which(almeida_sample_info$Study == study), ]
  
  abun_table <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = sample_info$Run)
  
  sample_info <- sample_info[which(sample_info$Run %in% colnames(abun_table)), ]
  group1_samples <- sample_info[which(sample_info$Health.state == "Diseased"), "Run"]
  group2_samples <- sample_info[which(sample_info$Health.state == "Healthy"), "Run"]
  
  prepped_tree <- POMS::prep_tree(phy = almeida_tree, tips2keep = rownames(abun_table))
  
  metadata <- data.frame(samp = c(group1_samples,
                                  group2_samples),
                         group = c(rep("group1", length(group1_samples)),
                                   rep("group2", length(group2_samples))))
  
  taxa_specificity <- specificity_scores(abun_table = abun_table,
                                         meta_table = metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[prepped_tree$tip.label]
  
  for (f in names(func_tables)) {
  
    filt_func_table <-  POMS::filter_rare_table_cols(in_tab = func_tables[[f]][rownames(abun_table), ],
                                                     min_nonzero_count = 5,
                                                     min_nonzero_prop = 0.001)
    
    regress_out[[study]][[f]] <- genome_content_phylo_regress(y = taxa_specificity,
                                                              func =  filt_func_table,
                                                              in_tree = prepped_tree,
                                                              ncores = 20,
                                                              model_type = "BM")
    
    regress_out[[study]][[f]]$BH <- p.adjust(regress_out[[study]][[f]]$p, "BH")

  }

}


for(study in names(regress_out)) {
  print("=========")
  print(study)
  print("   ")
  
  for (func_type in names(func_tables)) {
    print(func_type)
    print(length(which(regress_out[[study]][[func_type]]$BH < 0.2)))
    print("     ")
  }
}

saveRDS(object = regress_out, file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_regress_specificity_output/combined_output.rds")
