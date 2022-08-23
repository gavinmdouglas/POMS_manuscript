rm(list = ls(all.names = TRUE))

# Compute Faiths PD for all sig. hits in POMS and Phylo. Regress. output.

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


almeida_regress_out <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_regress_specificity_output/combined_output.rds")
almeida_POMS_out <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/combined_output.rds")

studies <- c("ERP002061", "ERP012177", "ERP003612")

func_tables <- list()
func_tables[["kos"]] <- almeida_ko
func_tables[["pathways"]] <- almeida_pathways
func_tables[["modules"]] <- almeida_modules

almeida_output_raw <- list()
almeida_output_edges1_raw <- list()

counter <- 1
for (study in studies) {
  
  sample_info <- almeida_sample_info[which(almeida_sample_info$Study == study), ]
  
  abun_table <- subset_by_col_and_filt(in_tab = almeida_abun, col2keep = sample_info$Run)
  
  for (func_type in names(func_tables)) {

    filt_func_table <-  POMS::filter_rare_table_cols(in_tab = func_tables[[func_type]][rownames(abun_table), ],
                                                     min_nonzero_count = 5,
                                                     min_nonzero_prop = 0.001)
    
    POMS_hits <- rownames(almeida_POMS_out[[study]][[func_type]]$results)[which(almeida_POMS_out[[study]][[func_type]]$results$multinomial_corr < 0.25)]
    regress_hits <- rownames(almeida_regress_out[[study]][[func_type]])[which(almeida_regress_out[[study]][[func_type]]$BH < 0.25)]
    all_hits <- unique(c(POMS_hits, regress_hits))
    
    input_tree <- POMS::prep_tree(phy = almeida_tree, rownames(filt_func_table))
    pd_output <- picante::pd(samp = t(filt_func_table[, all_hits]), tree = input_tree)

    input_tree_edge1 <- input_tree
    input_tree_edge1$edge.length <- rep(1, length(input_tree_edge1$edge.length))
    pd_output_edge1 <- picante::pd(samp = t(filt_func_table[, all_hits]), tree = input_tree_edge1)
    
    if (study == "ERP002061") {
      study_clean <- "Obesity 1"
    } else if (study == "ERP003612") {
      study_clean <- "Obesity 2"
    } else if (study == "ERP012177") {
      study_clean <- "Colorectal cancer"
    } else {
      stop("Error!")
    }
    
    almeida_output_raw[[counter]] <- data.frame(Dataset = study_clean,
                                                Dataset_id = study,
                                                Approach = c(rep("POMS", length(POMS_hits)),
                                                              rep("Phylo. regress. (specificity)", length(regress_hits))),
                                                Datatype = func_type,
                                                Func = c(POMS_hits, regress_hits),
                                                PD = pd_output[c(POMS_hits, regress_hits), "PD"],
                                                num_encoders = pd_output[c(POMS_hits, regress_hits), "SR"])
    
    almeida_output_edges1_raw[[counter]] <- data.frame(Dataset = study_clean,
                                                       Dataset_id = study,
                                                        Approach = c(rep("POMS", length(POMS_hits)),
                                                                     rep("Phylo. regress. (specificity)", length(regress_hits))),
                                                        Datatype = func_type,
                                                        Func = c(POMS_hits, regress_hits),
                                                        PD = pd_output_edge1[c(POMS_hits, regress_hits), "PD"],
                                                        num_encoders = pd_output_edge1[c(POMS_hits, regress_hits), "SR"])
    
    counter <- counter + 1

  }
}

almeida_output <- do.call(rbind, almeida_output_raw)
almeida_output_edges1 <- do.call(rbind, almeida_output_edges1_raw)

write.table(x = almeida_output,
            file = "/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits/case.control_Almeida.2019.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(x = almeida_output_edges1,
            file = "/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits/case.control_Almeida.2019_edges1.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Same thing, but for TARA oceans.
rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(picante)
library(parallel)

# First compute Faith's phylogenetic diversity for all tested functions (based on the subset of MAGs used).
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

TARA_tree <- read.tree(file = "GToTree_output/GToTree_output_modified.tre")


TARA_ko <- POMS::filter_rare_table_cols(in_tab = TARA_ko, min_nonzero_count = 5, min_nonzero_prop = 0.001)
TARA_pathways <- POMS::filter_rare_table_cols(in_tab = TARA_pathways, min_nonzero_count = 5, min_nonzero_prop = 0.001)
TARA_modules <- POMS::filter_rare_table_cols(in_tab = TARA_modules, min_nonzero_count = 5, min_nonzero_prop = 0.001)

TARA_func <- list()
TARA_func[["ko"]] <- TARA_ko
TARA_func[["pathway"]] <- TARA_pathways
TARA_func[["module"]] <- TARA_modules


TARA_abun <- read.table(file = "NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/modified/TARA_abundance_min_mean_coverage1.tsv",
                        header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

TARA_sample_info <- read.table("Table_S1_sample_info.txt",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)


# Identify significant nodes outside of main POMS pipeline, based on Spearman correlations with sample metadata.
TARA_tree <- POMS::prep_tree(phy = TARA_tree, tips2keep = TARA_tree$tip.label)

TARA_POMS_out <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/results/TARA_POMS_out.rds")
TARA_regress_out <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/results/TARA_regress_out.rds")

TARA_output_raw <- list()
TARA_output_edges1_raw <- list()

counter <- 1

for (env_factor in c("Mean_Salinity", "PO4")) {

  for (func_type in names(TARA_func)) {
    
    filt_func_table <-  POMS::filter_rare_table_cols(in_tab = TARA_func[[func_type]][rownames(TARA_abun), ],
                                                     min_nonzero_count = 5,
                                                     min_nonzero_prop = 0.001)
    
    POMS_hits <- rownames(TARA_POMS_out[[env_factor]][[func_type]]$results)[which(TARA_POMS_out[[env_factor]][[func_type]]$results$multinomial_corr < 0.25)]
    regress_hits <- rownames(TARA_regress_out[[env_factor]][[func_type]])[which(TARA_regress_out[[env_factor]][[func_type]]$BH < 0.25)]
    all_hits <- unique(c(POMS_hits, regress_hits))
    
    input_tree <- POMS::prep_tree(phy = TARA_tree, rownames(filt_func_table))
    pd_output <- picante::pd(samp = t(filt_func_table[, all_hits]), tree = input_tree)
    
    input_tree_edge1 <- input_tree
    input_tree_edge1$edge.length <- rep(1, length(input_tree_edge1$edge.length))
    pd_output_edge1 <- picante::pd(samp = t(filt_func_table[, all_hits]), tree = input_tree_edge1)
    
    
    TARA_output_raw[[counter]] <- data.frame(Dataset = "Tara Oceans",
                                             Env_factor = env_factor,
                                             Approach = c(rep("POMS", length(POMS_hits)),
                                                          rep("Phylo. regress. (specificity)", length(regress_hits))),
                                             Datatype = func_type,
                                             Func = c(POMS_hits, regress_hits),
                                             PD = pd_output[c(POMS_hits, regress_hits), "PD"],
                                             num_encoders = pd_output[c(POMS_hits, regress_hits), "SR"])
    
    TARA_output_edges1_raw[[counter]] <- data.frame(Dataset = "Tara Oceans",
                                                    Env_factor = env_factor,
                                                    Approach = c(rep("POMS", length(POMS_hits)),
                                                                 rep("Phylo. regress. (specificity)", length(regress_hits))),
                                                    Datatype = func_type,
                                                    Func = c(POMS_hits, regress_hits),
                                                    PD = pd_output_edge1[c(POMS_hits, regress_hits), "PD"],
                                                    num_encoders = pd_output_edge1[c(POMS_hits, regress_hits), "SR"])
    
    counter <- counter + 1
    
  }

}
  
TARA_output <- do.call(rbind, TARA_output_raw)
TARA_output_edges1 <- do.call(rbind, TARA_output_edges1_raw)

write.table(x = TARA_output,
            file = "/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits/Tara_Oceans.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(x = TARA_output_edges1,
            file = "/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits/Tara_Oceans_edges1.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


