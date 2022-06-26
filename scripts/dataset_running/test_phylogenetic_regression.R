rm(list = ls(all.names = TRUE))

library(ape)
library(phylolm)
devtools::load_all(path = "~/github_repos/POMS/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/")

# Read in input files.
almeida_ko <- read.table(file = "functional_analyses/kegg_summary.csv.gz",
                         sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)

almeida_ko <- data.frame(t(almeida_ko), check.names = FALSE)

almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")

sim_input <- readRDS("/home/gavin/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep1000.rds")

prepped_almeida_tree <- POMS::prep_tree(phy = almeida_tree, tips2keep = rownames(sim_input$taxa_perturb_abun))

almeida_ko_subset <- almeida_ko[prepped_almeida_tree$tip.label, ]

almeida_ko_subset <- almeida_ko_subset[, -which(colSums(almeida_ko_subset) == 0)]

sim_input$taxa_perturb_abun <- sim_input$taxa_perturb_abun[prepped_almeida_tree$tip.label, ]

# For this test, get significant tips based on Wilcoxon relative abundance tests,
# just to quickly get output.
wilcoxon_output <- wilcoxon_2group_pvalues(intable = sim_input$taxa_perturb_abun,
                                           group1_samples = sim_input$group1,
                                           group2_samples = sim_input$group2,
                                           convert_relab = TRUE)

sig_tips <- rownames(wilcoxon_output)[which(wilcoxon_output$wilcox_p < 0.05)]

tip_info <- data.frame(tip_ids = prepped_almeida_tree$tip.label,
                       sig = 0)

rownames(tip_info) <- tip_info$tip_ids

tip_info[which(tip_info$tip_ids %in% sig_tips), "sig"] <- 1

phylolm_sig_vs_func <- function(func_id, in_tree, in_data, in_func) {
  
  if (! identical(rownames(in_func), rownames(in_data))) {
    stop("Rownames of function and data tables must be identical.") 
  }
  
  if (! identical(in_tree$tip.label, rownames(in_data))) {
    stop("Rownames of data table and tip labels must be identical.") 
  }
  
  in_data$func_present <- 0
  
  in_data[which(in_func[, func_id] > 0), "func_present"] <- 1
  
  summary(phylolm(sig ~ func_present, data = in_data, phy = in_tree, model = "BM"))$coefficients[2, "p.value"]
  
}


func_p_raw <- mclapply(colnames(almeida_ko_subset), phylolm_sig_vs_func, in_tree = prepped_almeida_tree, in_data = tip_info, in_func = almeida_ko_subset, mc.cores = 70)
names(func_p_raw) <- colnames(almeida_ko_subset)

func_p <- unlist(func_p_raw)

func_bh <- p.adjust(func_p, "BH")

which.min(func_bh[which(func_bh < 0.05)])

func_bh_sig <- func_bh[which(func_bh < 0.05)]

func_bh_sig_rank <- rank(func_bh_sig)[sim_input$func]
