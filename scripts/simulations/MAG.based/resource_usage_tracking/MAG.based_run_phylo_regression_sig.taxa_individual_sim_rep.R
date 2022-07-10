rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages(devtools::load_all(path = "/home/gavin/github_repos/POMS/"))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(source("/home/gavin/github_repos/POMS_manuscript/scripts/alt_tool_functions.R"))

parser <- ArgumentParser()

parser$add_argument('sim_info', metavar = 'sim_info', type = "character", nargs = 1,
                    help = 'Path to simulation info RDS file for a particular replicate.')

parser$add_argument('func_table', metavar = 'func_table', type = "character", nargs = 1,
                    help = 'Path to func abundance table rds.')

parser$add_argument('tree', metavar = 'treefile', type = "character", nargs = 1,
                    help = 'Path to tree rds.')

parser$add_argument('ncores', metavar = 'num_cores', type = "integer", nargs = 1,
                    help = 'Number of cores to use.')

parser$add_argument('out_rds', metavar = 'outfile', type = "character", nargs = 1,
                    help = 'Path to output file.')

args <- parser$parse_args()



# Read in needed prepped files.
func_tab <- readRDS(args$func_table)
tree <- readRDS(args$tree)
sim_info <- readRDS(args$sim_info)

rep_wilcox_output <- wilcoxon_2group_pvalues(intable = sim_info$taxa_perturb_abun,
                                             group1_samples = sim_info$group1,
                                             group2_samples = sim_info$group2,
                                             convert_relab = TRUE)

rep_wilcox_output <- rep_wilcox_output[tree$tip.label, ]

sig_taxa <- rep(0, nrow(rep_wilcox_output))
sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1

sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                     func =  func_tab,
                                                     in_tree = tree,
                                                     ncores = args$ncores,
                                                     model_type = "BM")

sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")

saveRDS(object = sig_taxa_regress_out, file = args$out_rds)
