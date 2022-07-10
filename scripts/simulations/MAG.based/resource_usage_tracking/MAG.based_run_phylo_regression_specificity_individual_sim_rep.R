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

rep_metadata <- data.frame(samp = c(sim_info$group1,
                                    sim_info$group2),
                           group = c(rep("group1", length(sim_info$group1)),
                                     rep("group2", length(sim_info$group2))))

taxa_specificity <- specificity_scores(abun_table = sim_info$taxa_perturb_abun,
                                       meta_table = rep_metadata,
                                       focal_var_level = "group1",
                                       var_colname = "group",
                                       sample_colname = "samp",
                                       silence_citation = TRUE)

taxa_specificity <- taxa_specificity$ess[tree$tip.label]

specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                        func =  func_tab,
                                                        in_tree = tree,
                                                        ncores = 1,
                                                        model_type = "BM")

htop$BH <- p.adjust(specificity_regress_out$p, "BH")

saveRDS(object = specificity_regress_out, file = args$out_rds)
