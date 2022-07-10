rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages(devtools::load_all(path = "/home/gavin/github_repos/POMS/"))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(parallel))


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

output <- POMS_pipeline(abun = sim_info$taxa_perturb_abun,
                        func = func_tab,
                        tree = tree,
                        group1_samples = sim_info$group1,
                        group2_samples = sim_info$group2,
                        ncores = args$ncores,
                        BSN_p_cutoff = 0.05,
                        BSN_correction = "none",
                        FSN_p_cutoff = 0.05,
                        FSN_correction = "none",
                        min_func_instances = 0,
                        min_func_prop = 0,
                        multinomial_correction = "BH",
                        detailed_output = FALSE,
                        verbose = FALSE)

saveRDS(object = output, file = args$out_rds)
