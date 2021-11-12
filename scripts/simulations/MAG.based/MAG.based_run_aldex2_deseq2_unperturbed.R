rm(list = ls(all.names = TRUE))

# Run deseq2 and aldex2 on unperturbed profiles
 
setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

ptm <- proc.time()

almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")
random_groups <- readRDS("almeida_healthy_random_groups.rds")

almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- subset_abun_table(in_abun = almeida_abun,
                                  col2keep = c(random_groups$group1[[1]], random_groups$group2[[1]]))

almeida_func_abun <- calc_func_abun(in_abun = almeida_abun, in_func = almeida_func_subset, ncores = 30)

unperturbed_alt.tools <- mclapply(1:1000, function(rep_i) {
  
  alt_tools_out <- run_alt.tools(func_abun_table = almeida_func_abun,
                                 group1_samples = random_groups$group1[[rep_i]],
                                 group2_samples = random_groups$group2[[rep_i]],
                                 tools_to_run = c("aldex2", "deseq2"))

  alt_tools_out[["func"]] <- "unperturbed"

  out_rds <- paste("MAG.based_unperturbed_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = "")
  
  saveRDS(object = alt_tools_out, file = out_rds)
  
  return("success")

}, mc.cores = 10)

print(proc.time() - ptm)

