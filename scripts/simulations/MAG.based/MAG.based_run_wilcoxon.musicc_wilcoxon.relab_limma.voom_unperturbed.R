rm(list = ls(all.names = TRUE))

# Run wilcoxon.musicc, wilcoxon.relab, and limma.voom for 1000 replicates (these are run together because they run very quickly)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1

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
                                 tools_to_run = c("limma.voom", "wilcoxon.musicc", "wilcoxon.relab"),
                                 USCGs = musicc_uscgs)

  alt_tools_out[["func"]] <- "unperturbed"

  return(alt_tools_out)

}, mc.cores = 30)


saveRDS(object = unperturbed_alt.tools, file = "MAG.based_wilcoxon.musicc_wilcoxon.relab_limma.voom_unperturbed_1000reps.rds")

print(proc.time() - ptm)
