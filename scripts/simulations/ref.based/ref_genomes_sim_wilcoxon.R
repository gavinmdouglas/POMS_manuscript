### Compare results identified when jittering an input table when run through
### basic Wilcoxon test on functions and with POMS.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/reference_genome_sim/")

source("/home/gavin/github_repos/POMS/balance_tree_functions.R")

# Run wilcoxon tests on all of these simulated tables for comparison.
ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                      sep="\t", header=TRUE, row.names=1)

func_sim_info_100rep <- readRDS("func_sim_info_100rep.rds")


func_rand_wilcoxon <- list()



for(cutoff in names(func_sim_info_100rep)) {
  
  BEZI_table <- read.table(file = paste("relabun_tables/BEZI_genome_abun_", cutoff, "_sigma1.tsv", sep=""),
                           header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")
  
  
  ref_func_cutoff <- ref_func[rownames(BEZI_table), ]
  
  if(length(which(colSums(ref_func_cutoff) == 0) > 0)) {
    ref_func_cutoff <- ref_func_cutoff[, -which(colSums(ref_func_cutoff) == 0)]
  }
  
  func_rand_wilcoxon[[cutoff]] <- mclapply(X = 1:100, FUN = function(x) {
    tmp_metagenome <- calc_func_abun(in_abun = func_sim_info_100rep[[cutoff]][[x]]$taxa_perturb_abun, in_func = ref_func_cutoff, ncores = 1)
    return(wilcoxon_2group_pvalues(tmp_metagenome, func_sim_info_100rep[[cutoff]][[x]]$group1, func_sim_info_100rep[[cutoff]][[x]]$group2))
  }
  , mc.cores=30)
  
}

saveRDS(object = func_rand_wilcoxon, file="func_rand_output_wilcoxon_100reps_BEZI_cutoffs.rds")
