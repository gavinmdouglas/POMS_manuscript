rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep = "\t", header = TRUE, row.names = 1)

func_sim_info <- readRDS("func_sim_info.rds")

func_rand_POMS <- list()

in_tree <- read.tree(file = "GToTree_output_subset.tre")

ref_func <- ref_func[in_tree$tip.label, ]

ref_func <- ref_func[ , -which(colSums(ref_func) == 0)]


for (cutoff in names(func_sim_info)) {
  
  ptm <- proc.time()
  
  BEZI_table <- read.table(file = paste("relabun_tables/BEZI_genome_abun_", cutoff, "_sigma1.tsv.gz", sep = ""),
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")
  
  
  ref_func_cutoff <- ref_func[rownames(BEZI_table), ]
  
  if (length(which(colSums(ref_func_cutoff) == 0) > 0)) {
    ref_func_cutoff <- ref_func_cutoff[, -which(colSums(ref_func_cutoff) == 0)]
  }
  
  func_rand_POMS[[cutoff]] <- mclapply(X = 1:1000, FUN = function(rep_i) {

    output <- POMS_pipeline(abun = func_sim_info[[cutoff]][[rep_i]]$taxa_perturb_abun,
                            func = ref_func_cutoff,
                            phylogeny = in_tree,
                            group1_samples = func_sim_info[[cutoff]][[rep_i]]$group1,
                            group2_samples = func_sim_info[[cutoff]][[rep_i]]$group2,
                            ncores = 1,
                            balance_p_cutoff = 0.05,
                            balance_correction = "none",
                            function_p_cutoff = 0.05,
                            function_correction = "none",
                            min_func_instances = 0,
                            min_func_prop = 0,
                            run_multinomial_test = TRUE,
                            multinomial_correction = "BH",
                            calc_node_dist = FALSE,
                            detailed_output = FALSE,
                            verbose = FALSE)
    
    return(list(func = func_sim_info[[cutoff]][[rep_i]]$func, output = output))
  }
  , mc.cores = 20)
  
  
  print(proc.time() - ptm)
  
}

saveRDS(object = func_rand_POMS, file="ref.based_POMS_sim_rand_func_1000reps.rds")
