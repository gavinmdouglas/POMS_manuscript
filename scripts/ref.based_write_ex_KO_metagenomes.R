### Compare results identified when jittering an input table when run through
### basic Wilcoxon test on functions and with POMS.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/reference_genome_sim/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")
devtools::load_all(path = "/home/gavin/github_repos/POMS/")

representative_KOs <- c("K00480", "K01845", "K06077")

# Run wilcoxon tests on all of these simulated tables for comparison.
ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep="\t", header=TRUE, row.names=1)

func_sim_info_100rep <- readRDS("func_sim_info_100rep.rds")


BEZI_table <- list()
ref_func_subset <- list()
metagenomes <- list()

for(cutoff in names(func_sim_info_100rep)) {
  
  BEZI_table[[cutoff]] <- read.table(file = paste("relabun_tables/BEZI_genome_abun_", cutoff, "_sigma1.tsv", sep=""),
                                      header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")
  
  
  ref_func_subset[[cutoff]] <- ref_func[rownames( BEZI_table[[cutoff]]), ]
  
  if(length(which(colSums(ref_func_subset[[cutoff]]) == 0) > 0)) {
    ref_func_subset[[cutoff]] <- ref_func_subset[[cutoff]][, -which(colSums(ref_func_subset[[cutoff]]) == 0)]
  }
  
  metagenomes[[cutoff]] <- list()
  
  cutoff_func_ids <- sapply(func_sim_info_100rep[[cutoff]], function(x) { x$func })
  
  representative_KOs_i <- which(cutoff_func_ids %in% representative_KOs)
  
  for(i in representative_KOs_i) {
    func_id <- cutoff_func_ids[i]
    metagenomes[[cutoff]][[func_id]] <- calc_func_abun(in_abun = func_sim_info_100rep[[cutoff]][[i]]$taxa_perturb_abun, in_func = ref_func_subset[[cutoff]], ncores = 60)
  }
  
}

saveRDS(object=metagenomes, file="/home/gavin/github_repos/POMS_manuscript/data/ref.genome_simulated_abundances/ex_KO_metagenome_tables.rds")
