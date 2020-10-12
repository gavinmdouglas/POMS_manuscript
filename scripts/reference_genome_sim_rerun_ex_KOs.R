### Explore what the distribution of pos / neg enrichments looks like when the same number of balances are set to be arbitrarily significant

rm(list=ls(all.names=TRUE))

library("ape")
library("parallel")

setwd("/home/gavin/projects/POMS/reference_genome_sim/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

# Run wilcoxon tests on all of these simulated tables for comparison.
ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep="\t", header=TRUE, row.names=1)

func_sim_info_100rep <- readRDS("func_sim_info_100rep.rds")

representative_ko_POMS <- list()

in_tree <- read.tree(file = "JGI_PICRUSt2_genomes/GToTree_output/GToTree_output_subset.tre")
in_tree <- multi2di(in_tree)

ref_func <- ref_func[in_tree$tip.label, ]

ref_func <- ref_func[ , -which(colSums(ref_func) == 0)]

representative_KOs <- c("K00480", "K01845", "K06077")

for(cutoff in names(func_sim_info_100rep)) {

  print(cutoff)
  
  representative_ko_POMS[[cutoff]] <- list()
  
  BEZI_table <- read.table(file = paste("relabun_tables/BEZI_genome_abun_", cutoff, "_sigma1.tsv", sep=""),
                           header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")
  
  
  ref_func_cutoff <- ref_func[rownames(BEZI_table), ]
  
  if(length(which(colSums(ref_func_cutoff) == 0) > 0)) {
    ref_func_cutoff <- ref_func_cutoff[, -which(colSums(ref_func_cutoff) == 0)]
  }
  
  focal_genes <- sapply(func_sim_info_100rep[[cutoff]], function(x) { x$func })
  representative_KOs_i <- which(focal_genes %in% representative_KOs)
  if(length(representative_KOs_i) != 3) { stop("All three representative genes should be present.") }

  for(rep in representative_KOs_i) {
    
    print(focal_genes[rep])
    
    POMS_out <- two_group_balance_tree_pipeline(abun=func_sim_info_100rep[[cutoff]][[rep]]$taxa_perturb_abun,
                                                func=ref_func_cutoff,
                                                phylogeny=in_tree,
                                                group1_samples = func_sim_info_100rep[[cutoff]][[rep]]$group1,
                                                group2_samples = func_sim_info_100rep[[cutoff]][[rep]]$group2,
                                                ncores=20,
                                                balance_p_cutoff = 0.05,
                                                balance_correction = "BY",
                                                function_p_cutoff = 0.05,
                                                function_correction = "none",
                                                min_func_instances=0,
                                                min_func_prop=0,
                                                detailed_output = TRUE,
                                                verbose=FALSE,
                                                calc_node_dist=FALSE)
    
    
    POMS_out_null <- POM_pseudo_null(num_sig_nodes=length(POMS_out$sig_nodes),
                                        funcs_per_node=POMS_out$funcs_per_node,
                                        num_null_rep=1000,
                                        ncores=20)
    
    
    
    POMS_out_null_p <- pseudo_null_pvalues(POMS_out_null, POMS_out$df)
    
    
    func_id <- func_sim_info_100rep[[cutoff]][[rep]]$func
    
    representative_ko_POMS[[cutoff]][[func_id]] <- list(output=POMS_out, pseudo_null_p=POMS_out_null_p)
  
  }
}

saveRDS(object = representative_ko_POMS, file = "/home/gavin/github_repos/POMS_manuscript/data/ref.genome_sim_summaries/ref.genome_sim_ex_KOs.rds")
