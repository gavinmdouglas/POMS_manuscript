
rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/reference_genome_sim/")

source("/home/gavin/github_repos/POMS/balance_tree_functions.R")

library(parallel)

# Read in pre-determined random sample groupings.
random_groups <- readRDS("ref_genomes_random_groups.rds")

BEZI_tables <- list()

BEZI_tables[["mu0.1_nu_0.5"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.5_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.9"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.9_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.99_sigma1.tsv",
                                             header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.01_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.01_nu_0.99_sigma1.tsv",
                                             header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

in_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                      sep="\t", header=TRUE, row.names=1)

in_tree <- read.tree(file = "JGI_PICRUSt2_genomes/GToTree_output/GToTree_output_subset.tre")



cutoff_func <- list()
cutoff_tree <- list()

for(cutoff in names(BEZI_tables)) {
  
  cutoff_func[[cutoff]] <- in_func[rownames(BEZI_tables[[cutoff]]), ]
  
  if(length(which(colSums(cutoff_func[[cutoff]]) == 0)) > 0) {
    cutoff_func[[cutoff]] <- cutoff_func[[cutoff]][, -which(colSums(cutoff_func[[cutoff]]) == 0)]
  }

  tips2drop <- in_tree$tip.label[which(! in_tree$tip.label %in% rownames(BEZI_tables[[cutoff]]))]

  if(length(tips2drop) > 0) {
    cutoff_tree[[cutoff]] <- multi2di(drop.tip(phy = in_tree, tip = tips2drop, trim.internal = TRUE))
  } else {
    cutoff_tree[[cutoff]] <- multi2di(in_tree)
  }
  
}


strictest_cutoff <- "mu0.01_nu_0.99"
set.seed(141515)
in_func_no_rare <- filter_rare_table_cols(in_tab = cutoff_func[[strictest_cutoff]], min_nonzero_count = 5, min_nonzero_prop = 0.001, verbose=TRUE)
random_funcs <- sample(colnames(in_func_no_rare), size = 1000)
write.table(x =random_funcs, file="1000_random_KOs.txt", col.names = FALSE, row.names=FALSE, quote=FALSE)


taxa_sim_info <- list()
func_sim_info <- list()

for(cutoff in names(BEZI_tables)) {
  
  taxa_sim_info[[cutoff]] <- mclapply(1:1000,
  
                                                            function(rep) {
  
                                                              group1_samples <- random_groups$group1[[rep]]
                                                              group2_samples <- random_groups$group2[[rep]]
  
                                                              func <- random_funcs[rep]
  
                                                              orig_contributors <- rownames(cutoff_func[[cutoff]])[which(cutoff_func[[cutoff]][, func] > 0)]
  
                                                              ran_contributors <- sample(rownames(cutoff_func[[cutoff]]), length(orig_contributors))
  
                                                              tmp_abun <- BEZI_tables[[cutoff]]
  
                                                              missing_group1 <- which(! group1_samples %in% colnames(tmp_abun))
                                                              missing_group2 <- which(! group2_samples %in% colnames(tmp_abun))
  
                                                              if(length(missing_group1) > 0) { group1_samples <- group1_samples[-missing_group1] }
                                                              if(length(missing_group2) > 0) { group2_samples <- group2_samples[-missing_group2] }
  
                                                              tmp_abun[ran_contributors, group1_samples] <- (tmp_abun[ran_contributors, group1_samples] + 1) * 10
  
                                                              return(list(taxa_perturb_abun=tmp_abun,
                                                                          orig_contributors=orig_contributors,
                                                                          ran_contributors=ran_contributors,
                                                                          func=func,
                                                                          group1=group1_samples,
                                                                          group2=group2_samples,
                                                                          selected_group="group1"))
                                                            }, mc.cores = 30)


    func_sim_info[[cutoff]] <- mclapply(1:1000,
                                 
                                                   function(rep) {
                                                     
                                                     group1_samples <- random_groups$group1[[rep]]
                                                     group2_samples <- random_groups$group2[[rep]]
                                                     
                                                     func <- random_funcs[rep]
                                                     
                                                     orig_contributors <- rownames(cutoff_func[[cutoff]])[which(cutoff_func[[cutoff]][, func] > 0)]
                                                     
                                                     tmp_abun <- BEZI_tables[[cutoff]]
                                                     
                                                     missing_group1 <- which(! group1_samples %in% colnames(tmp_abun))
                                                     missing_group2 <- which(! group2_samples %in% colnames(tmp_abun))
                                                     
                                                     if(length(missing_group1) > 0) { group1_samples <- group1_samples[-missing_group1] }
                                                     if(length(missing_group2) > 0) { group2_samples <- group2_samples[-missing_group2] }
                                                     
                                                     tmp_abun[orig_contributors, group1_samples] <- (tmp_abun[orig_contributors, group1_samples] + 1) * 10
                                                     
                                                     return(list(taxa_perturb_abun=tmp_abun,
                                                                 orig_contributors=orig_contributors,
                                                                 func=func,
                                                                 group1=group1_samples,
                                                                 group2=group2_samples,
                                                                 selected_group="group1"))
                                                   }, mc.cores = 30)

}

saveRDS(object = taxa_sim_info, file = "taxa_sim_info.rds")
saveRDS(object = func_sim_info, file = "func_sim_info.rds")

taxa_sim_info_100rep <- taxa_sim_info
func_sim_info_100rep <- func_sim_info

for(cutoff in names(taxa_sim_info_100rep)) {
  taxa_sim_info_100rep[[cutoff]] <- taxa_sim_info_100rep[[cutoff]][1:100]
  func_sim_info_100rep[[cutoff]] <- func_sim_info_100rep[[cutoff]][1:100]
}

saveRDS(object = taxa_sim_info_100rep, file = "taxa_sim_info_100rep.rds")
saveRDS(object = func_sim_info_100rep, file = "func_sim_info_100rep.rds")

