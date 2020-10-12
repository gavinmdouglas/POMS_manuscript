rm(list=ls(all.names=TRUE))

library(parallel)

setwd("/home/gavin/projects/POMS/reference_genome_sim/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")
devtools::load_all(path = "/home/gavin/github_repos/POMS/")

representative_KOs <- c("K00480", "K01845", "K06077")

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/MUSiCC_KEGG_single_copy_genes.txt", stringsAsFactors = FALSE)$V1

# Run wilcoxon tests on all of these simulated tables for comparison.
ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep="\t", header=TRUE, row.names=1)

func_sim_info_100rep <- readRDS("func_sim_info_100rep.rds")

metagenomes <- readRDS(file="/home/gavin/github_repos/POMS_manuscript/data/ref.genome_simulated_abundances/ex_KO_metagenome_tables.rds")


BEZI_table <- list()
ref_func_subset <- list()
groupings <- list()

for(cutoff in names(func_sim_info_100rep)) {
  
  BEZI_table[[cutoff]] <- read.table(file = paste("relabun_tables/BEZI_genome_abun_", cutoff, "_sigma1.tsv", sep=""),
                                      header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")
  
  
  ref_func_subset[[cutoff]] <- ref_func[rownames( BEZI_table[[cutoff]]), ]
  
  if(length(which(colSums(ref_func_subset[[cutoff]]) == 0) > 0)) {
    ref_func_subset[[cutoff]] <- ref_func_subset[[cutoff]][, -which(colSums(ref_func_subset[[cutoff]]) == 0)]
  }
  
  groupings[[cutoff]] <- list()
  
  cutoff_func_ids <- sapply(func_sim_info_100rep[[cutoff]], function(x) { x$func })
  
  representative_KOs_i <- which(cutoff_func_ids %in% representative_KOs)
  
  for(i in representative_KOs_i) {
    func_id <- cutoff_func_ids[i]
    groupings[[cutoff]][[func_id]] <- list(group1=func_sim_info_100rep[[cutoff]][[i]]$group1, group2=func_sim_info_100rep[[cutoff]][[i]]$group2)
  }

}


return_filt_tab <- function(x, in_list) {
  
  orig_tab <- in_list[[x]]
  rounded_tab <- floor(orig_tab)
  
  if(length(which(colSums(rounded_tab) == 0) > 0)) {
    orig_tab <- orig_tab[, -which(colSums(orig_tab) == 0)]
  }
  
  if(length(which(rowSums(rounded_tab) == 0) > 0)) {
    orig_tab <- orig_tab[-which(rowSums(orig_tab) == 0), ]
  }
  
  return(orig_tab)
}


metagenomes_filt <- list()
for(cutoff in names(metagenomes)) {
  metagenomes_filt[[cutoff]] <- lapply(representative_KOs, return_filt_tab, in_list=metagenomes[[cutoff]])
  names(metagenomes_filt[[cutoff]]) <- representative_KOs
}

#alt_tools <- c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc")

# Leave out deseq2 for now (prob. wont be in manuscript anyway)
alt_tools <- c("aldex2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc")

categories_to_test <- list()

for(cutoff in names(metagenomes)) {
  for(tool in alt_tools) {
    for(gene in representative_KOs) {
      category <- paste(cutoff, tool, gene, sep=" ")
      categories_to_test[[category]] <- list()
      categories_to_test[[category]]$cutoff <- cutoff
      categories_to_test[[category]]$tool <- tool
      categories_to_test[[category]]$gene <- gene
    }
  }
}
      
      
ex_KOs_DA_out <- mclapply(names(categories_to_test),
                          function(x) {
                            
                            cutoff <- categories_to_test[[x]]$cutoff
                            tool <- categories_to_test[[x]]$tool
                            gene <- categories_to_test[[x]]$gene
                            
                            group1_samples_subset <- groupings[[cutoff]][[gene]]$group1[which(groupings[[cutoff]][[gene]]$group1 %in% colnames(metagenomes_filt[[cutoff]][[gene]]))]
                            group2_samples_subset <- groupings[[cutoff]][[gene]]$group2[which(groupings[[cutoff]][[gene]]$group2 %in% colnames(metagenomes_filt[[cutoff]][[gene]]))]
                            
                            if(tool == "aldex2") { 
                                return(run_2group_ALDEx2(in_table = metagenomes_filt[[cutoff]][[gene]],
                                                         group1_samples = group1_samples_subset,
                                                         group2_samples = group2_samples_subset,
                                                         divide_sum=1e6))
                            
                            } else if(tool == "deseq2") {
                                return(deseq2_default_two_groups(table = metagenomes_filt[[cutoff]][[gene]],
                                                                 group1 = group1_samples_subset,
                                                                 group2 = group2_samples_subset,
                                                                 divide_sum=1e6))
                            
                            } else if(tool == "limma.voom") {
                                return(limma_voom_two_group_TMM(table = metagenomes_filt[[cutoff]][[gene]],
                                                                group1 = group1_samples_subset,
                                                                group2 = group2_samples_subset))
                            
                            } else if(tool == "wilcoxon.relab") {
                                return(wilcoxon_2group_pvalues(intable = metagenomes_filt[[cutoff]][[gene]],
                                                               group1_samples = group1_samples_subset,
                                                               group2_samples = group2_samples_subset,
                                                               convert_relab = TRUE))
                            
                            } else if(tool == "wilcoxon.musicc") {
                                uscg_set <- musicc_uscgs[which(musicc_uscgs %in% rownames(metagenomes_filt[[cutoff]][[gene]]))]
                                tmp_metagenome_musicc <- data.frame(sweep(metagenomes_filt[[cutoff]][[gene]], 2, colMedians(as.matrix(metagenomes_filt[[cutoff]][[gene]][uscg_set, ])), `/`))
                                return(wilcoxon_2group_pvalues(intable = tmp_metagenome_musicc,
                                                               group1_samples = group1_samples_subset,
                                                               group2_samples = group2_samples_subset,
                                                               convert_relab = FALSE))
                                
                            } else { stop(paste("tool not matched:", x)) }
                            
                          
                          }, mc.cores = 30)

      

# Write out all DA output results

names(ex_KOs_DA_out) <- categories_to_test

saveRDS(object = ex_KOs_DA_out, file = "/home/gavin/github_repos/POMS_manuscript/data/ref.based_ex_KOs_DA_out.rds")

