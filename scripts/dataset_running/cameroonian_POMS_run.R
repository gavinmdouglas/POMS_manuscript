### Explore what the distribution of pos / neg enrichments looks like when the same number of balances are set to be arbitrarily significant

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/MAGs/Almeida2019/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

cameroon_abun <- read.table(file = "/home/gavin/projects/POMS/MIDAS_running/cameroon_merged_species/count_reads.txt", header=TRUE, sep="\t", row.names=1)
cameroon_abun <- cameroon_abun[-which(rowSums(cameroon_abun) == 0), ]

cameroon_sample_info <- read.table("/home/gavin/projects/POMS/MIDAS_running/cameroon_MIDAS_metadata.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")
cameroon_group1 <- cameroon_sample_info$sample[which(cameroon_sample_info$env == "Yes")]
cameroon_group2 <- cameroon_sample_info$sample[which(cameroon_sample_info$env == "No")]

in_func <- read.table(file = "/home/gavin/projects/POMS/MIDAS_running/midas_db_v1.2_FIGFAM.tsv",
                      sep="\t", header=TRUE, check.names = FALSE, row.names=1)
in_func <- in_func[, rownames(cameroon_abun), ]
# Remove function not in at least 0.1% of genomes.
func_nonzero_percent <- (rowSums(in_func > 0) / ncol(in_func)) * 100
in_func <- in_func[-which(func_nonzero_percent < 0.1), ]
in_func <- as.data.frame(t(in_func))


# Change tree tip labels to be full names
midas_tree <- read.tree(file = "/home/gavin/local/database/midas_db_v1.2/species_tree.newick")
strain_full_ids <- rownames(in_func)
names(strain_full_ids) <- gsub("^.*_", "", strain_full_ids)
midas_tree <- drop.tip(midas_tree, tip=which(! midas_tree$tip.label %in% names(strain_full_ids)))
midas_tree$tip.label <- as.character(strain_full_ids[midas_tree$tip.label])

cameroon_POMS_out <- two_group_balance_tree_pipeline(abun=cameroon_abun,
                                                         func=in_func,
                                                         phylogeny=midas_tree,
                                                         group1_samples = cameroon_group1,
                                                         group2_samples = cameroon_group2,
                                                         ncores=50,
                                                         calc_node_dist=FALSE,
                                                         detailed_output = TRUE,
                                                         verbose=TRUE,
                                                          #NOTE NO MULTIPLE-TEST CORRECTION PERFORMED ON BALANCES:
                                                         balance_correction = "none")

cameroon_POMS_out_null <- POM_pseudo_null(num_sig_nodes=length(cameroon_POMS_out$sig_nodes),
                                         funcs_per_node=cameroon_POMS_out$funcs_per_node,
                                         num_null_rep=1000,
                                         ncores=60)


cameroon_POMS_out_null_p <- pseudo_null_pvalues(cameroon_POMS_out_null, cameroon_POMS_out$df)


# Only sig figfam based on this approach: FIG00003555

superoxide <- read.table("/home/gavin/tmp/superoxide_figam.txt", header=FALSE, sep="\t")

cameroon_POMS_out$df[as.character(superoxide$V1), ]


saveRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/POMS_vs_phylogenize/cameroonian_POMS_out.rds",
        object = list(df=cameroon_POMS_out, null=cameroon_POMS_out_null, null_sig=which(p.adjust(cameroon_POMS_out_null_p, "BY") < 0.05)))
