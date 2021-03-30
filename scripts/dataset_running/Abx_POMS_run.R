rm(list=ls(all.names=TRUE))

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

Abx_abun <- read.table(file = "/home/gavin/projects/POMS/MIDAS_running/Abx_merged_species/count_reads_collapsed.txt", header=TRUE, sep="\t", row.names=1)

Abx_sample_info <- read.table("/home/gavin/projects/POMS/MIDAS_running/Abx_MIDAS_metadata.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")
Abx_group1 <- Abx_sample_info$sample[which(Abx_sample_info$env == "seven_days")]
Abx_group2 <- Abx_sample_info$sample[which(Abx_sample_info$env == "baseline")]

in_func <- read.table(file = "/home/gavin/projects/POMS/MIDAS_running/midas_db_v1.2_FIGFAM.tsv",
                      sep="\t", header=TRUE, check.names = FALSE, row.names=1)
in_func <- in_func[, rownames(Abx_abun), ]
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

Abx_POMS_out <- two_group_balance_tree_pipeline(abun=Abx_abun,
                                                         func=in_func,
                                                         phylogeny=midas_tree,
                                                         group1_samples = Abx_group1,
                                                         group2_samples = Abx_group2,
                                                         ncores=60,
                                                         calc_node_dist=FALSE,
                                                         detailed_output = TRUE,
                                                         verbose=TRUE,
                                                         balance_correction = "none")

Abx_POMS_out_null <- POM_pseudo_null(num_sig_nodes=length(Abx_POMS_out$sig_nodes),
                                         funcs_per_node=Abx_POMS_out$funcs_per_node,
                                         num_null_rep=1000,
                                         ncores=60)


Abx_POMS_out_null_p <- pseudo_null_pvalues(Abx_POMS_out_null, Abx_POMS_out$df)


figfam_descrip <- read.table("/home/gavin/local/database/midas_db_v1.2/ontologies/figfam.txt",
                             sep="\t", header=FALSE, stringsAsFactors = FALSE, quote="")

rownames(figfam_descrip) <- figfam_descrip$V1

all_sig <- names(Abx_POMS_out_null_p)[which(Abx_POMS_out_null_p < 0.05)]

all_sig_df <- Abx_POMS_out$df[all_sig, ]


up_sig <- rownames(all_sig_df)[which(all_sig_df$num_sig_nodes_pos_enrich > all_sig_df$num_sig_nodes_neg_enrich)]
down_sig <- rownames(all_sig_df)[which(all_sig_df$num_sig_nodes_pos_enrich < all_sig_df$num_sig_nodes_neg_enrich)]


figfam_descrip[up_sig, ]
figfam_descrip[down_sig, ]


saveRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/POMS_vs_phylogenize/Abx_POMS_out.rds",
        object = list(df=Abx_POMS_out, null=Abx_POMS_out_null, null_sig=which(p.adjust(Abx_POMS_out_null_p, "BY") < 0.05)))
