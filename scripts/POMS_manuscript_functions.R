library(ape)
library(parallel)

threeWayVennWrapper <- function(set1, set2, set3,
                                labels=c("cat1", "cat2", "cat3"),
                                colours=c("#009E73", "#E69F00", "#56B4E9")) {
  
  set1_count <- length(set1)
  set2_count <- length(set2)
  set3_count <- length(set3)
  set1_2_count <- length(which(set1 %in% set2))
  set2_3_count <- length(which(set2 %in% set3))
  set1_3_count <- length(which(set1 %in% set3))
  set1_2_3_count_TMP <- set1[which(set1 %in% set2)]
  set1_2_3_count <- length(set1_2_3_count_TMP[which(set1_2_3_count_TMP %in% set3)])
  
  grid.newpage()
  
  venn_out <- draw.triple.venn(area1=set1_count,
                               area2=set2_count,
                               area3=set3_count,
                               n12 = set1_2_count,
                               n23 = set2_3_count,
                               n13 = set1_3_count,
                               n123 = set1_2_3_count,
                               category = labels,
                               scaled=TRUE,
                               fill = colours,
                               cex=rep(2, 7),
                               cat.cex=rep(2, 3))
  
  return(venn_out)
  
}




sim_jaccard_wrapper <- function(focal_func, sig_funcs, func_table) {
  
  focal_pos_mags <- which(func_table[, focal_func] > 0)
  
  sig_func_jaccards <- c()
  
  for (sig_func in sig_funcs) {
    
    if (sig_func == focal_func) { next }
    
    sig_func_pos_mags <- which(func_table[, sig_func] > 0)
    num_intersecting_mags <- length(which(focal_pos_mags %in% sig_func_pos_mags))
    num_focal_only_mags <- length(which(!focal_pos_mags %in% sig_func_pos_mags))
    num_sig_only_mags <- length(which(!sig_func_pos_mags %in% focal_pos_mags))
    
    sig_func_jaccard <- num_intersecting_mags / (num_intersecting_mags + num_focal_only_mags + num_sig_only_mags)
    sig_func_jaccards <- c(sig_func_jaccards, sig_func_jaccard)
  }
  
  if (length(sig_func_jaccards) > 0) {
    return(mean(sig_func_jaccards))
  } else {
    return(NA)
  }
}


simulation_summaries <- function(POMS_sims,func_table, alt_tool_sims=NULL, focal_func_present = TRUE,
                                 skip_alt_tools=FALSE, num_cores=1, sig_cutoffs=c(0.25, 0.05, 0.001)) {
  
  num_tested_functions <- ncol(func_table)
  
  if (length(alt_tool_sims) != length(POMS_sims)) {
    stop("Different number of sim reps between POMS and alt tools lists.") 
  }
  
  if ( (!skip_alt_tools) && (is.null(alt_tool_sims)) ) {
    stop("Option skip_alt_tools is FALSE by no alt_tool_sums object specified.") 
  }
  
  sim_summaries <- data.frame(matrix(NA, ncol = 0, nrow = length(POMS_sims)))
  
  if (focal_func_present) {
    sim_summaries[, "focal_func"] <- sapply(POMS_sims, function(x) { return(x$func) })
    sim_summaries[, "num_focal_pos_mags"] <- sapply(POMS_sims, function(x) { return(length(which(func_table[, x$func] > 0)))})
  }
  
  sim_summaries[, "POMS_num_sig_nodes"] <- sapply(POMS_sims, function(x) { return(length(x$output$sig_nodes)) })
  
  for (sig_cutoff in sig_cutoffs) {
    
    if (focal_func_present) {
      POMS_raw_summary <- as.data.frame(do.call(rbind, mclapply(POMS_sims, POMS_sim_rep_summary, sig_cutoff = sig_cutoff, func_table = func_table, mc.cores = num_cores)))
      colnames(POMS_raw_summary) <- paste(colnames(POMS_raw_summary), as.character(sig_cutoff), sep = "_")
      sim_summaries <- cbind(sim_summaries, POMS_raw_summary)
    } else {
      POMS_col <- paste("POMS_sig", as.character(sig_cutoff), sep = "_")
      sim_summaries[, POMS_col] <- sapply(POMS_sims, function(x) { return(length(which(x$output$df$multinomial_corr < sig_cutoff)) / num_tested_functions) })
    }
    
    for (alt_tool in names(alt_tool_sims[[1]])) {

      if (alt_tool == "func") { next }
      
      if (focal_func_present) {
        
        rep_summary_raw <- as.data.frame(do.call(rbind, mclapply(1:length(alt_tool_sims), 
                                                         function(i) {
                                                           
                                                           return(alt_sim_rep_summary(alt.tool_sim_rep = alt_tool_sims[[i]][[alt_tool]],
                                                                                      focal_func = sim_summaries[i, "focal_func"],
                                                                                      func_table = func_table,
                                                                                      sig_cutoff = sig_cutoff,
                                                                                      tool_name = alt_tool))
                                                           
                                                         }, mc.cores = num_cores)))
        
        colnames(rep_summary_raw) <- paste(colnames(rep_summary_raw), as.character(sig_cutoff), sep = "_")
        sim_summaries <- cbind(sim_summaries, rep_summary_raw)
        
      } else {
        alt_tool_col <- paste(alt_tool, as.character(sig_cutoff), sep = "_")
        sim_summaries[, alt_tool_col] <- sapply(alt_tool_sims, function(x) { return(length(which(x[[alt_tool]]$BH_corr_p < sig_cutoff)) / num_tested_functions) })
      }
      
    }
    
  }
  
  return(sim_summaries)
  
}


POMS_sim_rep_summary <- function(POMS_sim_rep, sig_cutoff, func_table) {
  
  focal_func <- POMS_sim_rep$func
  sig_funcs <- rownames(POMS_sim_rep$output$df)[which(POMS_sim_rep$output$df$multinomial_corr < sig_cutoff)]
  sig_funcs_P <- POMS_sim_rep$output$df[sig_funcs, "multinomial_corr"]
  
  prop_sig <- length(sig_funcs) / ncol(func_table)
  
  if (focal_func %in% sig_funcs) {
    focal_rank <- rank(sig_funcs_P)[which(sig_funcs == focal_func)]
  } else {
    focal_rank <- NA  
  }
  
  focal_rel_rank <- focal_rank / length(sig_funcs)
  focal_jaccard <- sim_jaccard_wrapper(focal_func = focal_func, sig_funcs = sig_funcs, func_table = func_table)
  
  stats <- c(prop_sig, focal_rank, focal_rel_rank, focal_jaccard)
  names(stats) <- c("POMS_sig", "POMS_rank", "POMS_rel_rank", "POMS_jaccard")
  
  return(stats)
  
}


alt_sim_rep_summary <- function(alt.tool_sim_rep, focal_func, func_table, sig_cutoff, tool_name) {
  
  corr_P <- alt.tool_sim_rep$BH_corr_p
  
  sig_funcs <- rownames(alt.tool_sim_rep)[which(corr_P < sig_cutoff)]
  sig_funcs_P <- corr_P[which(corr_P < sig_cutoff)]
  
  prop_sig <- length(sig_funcs) / ncol(func_table)
  
  if (focal_func %in% sig_funcs) {
    focal_rank <- rank(sig_funcs_P)[which(sig_funcs == focal_func)]
  } else {
    focal_rank <- NA  
  }
  
  focal_rel_rank <- focal_rank / length(sig_funcs)
  focal_jaccard <- sim_jaccard_wrapper(focal_func = focal_func, sig_funcs = sig_funcs, func_table = func_table)
  
  stats <- c(prop_sig, focal_rank, focal_rel_rank, focal_jaccard)
  names(stats) <- c("sig", "rank", "rel_rank", "jaccard")
  names(stats) <- paste(tool_name, names(stats), sep = "_")
  
  return(stats)
  
}


# Modified version of pipeline that can work with samples with continuous data rather than discrete groups.
POMS_pipeline_continuous <- function(abun,
                                     func,
                                     phylogeny,
                                     group1_samples,
                                     group2_samples,
                                     ncores=1,
                                     pseudocount=1,
                                     significant_nodes=NULL,
                                     tested_balances=NULL,
                                     nodes_dir=NULL,
                                     min_num_tips=10,
                                     min_func_instances=10,
                                     min_func_prop=0.001,
                                     multinomial_min_sig=5,
                                     balance_p_cutoff = 0.05,
                                     balance_correction = "none",
                                     function_p_cutoff = 0.05,
                                     function_correction = "none",
                                     func_descrip_infile = NULL,
                                     run_multinomial_test=TRUE,
                                     multinomial_correction="BH",
                                     calc_node_dist=FALSE,
                                     detailed_output=FALSE,
                                     verbose=FALSE) {
  
  if(verbose) { message("Prepping input phylogeny.") }
  phylogeny <- prep_tree(phy=phylogeny, tips2keep=rownames(abun))
  
  if((min_func_instances > 0) || (min_func_prop > 0)) {
    func <- filter_rare_table_cols(func,  min_func_instances, min_func_prop, verbose)
  }
  
  if(is.null(significant_nodes)) {
    
    if(verbose) { message("Calculating balances.") }
    
    calculated_balances <- compute_tree_node_balances(abun=abun, phylogeny=phylogeny, ncores=ncores, min_num_tips = min_num_tips)
    
    if(verbose) { message("Identified ", length(calculated_balances$balances), " of ", length(phylogeny$node.label), " nodes as non-negligible and will be used for analyses.") }
    
    if(verbose) { message("Running Wilcoxon tests to test for differences in balances between groups at each non-negligible node.") }
    
    pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances, group1_samples, group2_samples, corr_method=balance_correction, skip_wilcoxon=FALSE)
    
    sig_nodes <- names(calculated_balances$balances)[which(pairwise_node_out$wilcox_corrected_p < balance_p_cutoff)]
    
  } else {
    
    sig_nodes <- significant_nodes
    
    if(any(! sig_nodes %in% phylogeny$node.label)) { stop("Not all sig. nodes are not found in phylogeny.")}
    
    calculated_balances <- list()
    calculated_balances$balances <- tested_balances
    calculated_balances$features <- parallel::mclapply(names(tested_balances),
                                                       lhs_rhs_asvs,
                                                       tree=phylogeny,
                                                       get_node_index=TRUE,
                                                       mc.cores=ncores)
    
    names(calculated_balances$features) <- names(tested_balances)
    
    negligible_nodes_i <- which(! names(phylogeny$node.label) %in% names(tested_balances))
    if(length(negligible_nodes_i) > 0) {
      calculated_balances$negligible_nodes <- phylogeny$node.label[negligible_nodes_i]
    } else {
      calculated_balances$negligible_nodes <- c()
    }
    
    if(is.null(nodes_dir)) {
      pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances, group1_samples, group2_samples, skip_wilcoxon=TRUE)
    } else {
      pairwise_node_out <- list(mean_direction = nodes_dir)
    }
  }
  
  if(verbose) { message("Identifying enriched functions at all non-negligible nodes.") }
  
  all_balances_enriched_funcs <- mclapply(names(calculated_balances$balances),
                                          function(x) {
                                            return(node_func_fisher(node = x,
                                                                    in_tree = phylogeny,
                                                                    in_func = func,
                                                                    higher_group=pairwise_node_out$mean_direction[x],
                                                                    pseudocount=1,
                                                                    multiple_test_corr=function_correction))
                                          },
                                          mc.cores = ncores)
  
  names(all_balances_enriched_funcs) <- names(calculated_balances$balances)
  
  if(verbose) { message("Summarizing significant functions across nodes.") }
  
  func_summaries <- summarize_node_enrichment(all_balances_enriched_funcs, sig_nodes, function_p_cutoff)
  
  # Get single DF summarizing the key metrics and print this out.
  all_func_id <- c()
  for(balance in names(all_balances_enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(all_balances_enriched_funcs[[balance]]))
  }
  
  all_func_id <- all_func_id[-which(duplicated(all_func_id))]
  
  if(verbose) { message("Creating results dataframe.") }
  
  summary_df <- data.frame(matrix(NA, nrow=length(all_func_id), ncol=4))
  
  rownames(summary_df) <- all_func_id
  
  colnames(summary_df) <- c("num_nodes_enriched",
                            "num_sig_nodes_group1_enrich",
                            "num_sig_nodes_group2_enrich",
                            "num_nonsig_nodes_enrich")
  
  rownames(summary_df) <- all_func_id
  
  if(verbose) { message("Creating results dataframe.") }
  
  
  if(calc_node_dist) {
    
    if(verbose) { message("Calculating inter-node distance") }    
    
    phylogeny_node_dists <- dist.nodes(phylogeny)
    
    summary_df$mean_internode_dist_neg_enrich <- NA
    summary_df$max_internode_dist_neg_enrich <- NA
    summary_df$mean_internode_dist_pos_enrich <- NA
    summary_df$max_internode_dist_pos_enrich <- NA
    
  }
  
  if(run_multinomial_test) {
    
    if(verbose) { message("Will run multinomial test on every function (that meet the multinomial_min_sig cut-off).") } 
    
    prop_sig_node_balances <- length(sig_nodes) / length(calculated_balances$balances)
    
    multinomial_exp_prop <- c(prop_sig_node_balances * 0.5, prop_sig_node_balances * 0.5, 1 - prop_sig_node_balances)
    
    names(multinomial_exp_prop) <- c("exp_sig_nodes_group1_enrich_prop", "exp_sig_nodes_group2_enrich_prop", "exp_nonsig_nodes_enrich_prop")
    
    summary_df$multinomial_p <- NA
  }
  
  
  for(func_id in all_func_id) {
    
    summary_df[func_id, c("num_sig_nodes_group1_enrich",
                          "num_sig_nodes_group2_enrich",
                          "num_nonsig_nodes_enrich")] <- c(length(func_summaries[[func_id]]$positive_nodes),
                                                           length(func_summaries[[func_id]]$negative_nodes),
                                                           length(func_summaries[[func_id]]$enriched_nonsig_nodes))
    
    
    summary_df[func_id, "num_nodes_enriched"] <- sum(as.numeric(summary_df[func_id, c("num_sig_nodes_group1_enrich",
                                                                                      "num_sig_nodes_group2_enrich",
                                                                                      "num_nonsig_nodes_enrich")]))
    
    all_nodes_present <- c(func_summaries[[func_id]]$nonenriched_sig_nodes,
                           func_summaries[[func_id]]$positive_nodes,
                           func_summaries[[func_id]]$negative_nodes,
                           func_summaries[[func_id]]$nonenriched_nonsig_nodes,
                           func_summaries[[func_id]]$enriched_nonsig_nodes)
    
    if(max(table(all_nodes_present)) > 1) {
      stop("Node categorized into at least 2 mutually exclusive groups.")
    }
    
    if(run_multinomial_test) {
      
      observed_counts <- as.numeric(summary_df[func_id, c("num_sig_nodes_group1_enrich",
                                                          "num_sig_nodes_group2_enrich",
                                                          "num_nonsig_nodes_enrich")])
      
      if((length(sig_nodes) > 0) && (summary_df[func_id, "num_nodes_enriched"] >= multinomial_min_sig) && (prop_sig_node_balances != 1)) {
        summary_df[func_id, "multinomial_p"] <- XNomial::xmulti(obs=observed_counts,
                                                                expr=multinomial_exp_prop, detail=0)$pProb 
      }
      
    }
    
    if(calc_node_dist) {
      
      summary_df[func_id, c("mean_internode_dist_group1_enrich",
                            "max_internode_dist_group1_enrich",
                            "mean_internode_dist_group2_enrich",
                            "max_internode_dist_group2_enrich")] <- c(internode_mean_max_dist(phy = phylogeny, dist_matrix = phylogeny_node_dists,
                                                                                              node_labels = func_summaries[[func_id]]$positive_nodes),
                                                                      internode_mean_max_dist(phy = phylogeny, dist_matrix = phylogeny_node_dists,
                                                                                              node_labels = func_summaries[[func_id]]$negative_nodes))
    }
    
  }
  
  if((run_multinomial_test) && (multinomial_correction != "none")) {
    summary_df$multinomial_corr <- p.adjust(summary_df$multinomial_p, multinomial_correction)
  }
  
  if(! is.null(func_descrip_infile)) {
    if(verbose) { message("Adding function descriptions to output.") }
    func_descrip <- read.table(func_descrip_infile,
                               header=FALSE, sep="\t", row.names=1, stringsAsFactors = FALSE, quote="")
    summary_df$description <- func_descrip[rownames(summary_df), "V2"]
  } else {
    if(verbose) { message("Function description mapfile not specified (func_descrip_infile argument), so no descriptions will be added.") }    
  }
  
  results <- list(balances_info=calculated_balances,
                  sig_nodes=sig_nodes,
                  df=summary_df)
  
  if(run_multinomial_test) {
    results[["multinomial_exp_prop"]] <- multinomial_exp_prop
  }
  
  if(detailed_output) {
    results[["balance_comparisons"]] <- pairwise_node_out
    results[["funcs_per_node"]] <- all_balances_enriched_funcs
    results[["out_list"]] <- func_summaries
    results[["tree"]] <- phylogeny
  }
  
  return(results)
  
}


node_balances_spearman_cor <- function(node_balances, sample_info, sample_var) {
  
  if (length(which(is.na(sample_info[, sample_var]))) > 0) {
    sample_info <- sample_info[-which(is.na(sample_info[, sample_var])), , drop = FALSE]
  }
  
  node_spearman_cor <- data.frame(matrix(NA, nrow = length(node_balances$balances), ncol = 2))
  colnames(node_spearman_cor) <- c("rho", "p")
  rownames(node_spearman_cor) <- names(node_balances$balances)
  
  for (node in names(node_balances$balances)) {
    
    cor_out <- cor.test(node_balances$balances[[node]][rownames(sample_info)],
                        sample_info[ , sample_var],
                        method = "spearman", exact = FALSE)
    
    node_spearman_cor[node, ] <- c(cor_out$estimate, cor_out$p.value) 
  }
  
  return(node_spearman_cor)
}