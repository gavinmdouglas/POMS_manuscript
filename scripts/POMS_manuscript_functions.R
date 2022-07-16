library(ape)
library(parallel)

calc_func_abun_crossprod <- function(in_abun, in_func) {

  # Check that all rows are found in function table.
  if(length(which(! rownames(in_abun) %in% rownames(in_func))) > 0) {
    stop("Stoppings - some rows in abundance table not found in function table.")
  }
  
  in_func <- in_func[rownames(in_abun), ]
  
  return(data.frame(t(crossprod(as.matrix(in_abun), as.matrix(in_func)))))
  
}

subset_abun_table <- function(in_abun, col2keep) {

  in_abun <- in_abun[, which(colnames(in_abun) %in% col2keep)]

  missing_rows <- which(rowSums(in_abun) == 0)
  missing_samples <- which(colSums(in_abun) == 0)

  if(length(missing_rows) > 0) {
    in_abun <- in_abun[-missing_rows, ]
  }

  if(length(missing_samples) > 0) {
    in_abun <- in_abun[, -missing_samples]
  }

  return(in_abun)
}

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



simulation_summaries <- function(POMS_sims, func_table, regress_sims = NULL, alt_tool_sims=NULL, focal_func_present = TRUE,
                                 skip_alt_tools=FALSE, num_cores=1, sig_cutoffs=c(0.25, 0.05, 0.001)) {
  
  num_tested_functions <- ncol(func_table)
  
  if ((! skip_alt_tools) && (length(alt_tool_sims) != length(POMS_sims))) {
    stop("Different number of sim reps between POMS and alt tools output lists.")
  }
  
  if (length(regress_sims) != length(POMS_sims)) {
    stop("Different number of sim reps between POMS and regress output lists.")
  }
  
  if ( (!skip_alt_tools) && (is.null(alt_tool_sims)) ) {
    stop("Option skip_alt_tools is FALSE by no alt_tool_sums object specified.") 
  }
  
  sim_summaries <- data.frame(matrix(NA, ncol = 0, nrow = length(POMS_sims)))
  
  if (focal_func_present) {
    sim_summaries[, "focal_func"] <- sapply(POMS_sims, function(x) { return(x$func) })
    sim_summaries[, "num_focal_pos_mags"] <- sapply(POMS_sims, function(x) { return(length(which(func_table[, x$func] > 0)))})
  }
  
  sim_summaries[, "POMS_num_BSNs"] <- sapply(POMS_sims, function(x) { return(length(x$output$BSNs)) })
  
  for (sig_cutoff in sig_cutoffs) {
    
    if (focal_func_present) {
      POMS_raw_summary <- as.data.frame(do.call(rbind, mclapply(POMS_sims, POMS_sim_rep_summary, sig_cutoff = sig_cutoff, func_table = func_table, mc.cores = num_cores)))
      colnames(POMS_raw_summary) <- paste(colnames(POMS_raw_summary), as.character(sig_cutoff), sep = "_")
      sim_summaries <- cbind(sim_summaries, POMS_raw_summary)
      
      regress_raw_summary <- as.data.frame(do.call(rbind, mclapply(regress_sims, regress_sim_rep_summary, sig_cutoff = sig_cutoff, func_table = func_table, mc.cores = num_cores)))
      colnames(regress_raw_summary) <- paste(colnames(regress_raw_summary), as.character(sig_cutoff), sep = "_")
      sim_summaries <- cbind(sim_summaries, regress_raw_summary)
      
    } else {
      
      POMS_col <- paste("POMS_sig", as.character(sig_cutoff), sep = "_")
      sim_summaries[, POMS_col] <- sapply(POMS_sims, function(x) { return(length(which(x$output$results$multinomial_corr < sig_cutoff)) / num_tested_functions) })
      
      regress_specificity_col <- paste("regress_specificity_sig", as.character(sig_cutoff), sep = "_")
      sim_summaries[, regress_specificity_col] <- sapply(regress_sims, function(x) {
                                          if ((length(x$specificity_regress) > 1) || (x$specificity_regress[1] != "Failed")) {
                                            return(length(which(x$specificity_regress$BH < sig_cutoff)) / num_tested_functions)
                                          } else {
                                            return(NA)
                                          }
                                          
        })
      
      regress_sig.taxa_col <- paste("regress_sig.taxa_sig", as.character(sig_cutoff), sep = "_")
      sim_summaries[, regress_sig.taxa_col] <- sapply(regress_sims, function(x) { return(length(which(x$sig_taxa_regress$BH < sig_cutoff)) / num_tested_functions) })
      
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


regress_sim_rep_summary <- function(regress_sim_rep, sig_cutoff, func_table) {
  
  focal_func <- regress_sim_rep$func
  
  if (is.null(focal_func)) { stop("Focal func not found") }
  
  # specificity
  if ((length(regress_sim_rep$specificity_regress) > 1) || (regress_sim_rep$specificity_regress[1] != "Failed")) {
    sig_specificity_funcs <- rownames(regress_sim_rep$specificity_regress)[which(regress_sim_rep$specificity_regress$BH < sig_cutoff)]
    sig_specificity_funcs_P <- regress_sim_rep$specificity_regress[sig_specificity_funcs, "BH"]
    
    specificity_func_prop_sig <- length(sig_specificity_funcs) / ncol(func_table)

    if (focal_func %in% sig_specificity_funcs) {
      specificity_focal_rank <- rank(sig_specificity_funcs_P)[which(sig_specificity_funcs == focal_func)]
      specificity_focal_rel_rank <- specificity_focal_rank / length(sig_specificity_funcs)
    } else {
      specificity_focal_rank <- NA
      specificity_focal_rel_rank <- NA
    }

    specificity_focal_jaccard <- sim_jaccard_wrapper(focal_func = focal_func, sig_funcs = sig_specificity_funcs, func_table = func_table)
  
  } else {

    specificity_func_prop_sig <- NA
    specificity_focal_rank <- NA
    specificity_focal_rel_rank <- NA
    specificity_focal_jaccard <- NA

  }
  
  # sig_taxa
  sig_sig_taxa_funcs <- rownames(regress_sim_rep$sig_taxa_regress)[which(regress_sim_rep$sig_taxa_regress$BH < sig_cutoff)]
  sig_sig_taxa_funcs_P <- regress_sim_rep$sig_taxa_regress[sig_sig_taxa_funcs, "BH"]
  
  sig_taxa_func_prop_sig <- length(sig_sig_taxa_funcs) / ncol(func_table)
  
  if (focal_func %in% sig_sig_taxa_funcs) {
    sig_taxa_focal_rank <- rank(sig_sig_taxa_funcs_P)[which(sig_sig_taxa_funcs == focal_func)]
    sig_taxa_focal_rel_rank <- sig_taxa_focal_rank / length(sig_sig_taxa_funcs)
  } else {
    sig_taxa_focal_rank <- NA 
    sig_taxa_focal_rel_rank <- NA
  }
  
  sig_taxa_focal_jaccard <- sim_jaccard_wrapper(focal_func = focal_func, sig_funcs = sig_sig_taxa_funcs, func_table = func_table)
  
  stats <- c(specificity_func_prop_sig, specificity_focal_rank, specificity_focal_rel_rank, specificity_focal_jaccard,
             sig_taxa_func_prop_sig, sig_taxa_focal_rank, sig_taxa_focal_rel_rank, sig_taxa_focal_jaccard)
  
  names(stats) <- c("regress_specificity_func_prop_sig", "regress_specificity_rank",
                    "regress_specificity_rel_rank", "regress_specificity_jaccard",
                    "regress_sig_taxa_func_prop_sig", "regress_sig_taxa_rank",
                    "regress_sig_taxa_rel_rank", "regress_sig_taxa_jaccard")
  
  return(stats)
  
}

POMS_sim_rep_summary <- function(POMS_sim_rep, sig_cutoff, func_table) {
  
  focal_func <- POMS_sim_rep$func
  sig_funcs <- rownames(POMS_sim_rep$output$results)[which(POMS_sim_rep$output$results$multinomial_corr < sig_cutoff)]
  sig_funcs_P <- POMS_sim_rep$output$results[sig_funcs, "multinomial_corr"]
  
  prop_sig <- length(sig_funcs) / ncol(func_table)
  
  if (focal_func %in% sig_funcs) {
    focal_rank <- rank(sig_funcs_P)[which(sig_funcs == focal_func)]
    focal_rel_rank <- focal_rank / length(sig_funcs)
  } else {
    focal_rank <- NA
    focal_rel_rank <- NA
  }
  
  
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

