library(ALDEx2)
library(ape)
library(DESeq2)
library(edgeR)
library(limma)
library(parallel)

# Run other differential abundance tools.
run_2group_ALDEx2 <- function(in_table, group1_samples, group2_samples, divide_sum=1) {
  in_table <- floor(in_table[, c(group1_samples, group2_samples)] / divide_sum)
  return(aldex(reads = in_table, conditions=c(rep("group1", length(group1_samples)), rep("group2", length(group2_samples)))))
}




wilcoxon_2group_pvalues <- function(intable, group1_samples, group2_samples, convert_relab=FALSE) {
  
  group1_intable <- intable[, group1_samples]
  group2_intable <- intable[, group2_samples]
  
  if(convert_relab) {
    group1_intable <- data.frame(t(sweep(x = group1_intable, MARGIN = 2, STATS = colSums(group1_intable), FUN = '/')), check.names=FALSE)
    group2_intable <- data.frame(t(sweep(x = group2_intable, MARGIN = 2, STATS = colSums(group2_intable), FUN = '/')), check.names=FALSE)
  } else {
    group1_intable <- data.frame(t(group1_intable))
    group2_intable <- data.frame(t(group2_intable))
  }
  
  wilcox_out_df <- data.frame(matrix(NA, nrow=ncol(group1_intable), ncol=4))
  rownames(wilcox_out_df) <- colnames(group1_intable)
  colnames(wilcox_out_df) <- c("mean_group1", "mean_group2", "wilcox_W", "wilcox_p")
  
  for(feature in colnames(group1_intable)) {
    
    raw_wilcoxon_out <- wilcox.test(group1_intable[, feature], group2_intable[, feature])
    wilcox_out_df[feature, c("wilcox_W", "wilcox_p")] <- c(raw_wilcoxon_out$statistic, raw_wilcoxon_out$p.value)
    
    wilcox_out_df[feature, "mean_group1"] <- mean(as.numeric(group1_intable[, feature]))
    wilcox_out_df[feature, "mean_group2"] <- mean(as.numeric(group2_intable[, feature]))
  }
  
  return(wilcox_out_df)
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

calc_func_abun <- function(in_abun, in_func, ncores=1) {
  
  out_df <- data.frame(matrix(NA, nrow=ncol(in_func), ncol=ncol(in_abun)))
  colnames(out_df) <- colnames(in_abun)
  rownames(out_df) <- colnames(in_func)
  
  # Check that all rows are found in function table.
  if(length(which(! rownames(in_abun) %in% rownames(in_func))) > 0) {
    stop("Stoppings - some rows in abundance table not found in function table.")
  }
  
  in_func <- in_func[rownames(in_abun), ]
  
  out_sample_func_abun <- mclapply(colnames(in_abun), function(x) { return(colSums(in_abun[, x] * in_func)) }, mc.cores=ncores)
  names(out_sample_func_abun) <- colnames(in_abun)
  
  for(sample in colnames(in_abun)) {
    out_df[, sample] <- out_sample_func_abun[[sample]]
  }
  
  return(out_df)
  
}


deseq2_default_two_groups <- function(table, group1, group2, alpha.set=0.05, divide_sum=1) {
  
  group1_subset <- group1[which(group1 %in% colnames(table))]
  group2_subset <- group2[which(group2 %in% colnames(table))]
  
  # Create metadata df.
  metadata_tab <- data.frame(matrix(NA, nrow=(length(group1_subset) + length(group2_subset)), ncol=1))
  rownames(metadata_tab) <- c(group1_subset, group2_subset)
  colnames(metadata_tab) <- c("group")
  metadata_tab[group1_subset, "group"] <- "group1"
  metadata_tab[group2_subset, "group"] <- "group2"
  metadata_tab$group <- as.factor(metadata_tab$group)
  
  # Input count df needs to have columns in same order as metadata rows.
  table_subset <- floor(table[, rownames(metadata_tab)] / divide_sum)
  
  dds <- DESeqDataSetFromMatrix(countData = table_subset,
                                colData = metadata_tab,
                                design = ~ group)
  default_deseq2 <- DESeq(dds)
  default_deseq2_results <- results(default_deseq2, alpha=alpha.set)
  
  return(default_deseq2_results)
}


limma_voom_two_group_TMM <- function(table, group1, group2) {
  
  group1_subset <- group1[which(group1 %in% colnames(table))]
  group2_subset <- group2[which(group2 %in% colnames(table))]
  
  # Create metadata df.
  metadata_tab <- data.frame(matrix(NA, nrow=(length(group1_subset) + length(group2_subset)), ncol=1))
  rownames(metadata_tab) <- c(group1_subset, group2_subset)
  colnames(metadata_tab) <- c("group")
  metadata_tab[group1_subset, "group"] <- "group1"
  metadata_tab[group2_subset, "group"] <- "group2"
  metadata_tab$group <- as.factor(metadata_tab$group)
  
  counts <- floor(as.matrix(table))
  
  dge <- DGEList(counts = counts)
  
  ### Check if upper quartile method works for selecting reference
  upper_quartile_norm_test <- calcNormFactors(dge, method="upperquartile")
  
  summary_upper_quartile <- summary(upper_quartile_norm_test$samples$norm.factors)[3]
  if(is.na(summary_upper_quartile)){
    message("Upper Quartile reference selection failed will use find sample with largest sqrt(read_depth) to use as reference")
    Ref_col <- which.max(colSums(sqrt(counts)))
    dge_norm <- calcNormFactors(dge, method = "TMM", refColumn = Ref_col)
  }else{
    dge_norm <- calcNormFactors(dge, method="TMM")
  }
  
  model_matrix <- model.matrix(as.formula("~ group"), metadata_tab)
  
  voom_out <- voom(dge_norm, model_matrix, plot=FALSE)
  
  voom_out_fit <- lmFit(voom_out, model_matrix)
  voom_out_fit_eBayes <- eBayes(voom_out_fit)
  
  return(topTable(voom_out_fit_eBayes, coef = 2, n = nrow(dge_norm), sort.by="none"))
  
}

run_alt.tools <- function(func_abun_table, group1_samples, group2_samples, USCGs=NULL,
                          tools_to_run=c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc")) {
  
  if (length(which(!tools_to_run %in% c("aldex2", "deseq2", "limma.voom", "wilcoxon.relab", "wilcoxon.musicc"))) > 0) {
    stop("At least one tool specified is not in the expected set.")
  }
  
  if (length(tools_to_run) == 0) {
    stop("At least one tool to run needs to be given.")
  }
  
  if (("wilcoxon.musicc" %in% tools_to_run) && (is.null(USCGs))) {
    stop("Must set USCGs argument to run wilcoxon.musicc")
  }
  
  if (("aldex2" %in% tools_to_run) || ("deseq2" %in% tools_to_run)) {
    func_abun_table_ceil <- ceiling(func_abun_table)
  }
  
  DA_out <- list()
  
  if ("aldex2" %in% tools_to_run) {
    DA_out[["aldex2"]] <- run_2group_ALDEx2(in_table = func_abun_table_ceil,
                                            group1_samples = group1_samples,
                                            group2_samples = group2_samples)
    
    DA_out[["aldex2"]]$BH_corr_p <- DA_out[["aldex2"]]$wi.eBH
  }
  
  if ("deseq2" %in% tools_to_run) {
    DA_out[["deseq2"]] <- deseq2_default_two_groups(table = func_abun_table_ceil,
                                                    group1 = group1_samples,
                                                    group2 = group2_samples)
    
    DA_out[["deseq2"]]$BH_corr_p <- p.adjust(DA_out[["deseq2"]]$pvalue, "BH")
  }
  
  if ("limma.voom" %in% tools_to_run) {
    DA_out[["limma.voom"]] <- limma_voom_two_group_TMM(table = func_abun_table,
                                                       group1 = group1_samples,
                                                       group2 = group2_samples)
    
    DA_out[["limma.voom"]]$BH_corr_p <- p.adjust(DA_out[["limma.voom"]]$P.Value, "BH")
  }
  
  if ("wilcoxon.relab" %in% tools_to_run) {
    DA_out[["wilcoxon.relab"]] <- wilcoxon_2group_pvalues(intable = func_abun_table,
                                                          group1_samples = group1_samples,
                                                          group2_samples = group2_samples,
                                                          convert_relab = TRUE)
    
    DA_out[["wilcoxon.relab"]]$BH_corr_p <- p.adjust(DA_out[["wilcoxon.relab"]]$wilcox_p, "BH")
  }
  
  if ("wilcoxon.musicc" %in% tools_to_run) {
    
    dataset_uscg_set <- USCGs[which(USCGs %in% rownames(func_abun_table))]
    func_abun_table_musicc <- data.frame(sweep(func_abun_table, 2, colMedians(as.matrix(func_abun_table[dataset_uscg_set, ])), `/`))
    
    DA_out[["wilcoxon.musicc"]] <- wilcoxon_2group_pvalues(intable = func_abun_table_musicc,
                                                           group1_samples = group1_samples,
                                                           group2_samples = group2_samples,
                                                           convert_relab = FALSE)
    
    DA_out[["wilcoxon.musicc"]]$BH_corr_p <- p.adjust(DA_out[["wilcoxon.musicc"]]$wilcox_p, "BH")
  }
  
  return(DA_out)
  
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
    
    if (!skip_alt_tools) {
      
      if ("aldex2" %in% names(alt_tool_sims[[1]])) { stop("Need to add code aldex2") }
      if ("deseq2" %in% names(alt_tool_sims[[1]])) { stop("Need to add code deseq2") }
      if ("limma.voom" %in% names(alt_tool_sims[[1]])) { stop("Need to add code limma.voom") }
      
      
      if ("wilcoxon.relab.out" %in% names(alt_tool_sims[[1]])) {
        
        if (focal_func_present) {
          wilcoxon.relab_raw <- as.data.frame(do.call(rbind, mclapply(1:length(alt_tool_sims), 
                                                                      function(i) {
                                                                        
                                                                        return(Wilcoxon_sim_rep_summary(alt.tool_sim_rep = alt_tool_sims[[i]]$wilcoxon.relab.out,
                                                                                                        focal_func = sim_summaries[i, "focal_func"],
                                                                                                        func_table = func_table,
                                                                                                        sig_cutoff = sig_cutoff,
                                                                                                        tool_name = "wilcoxon.relab"))
                                                                        
                                                                      }, mc.cores = num_cores)))
          
          colnames(wilcoxon.relab_raw) <- paste(colnames(wilcoxon.relab_raw), as.character(sig_cutoff), sep = "_")
          sim_summaries <- cbind(sim_summaries, wilcoxon.relab_raw)
          
        } else {
          wilcoxon.relab_col <- paste("wilcoxon.relab", as.character(sig_cutoff), sep = "_")
          sim_summaries[, wilcoxon.relab_col] <- sapply(alt_tool_sims, function(x) { return(length(which(p.adjust(x$wilcoxon.relab.out$wilcox_p, "BH") < sig_cutoff)) / num_tested_functions) })
        }
      }
      
      if ("wilcoxon.musicc.out" %in% names(alt_tool_sims[[1]])) {
        
        if (focal_func_present) {
          wilcoxon.musicc_raw <- as.data.frame(do.call(rbind, mclapply(1:length(alt_tool_sims), 
                                                                       function(i) {
                                                                         
                                                                         return(Wilcoxon_sim_rep_summary(alt.tool_sim_rep = alt_tool_sims[[i]]$wilcoxon.musicc.out,
                                                                                                         focal_func = sim_summaries[i, "focal_func"],
                                                                                                         func_table = func_table,
                                                                                                         sig_cutoff = sig_cutoff,
                                                                                                         tool_name = "wilcoxon.musicc"))
                                                                         
                                                                       }, mc.cores = num_cores)))
          
          colnames(wilcoxon.musicc_raw) <- paste(colnames(wilcoxon.musicc_raw), as.character(sig_cutoff), sep = "_")
          sim_summaries <- cbind(sim_summaries, wilcoxon.musicc_raw)
        } else {
          wilcoxon.musicc_col <- paste("wilcoxon.musicc", as.character(sig_cutoff), sep = "_")
          sim_summaries[, wilcoxon.musicc_col] <- sapply(alt_tool_sims, function(x) { return(length(which(p.adjust(x$wilcoxon.musicc.out$wilcox_p, "BH") < sig_cutoff)) / num_tested_functions) })
        }
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


Wilcoxon_sim_rep_summary <- function(alt.tool_sim_rep, focal_func, func_table, sig_cutoff, tool_name) {
  
  corr_P <- p.adjust(alt.tool_sim_rep$wilcox_p, "BH")
  
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

