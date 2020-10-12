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




POMS_jaccard_wrapper <- function(POMS_output, func_table, num_sig_nodes, num_other_nodes, num_rep) {
  
  POMS_focal_func_set <- c()
  POMS_mean_jaccard <- c()
  
  for(rep_index in c(1:num_rep)) {
    
    focal_func <- POMS_output[[rep_index]]$func
    
    focal_pos_mags <- which(func_table[, focal_func] > 0)
    
    pos_indices <- which(POMS_output[[rep_index]]$output$df$num_sig_nodes_pos_enrich >= num_sig_nodes & POMS_output[[rep_index]]$output$df$num_sig_nodes_neg_enrich <= num_other_nodes)
    neg_indices <- which(POMS_output[[rep_index]]$output$df$num_sig_nodes_neg_enrich >= num_sig_nodes & POMS_output[[rep_index]]$output$df$num_sig_nodes_pos_enrich <= num_other_nodes)
    
    pos_sig_func <- POMS_output[[rep_index]]$output$df[pos_indices, "func"]
    neg_sig_func <- POMS_output[[rep_index]]$output$df[neg_indices, "func"]
    
    if(focal_func %in% pos_sig_func) {
      POMS_focal_func_set <- c(POMS_focal_func_set, "pos")
    } else if(focal_func %in% neg_sig_func) {
      POMS_focal_func_set <- c(POMS_focal_func_set, "neg")
    } else {
      POMS_focal_func_set <- c(POMS_focal_func_set, "non_sig")
    }
    
    all_sig_func <- c(pos_sig_func, neg_sig_func)
    
    sig_func_jaccards <- c()
    for(sig_func in all_sig_func) {
      
      if(sig_func == focal_func) { next }
      
      sig_func_pos_mags <- which(func_table[, sig_func] > 0)
      num_intersecting_mags <- length(which(focal_pos_mags %in% sig_func_pos_mags))
      num_focal_only_mags <- length(which(! focal_pos_mags %in% sig_func_pos_mags))
      num_sig_only_mags <- length(which(! sig_func_pos_mags %in% focal_pos_mags))
      
      sig_func_jaccard <- num_intersecting_mags / (num_intersecting_mags + num_focal_only_mags + num_sig_only_mags)
      sig_func_jaccards <- c(sig_func_jaccards, sig_func_jaccard)
    }
    
    if(length(sig_func_jaccards) > 0) {
      POMS_mean_jaccard <- c(POMS_mean_jaccard, mean(sig_func_jaccards))
    } else {
      POMS_mean_jaccard <- c(POMS_mean_jaccard, NA)
    }
  }
  
  return(POMS_mean_jaccard)
}


wilcoxon_jaccard_wrapper <- function(POMS_output, wilcoxon_output, func_table, p_cutoff, num_rep) {
  
  wilcoxon_mean_jaccard <- c()
  
  for(rep_index in c(1:num_rep)) {
    
    f <- POMS_output[[rep_index]]$func
    
    f_pos_mags <- which(func_table[, f] > 0)
    
    wilcoxon_func <- colnames(func_table)[which(wilcoxon_output[[rep_index]] < p_cutoff)]
    
    wilcoxon_sig_func_jaccards <- c()
    
    for(wilcoxon_sig_func in wilcoxon_func) {
      if(wilcoxon_sig_func == f) { next }
      
      wilcoxon_sig_func_pos_mags <- which(func_table[, wilcoxon_sig_func] > 0)
      
      num_wilcoxon_intersecting_mags <- length(which(f_pos_mags %in% wilcoxon_sig_func_pos_mags))
      num_wilcoxon_f_only_mags <- length(which(! f_pos_mags %in% wilcoxon_sig_func_pos_mags))
      num_wilcoxon_sig_only_mags <- length(which(! wilcoxon_sig_func_pos_mags %in% f_pos_mags))
      
      wilcoxon_sig_func_jaccard <- num_wilcoxon_intersecting_mags / (num_wilcoxon_intersecting_mags + num_wilcoxon_f_only_mags + num_wilcoxon_sig_only_mags)
      wilcoxon_sig_func_jaccards <- c(wilcoxon_sig_func_jaccards, wilcoxon_sig_func_jaccard)
    }
    
    if(length(wilcoxon_sig_func_jaccards) > 0) {
      wilcoxon_mean_jaccard <- c(wilcoxon_mean_jaccard, mean(wilcoxon_sig_func_jaccards))
    } else {
      wilcoxon_mean_jaccard <- c(wilcoxon_mean_jaccard, NA)
    }
  }
  
  return(wilcoxon_mean_jaccard)
}

POMS_wilcoxon_rank_wrapper <- function(POMS_output, wilcoxon_output, func_table, num_rep) {
  
  func_ids <- c()
  all_num_pos_mags <- c()
  POMS_rel_ranks <- c()
  wilcoxon_rel_ranks <- c()
  
  for(rep in c(1:num_rep)) {
    
    f <- POMS_output[[rep]]$func
    
    func_ids <- c(func_ids, f)
    
    f_i <- which(colnames(func_table) == f)
    
    num_pos_mags <- length(which(func_table[, f] > 0))
    all_num_pos_mags <- c(all_num_pos_mags, num_pos_mags)
    
    wilcoxon_rank <- rank(wilcoxon_output[[rep]])[f_i]
    wilcoxon_rel_ranks <- c(wilcoxon_rel_ranks, wilcoxon_rank / length(wilcoxon_output[[rep]]))
    
    pos_num_sig_node_diff <- abs(POMS_output[[rep]]$output$df$num_sig_nodes_pos_enrich - POMS_output[[rep]]$output$df$num_sig_nodes_neg_enrich)
    names(pos_num_sig_node_diff) <- rownames(POMS_output[[rep]]$output$df)
    
    if(f %in% names(pos_num_sig_node_diff)) {
      f_POMS_i <- which(names(pos_num_sig_node_diff) == f)
      POMS_focal_rank <- rank(-pos_num_sig_node_diff)[f_POMS_i]
      POMS_rel_ranks <- c(POMS_rel_ranks, POMS_focal_rank / length(pos_num_sig_node_diff))
    } else {
      POMS_rel_ranks <- c(POMS_rel_ranks, NA)
    }
  }
  
  ranking_df <- data.frame(replicate=1:num_rep,
                           func=func_ids,
                           wilcoxon_rel_ranks=wilcoxon_rel_ranks,
                           POMS_rel_ranks=POMS_rel_ranks,
                           all_num_pos_mags=all_num_pos_mags)
  
  return(ranking_df)
  
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


deseq2_default_two_groups <- function(table, group1, group2, alpha.set=0.05, dataset_name="unknown", divide_sum=1) {
  
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
