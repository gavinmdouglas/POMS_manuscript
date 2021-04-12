library(ALDEx2)
library(DESeq2)
library(edgeR)
library(limma)

# Functions for running alternative differential abundance tools.
# These are packaged separately because if ALDEx2 is run with parallel loaded then it will overload the system.

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


deseq2_default_two_groups <- function(intable, group1, group2, alpha.set=0.05, divide_sum=1) {
  
  group1_subset <- group1[which(group1 %in% colnames(intable))]
  group2_subset <- group2[which(group2 %in% colnames(intable))]
  
  # Create metadata df.
  metadata_tab <- data.frame(matrix(NA, nrow = (length(group1_subset) + length(group2_subset)), ncol = 1))
  rownames(metadata_tab) <- c(group1_subset, group2_subset)
  colnames(metadata_tab) <- c("group")
  metadata_tab[group1_subset, "group"] <- "group1"
  metadata_tab[group2_subset, "group"] <- "group2"
  metadata_tab$group <- as.factor(metadata_tab$group)
  
  # Input count df needs to have columns in same order as metadata rows.
  table_subset <- floor(intable[, rownames(metadata_tab)] / divide_sum)
  
  dds <- DESeqDataSetFromMatrix(countData = table_subset,
                                colData = metadata_tab,
                                design = ~ group)
  default_deseq2 <- DESeq(dds, parallel = FALSE)
  default_deseq2_results <- results(default_deseq2, alpha=alpha.set, parallel = FALSE)
  
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
    DA_out[["deseq2"]] <- deseq2_default_two_groups(intable = func_abun_table_ceil,
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