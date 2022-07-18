rm(list = ls(all.names = TRUE))

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/")

TARA_POMS_out <- readRDS(file = "TARA_POMS_out.rds")

considered_variables <- c("Mean_Temperature", "Mean_Salinity", "Mean_Oxygen", "Mean_Nitrates", "NO2", "PO4")

combined_TARA_sig_df <- TARA_POMS_out$Chlorophyll.Sensor.s$pathway$results[-c(1:nrow(TARA_POMS_out$Chlorophyll.Sensor.s$pathway$results)), ]

for (sample_var in considered_variables) {
  
  for (func_type in c("ko", "pathway", "module")) {
  
    POMS_out <- TARA_POMS_out[[sample_var]][[func_type]]

    POMS_out_sig <- POMS_out$results[which(POMS_out$results$multinomial_corr < 0.25), ]

    if (nrow(POMS_out_sig) > 0) {
      POMS_out_sig$func <- rownames(POMS_out_sig)
      rownames(POMS_out_sig) <- NULL
      POMS_out_sig$variable <- sample_var
      POMS_out_sig$func_type <- func_type
      combined_TARA_sig_df <- rbind(combined_TARA_sig_df, POMS_out_sig)
    }
  }
}


write.table(x = combined_TARA_sig_df,
            file = "/home/gavin/github_repos/POMS_manuscript/display_items/Maintext_TARA_sig_func_RAW.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)




# Also get a similar table summarizing the phylogenetic regression results, so that they can be more sensibly compared.
TARA_regress_out <- readRDS(file = "TARA_regress_out.rds")

for (var in considered_variables) {
  for (func_level in names(TARA_regress_out[[var]])) {
    TARA_regress_out[[var]][[func_level]]$func <- rownames(TARA_regress_out[[var]][[func_level]])
    TARA_regress_out[[var]][[func_level]]$func_type <- func_level
    TARA_regress_out[[var]][[func_level]]$variable <- var
  }
}

TARA_regress_out_tmp <- list()

for (var in considered_variables) {
  TARA_regress_out_tmp[[var]] <- do.call(rbind, TARA_regress_out[[var]])
}

TARA_regress_combined <- do.call(rbind, TARA_regress_out_tmp)

TARA_regress_combined_sig <- TARA_regress_combined[which(TARA_regress_combined$BH < 0.25), ]


write.table(x = TARA_regress_combined_sig,
            file = "/home/gavin/github_repos/POMS_manuscript/key_datafiles/TARA_regress_results.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)

# Quick comparison
combined_TARA_sig_df$func_variable <- paste(combined_TARA_sig_df$func, "-", combined_TARA_sig_df$variable)
TARA_regress_combined_sig$func_variable <- paste(TARA_regress_combined_sig$func, "-", TARA_regress_combined_sig$variable)

combined_TARA_sig_df$func_variable[which(combined_TARA_sig_df$func_variable %in% TARA_regress_combined_sig$func_variable)]