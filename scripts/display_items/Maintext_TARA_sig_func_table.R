rm(list = ls(all.names = TRUE))

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/")

TARA_POMS_out <- readRDS(file = "TARA_POMS_out.rds")

combined_TARA_sig_df <- TARA_POMS_out$Chlorophyll.Sensor.s$pathway$df[-c(1:nrow(TARA_POMS_out$Chlorophyll.Sensor.s$pathway$df)), ]

for (sample_var in names(TARA_POMS_out)) {
  
  for (func_type in c("ko", "pathway", "module")) {
  
    POMS_out <- TARA_POMS_out[[sample_var]][[func_type]]

    POMS_out_sig <- POMS_out$df[which(POMS_out$df$multinomial_corr < 0.25), ]

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
