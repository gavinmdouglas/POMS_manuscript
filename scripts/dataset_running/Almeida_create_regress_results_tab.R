rm(list = ls(all.names = TRUE))

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/")

almeida_regress_out <- readRDS(file = "Almeida_2019_regress_specificity_output/combined_output.rds")

names(almeida_regress_out)[which(names(almeida_regress_out) == "ERP002061")] <- "ERP002061 | Obesity 1"
names(almeida_regress_out)[which(names(almeida_regress_out) == "ERP002061")] <- "ERP003612 | Obesity 2"
names(almeida_regress_out)[which(names(almeida_regress_out) == "ERP012177")] <- "ERP012177 | CRC"


for (dataset in names(almeida_regress_out)) {
  for (func_level in names(almeida_regress_out[[dataset]])) {
    almeida_regress_out[[dataset]][[func_level]]$func <- rownames(almeida_regress_out[[dataset]][[func_level]])
    almeida_regress_out[[dataset]][[func_level]]$func_type <- func_level
    almeida_regress_out[[dataset]][[func_level]]$dataset <- dataset
  }
}

almeida_regress_out_tmp <- list()

for (dataset in names(almeida_regress_out)) {
  almeida_regress_out_tmp[[dataset]] <- do.call(rbind, almeida_regress_out[[dataset]])
}

almeida_regress_combined <- do.call(rbind, almeida_regress_out_tmp)

almeida_regress_combined_sig <- almeida_regress_combined[which(almeida_regress_combined$BH < 0.25), ]

rownames(almeida_regress_combined_sig) <- NULL

write.table(x = almeida_regress_combined_sig,
            file = "/home/gavin/github_repos/POMS_manuscript/key_datafiles/almeida_regress_results.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
