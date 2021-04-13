### Prep KEGG mapfiles downloaded from API on April 12, 2021.

### First process the description mapping files.
### Should look like this:
### ko00010 Glycolysis / Gluconeogenesis
### ko00020 Citrate cycle (TCA cycle)

setwd("~/github_repos/POMS_manuscript/data/KEGG_mappings/rawfiles/")

KO_raw_descrip <- read.table("2021_04_12_KEGG_ko_descrip.tsv.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
path_raw_descrip <- read.table("2021_04_12_KEGG_pathway_descrip.tsv.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
module_raw_descrip <- read.table("2021_04_12_KEGG_module_descrip.tsv.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

KO_raw_descrip$V1 <- gsub("^ko:", "", KO_raw_descrip$V1)
path_raw_descrip$V1 <- gsub("^path:", "", path_raw_descrip$V1)
path_raw_descrip$V1 <- gsub("^map", "ko", path_raw_descrip$V1)
module_raw_descrip$V1 <- gsub("^md:", "", module_raw_descrip$V1)

### e.g.
### ko05340 K10887,K10603,K10989
KO_to_path_raw <- read.table("2021_04_12_KEGG_ko_pathway_links.tsv.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
KO_to_path_raw <- KO_to_path_raw[-grep("^path:map", KO_to_path_raw$V2), ]
KO_to_path_raw$V1 <- gsub("^ko:", "", KO_to_path_raw$V1)
KO_to_path_raw$V2 <- gsub("^path:", "", KO_to_path_raw$V2)

unique_KEGG_path <- unique(KO_to_path_raw$V2)

KO_to_path <- data.frame(matrix(NA, nrow = length(unique_KEGG_path), ncol = 2))
rownames(KO_to_path) <- unique_KEGG_path

for (pathway in unique_KEGG_path) {
  contributing_KOs <- unique(KO_to_path_raw[which(KO_to_path_raw$V2 ==  pathway), "V1"])
  contributing_KOs <- paste(contributing_KOs, collapse = ",")
  KO_to_path[pathway, ] <- c(pathway, contributing_KOs)
}

# Remove pathways that don't map to KEGG orthologs.
path_raw_descrip <- path_raw_descrip[-which(!path_raw_descrip$V1 %in% KO_to_path$X1), ]

KO_to_module_raw <- read.table("2021_04_12_KEGG_ko_module_links.tsv.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
KO_to_module_raw$V1 <- gsub("^ko:", "", KO_to_module_raw$V1)
KO_to_module_raw$V2 <- gsub("^md:", "", KO_to_module_raw$V2)

unique_KEGG_module <- unique(KO_to_module_raw$V2)

KO_to_module <- data.frame(matrix(NA, nrow = length(unique_KEGG_module), ncol = 2))
rownames(KO_to_module) <- unique_KEGG_module

for (module in unique_KEGG_module) {
  contributing_KOs <- unique(KO_to_module_raw[which(KO_to_module_raw$V2 ==  module), "V1"])
  contributing_KOs <- paste(contributing_KOs, collapse = ",")
  KO_to_module[module, ] <- c(module, contributing_KOs)
}

# Remove modules that don't map to KEGG orthologs.
module_raw_descrip <- module_raw_descrip[-which(!module_raw_descrip$V1 %in% KO_to_module$X1), ]


# Write out all files
write.table(KO_raw_descrip, file = "../prepped/2021_04_12_KEGG_ko_descrip.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(path_raw_descrip, file = "../prepped/2021_04_12_KEGG_pathway_descrip.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(module_raw_descrip, file = "../prepped/2021_04_12_KEGG_module_descrip.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(KO_to_path, file = "../prepped/2021_04_12_KEGG_ko_pathway_links.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(KO_to_module, file = "../prepped/2021_04_12_KEGG_ko_module_links.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
