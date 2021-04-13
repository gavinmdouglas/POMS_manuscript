# Get general KEGG pathway and module id mappings to genomes rather than organism-specific

rm(list = ls(all.names = TRUE))

library(stringr)

setwd("~/github_repos/POMS_manuscript/data/KEGG_mappings/")

KEGG_genome_info_clean <- read.table("prepped/2021_04_12_KEGG_genome_taxonomy.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")


KEGG_pathway_genome_links <- read.table("rawfiles/2021_04_12_KEGG_genome_pathway_links.tsv.gz",
                                        header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

KEGG_module_genome_links <- read.table("rawfiles/2021_04_12_KEGG_genome_module_links.tsv.gz",
                                        header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

KEGG_pathway_genome_links$V1 <- gsub(pattern = "gn:", replacement = "", x = KEGG_pathway_genome_links$V1)
KEGG_module_genome_links$V1 <- gsub(pattern = "gn:", replacement = "", x = KEGG_module_genome_links$V1)


KEGG_pathway_id_split_raw <- str_split(string = KEGG_pathway_genome_links[, "V2"], pattern = '[a-z]0', simplify = TRUE)

KEGG_pathway_genome_links$pathway_clean <- paste("ko0", KEGG_pathway_id_split_raw[, 2], sep = "")

write.table(KEGG_pathway_genome_links, file = "prepped/2021_04_12_KEGG_genome_pathway_links.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")



KEGG_module_id_split_raw <- str_split(string = KEGG_module_genome_links[, "V2"], pattern = '[a-z]_M', simplify = TRUE)

KEGG_module_genome_links$module_clean <- paste("M", KEGG_module_id_split_raw[, 2], sep = "")

write.table(KEGG_module_genome_links, file = "prepped/2021_04_12_KEGG_genome_module_links.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
