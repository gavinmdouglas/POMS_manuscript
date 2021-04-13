rm(list = ls(all.names = TRUE))

library(myTAI)
library(stringr)
library(taxize)

setwd("~/github_repos/POMS_manuscript/data/KEGG_mappings/")

KEGG_genome_info <- read.table("rawfiles/2021_04_12_KEGG_genome_descrip.tsv.gz",
                               header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")


KEGG_genome_info$V1 <- gsub("gn:", "", KEGG_genome_info$V1)

KEGG_genome_info_clean <- data.frame(matrix(NA, ncol = 10, nrow = nrow(KEGG_genome_info)))
colnames(KEGG_genome_info_clean) <- c("genome_id", "org", "taxid", "superkingdom", "class", "order", "family", "genus", "species", "strain")

weird_rows <- c()

org_ids <- c()
taxids <- c()

for (i in 1:nrow(KEGG_genome_info)) {
  
  info_split <- str_split(string = KEGG_genome_info$V2[i], pattern = ", |; ")[[1]]
  
  org_ids <- c(org_ids, info_split[1])
  
  number_only_indices <- grep('[a-zA-Z]', info_split, invert = TRUE)
  
  if (length(number_only_indices) == 1) {
    taxids <- c(taxids, info_split[number_only_indices])
  } else if (length(number_only_indices) == 0) {
    taxids <- c(taxids, NA)
  } else if (length(number_only_indices) > 1) {
    weird_rows <- c(weird_rows, i)
    taxids <- c(taxids, NA)
  }
  
}

KEGG_genome_info_clean$org <- org_ids
KEGG_genome_info_clean$taxid <- taxids
KEGG_genome_info_clean$genome_id <- KEGG_genome_info$V1

all_taxa_summary <- ncbi_get_taxon_summary(id = taxids)

for (j in 1:nrow(all_taxa_summary)) {
  taxon_summary <- taxonomy(organism = all_taxa_summary[j, "name"], db = "ncbi", output = "classification")
  
  taxon_row <- which(KEGG_genome_info_clean$taxid == all_taxa_summary[j, "uid"])
  
  for (rank in taxon_summary$rank) {
    if (rank %in% colnames(KEGG_genome_info_clean)) {
      KEGG_genome_info_clean[taxon_row, rank] <- paste(taxon_summary[which(taxon_summary$rank == rank), "name"], collapse = "|")
    }
  }

}


# Look specifically at cases where superkingdom was NA (sometimes they just need to be re-run). Restrict this to cases where there was actually a taxid.

NA_superkingdom_rows <- which(!is.na(KEGG_genome_info_clean$taxid) & is.na(KEGG_genome_info_clean$superkingdom))

for (j in NA_superkingdom_rows) {
  taxon_summary <- taxonomy(organism = all_taxa_summary[j, "name"], db = "ncbi", output = "classification")
  
  tax_split <- str_split(all_taxa_summary[j, "name"], " ")[[1]]
  
  while ((!"rank" %in% colnames(taxon_summary)) && (length(tax_split) > 1)) {
   
    tax_split <- tax_split[1:(length(tax_split) - 1)]
    
    shortened_tax <- paste(tax_split, collapse = " ")
    
    taxon_summary <- taxonomy(organism = shortened_tax, db = "ncbi", output = "classification")
    
  }
  
  taxon_row <- which(KEGG_genome_info_clean$taxid == all_taxa_summary[j, "uid"])
  
  for (rank in taxon_summary$rank) {
    if (rank %in% colnames(KEGG_genome_info_clean)) {
      KEGG_genome_info_clean[taxon_row, rank] <- paste(taxon_summary[which(taxon_summary$rank == rank), "name"], collapse = "|")
    }
  }
    
}

write.table(KEGG_genome_info_clean, file = "prepped/2021_04_12_KEGG_genome_taxonomy.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
