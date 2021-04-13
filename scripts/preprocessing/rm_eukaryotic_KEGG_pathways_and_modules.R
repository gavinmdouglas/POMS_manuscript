# Filter out pathways and modules that are exclusively found in eukaryotes.

rm(list = ls(all.names = TRUE))

library(stringr)

setwd("~/github_repos/POMS_manuscript/data/KEGG_mappings/")

KEGG_genome_info_clean <- read.table("prepped/2021_04_12_KEGG_genome_taxonomy.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")


KEGG_pathway_genome_links <- read.table("prepped/2021_04_12_KEGG_genome_pathway_links.tsv.gz",
                                        header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

KEGG_pathway_ko <- read.table("prepped/2021_04_12_KEGG_ko_pathway_links.tsv.gz",
                              header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)



KEGG_module_genome_links <- read.table("prepped/2021_04_12_KEGG_genome_module_links.tsv.gz",
                                        header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

KEGG_module_ko <- read.table("prepped/2021_04_12_KEGG_ko_module_links.tsv.gz",
                              header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)



unique_pathways <- KEGG_pathway_genome_links$V3[-which(duplicated(KEGG_pathway_genome_links$V3))]

pathways_to_keep <- c()
pathways_to_rm <- c()

for (pathway in unique_pathways) {
  
  contributing_org <- KEGG_pathway_genome_links[which(KEGG_pathway_genome_links$V3 == pathway), "V1"]
  
  contributing_superkingdom <- KEGG_genome_info_clean[which(KEGG_genome_info_clean$org %in% contributing_org), "superkingdom"]
  
  if ((!"Bacteria" %in%  contributing_superkingdom) & (!"Archaea" %in%  contributing_superkingdom)  & (!"Eukaryota" %in%  contributing_superkingdom)) {
    stop("Pathway not in any of the three superkingdoms")
  }
  
  if ("Bacteria" %in%  contributing_superkingdom | "Archaea" %in%  contributing_superkingdom) {
    pathways_to_keep <- c(pathways_to_keep, pathway)
  } else {
    pathways_to_rm <- c(pathways_to_rm, pathway)
  }

}


KEGG_pathway_ko_pro <- KEGG_pathway_ko[pathways_to_keep, , drop = FALSE]
write.table(KEGG_pathway_ko_pro, file = "prepped/2021_04_12_KEGG_ko_pathway_links_prokaryotic.tsv", col.names = FALSE, row.names = TRUE, quote = FALSE, sep = "\t")


unique_modules <- KEGG_module_genome_links$V3[-which(duplicated(KEGG_module_genome_links$V3))]

modules_to_keep <- c()
modules_to_rm <- c()

for (module in unique_modules) {
  
  contributing_org <- KEGG_module_genome_links[which(KEGG_module_genome_links$V3 == module), "V1"]
  
  contributing_superkingdom <- KEGG_genome_info_clean[which(KEGG_genome_info_clean$org %in% contributing_org), "superkingdom"]
  
  if ((!"Bacteria" %in%  contributing_superkingdom) & (!"Archaea" %in%  contributing_superkingdom)  & (!"Eukaryota" %in%  contributing_superkingdom)) {
    stop("Module not in any of the three superkingdoms")
  }
  
  if ("Bacteria" %in%  contributing_superkingdom | "Archaea" %in%  contributing_superkingdom) {
    modules_to_keep <- c(modules_to_keep, module)
  } else {
    modules_to_rm <- c(modules_to_rm, module)
  }
  
}

KEGG_module_ko_pro <- KEGG_module_ko[modules_to_keep, , drop = FALSE]
write.table(KEGG_module_ko_pro, file = "prepped/2021_04_12_KEGG_ko_module_links_prokaryotic.tsv", col.names = FALSE, row.names = TRUE, quote = FALSE, sep = "\t")
