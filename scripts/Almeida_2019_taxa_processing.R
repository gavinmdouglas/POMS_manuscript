rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/MAGs/Almeida2019/")

almeida_taxa_umgs <- read.table("taxonomy/taxonomy_umgs.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_umgs) <- almeida_taxa_umgs$MAG_ID
almeida_taxa_umgs <- almeida_taxa_umgs[, -which(colnames(almeida_taxa_umgs) %in% c("UMGS_ID", "MAG_ID"))]

almeida_taxa_hgr <- read.table("taxonomy/taxonomy_hgr.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_hgr) <- almeida_taxa_hgr$Genome
almeida_taxa_hgr <- almeida_taxa_hgr[, -which(colnames(almeida_taxa_hgr) == "Genome")]

almeida_taxa <- rbind(almeida_taxa_hgr, almeida_taxa_umgs)

saveRDS(object = almeida_taxa, file = "/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_POMS_output/taxa_table.rds")