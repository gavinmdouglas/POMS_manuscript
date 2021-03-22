rm(list=ls(all.names=TRUE))

setwd("~/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset")

almeida_taxa_umgs <- read.table("taxonomy/taxonomy_umgs.tab.gz", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_umgs) <- almeida_taxa_umgs$MAG_ID
almeida_taxa_umgs <- almeida_taxa_umgs[, -which(colnames(almeida_taxa_umgs) %in% c("UMGS_ID", "MAG_ID"))]

almeida_taxa_hgr <- read.table("taxonomy/taxonomy_hgr.tab.gz", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_hgr) <- almeida_taxa_hgr$Genome
almeida_taxa_hgr <- almeida_taxa_hgr[, -which(colnames(almeida_taxa_hgr) == "Genome")]

almeida_taxa <- rbind(almeida_taxa_hgr, almeida_taxa_umgs)

saveRDS(object = almeida_taxa, file = "~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/taxa_table.rds")

