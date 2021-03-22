### Take in a cohort of healthy samples and randomly subsample them into two groups.
### The expectation should be that overall there will be no significant differences since
### the groupings are totally random, but the co-occurence patterns and compositionality
### will be realistic.

rm(list=ls(all.names=TRUE))

### First identify all healthy adult samples that should be used for this analysis.
setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt.gz",
                                  header=TRUE, sep="\t", row.names=1, comment.char = "", quote="", stringsAsFactors = FALSE)

almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Health.state == "Healthy"), ]
almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Antibiotics == "No"), ]
almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Age.group == "Adult"), ]

# Only retain samples from studies with at least 40 eligible samples.
study_summaries <- table(almeida_sample_info$Study)

almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study %in% names(which(study_summaries >= 40))), ]

# Remove one sample to get even number of samples.
almeida_sample_info <- almeida_sample_info[-1, ]

# Run 1000 replicates sampling all samples into either group1 or group2.
# Set seed and save these groupings as list in RDS for easy reproducibility.
set.seed(2751651)

almeida_healthy_random_groups <- list()
almeida_healthy_random_groups$group1 <- list()
almeida_healthy_random_groups$group2 <- list()

for(rep in 1:1000) {
  samples <- sample(rownames(almeida_sample_info))
  almeida_healthy_random_groups$group1[[rep]] <- samples[1:352]
  almeida_healthy_random_groups$group2[[rep]] <- samples[353:704]
}

saveRDS(object = almeida_healthy_random_groups, file = "~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/almeida_healthy_random_groups.rds")
