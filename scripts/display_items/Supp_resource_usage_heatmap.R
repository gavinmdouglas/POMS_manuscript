rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

resource_usage <- read.table("/home/gavin/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/resource_tracking_out_combined.tsv",
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)


resource_usage$Approach <- NA
resource_usage[grep("POMS", resource_usage$filename), "Approach"] <- "POMS"
resource_usage[grep("regress_specificity", resource_usage$filename), "Approach"] <- "Phylogenetic regression\n(specificity)"
resource_usage[grep("regress_sig.taxa", resource_usage$filename), "Approach"] <- "Phylogenetic regression\n(significant taxa)"


resource_usage$num_mags <- NA
resource_usage[grep("MAGs50", resource_usage$filename), "num_mags"] <- 50
resource_usage[grep("MAGs100", resource_usage$filename), "num_mags"] <- 100
resource_usage[grep("MAGs250", resource_usage$filename), "num_mags"] <- 250
resource_usage[grep("MAGs500", resource_usage$filename), "num_mags"] <- 500
resource_usage[grep("MAGs750", resource_usage$filename), "num_mags"] <- 750
resource_usage[grep("MAGs1000", resource_usage$filename), "num_mags"] <- 1000
resource_usage[grep("MAGs1250", resource_usage$filename), "num_mags"] <- 1250
resource_usage[grep("MAGs1595", resource_usage$filename), "num_mags"] <- 1595


resource_usage$replicate <- NA
resource_usage[grep("_rep1_", resource_usage$filename), "replicate"] <- 1
resource_usage[grep("_rep2_", resource_usage$filename), "replicate"] <- 2
resource_usage[grep("_rep3_", resource_usage$filename), "replicate"] <- 3

resource_usage$ncores <- NA
resource_usage[grep("_ncores1_", resource_usage$filename), "ncores"] <- 1
resource_usage[grep("_ncores5_", resource_usage$filename), "ncores"] <- 5
resource_usage[grep("_ncores10_", resource_usage$filename), "ncores"] <- 10

rownames(resource_usage) <- resource_usage$filename
resource_usage <- resource_usage[, -1]

resource_usage_averaged <- aggregate(formula = . ~ Approach + num_mags + ncores, data = resource_usage, FUN = mean)
resource_usage_averaged <- resource_usage_averaged[, which(colnames(resource_usage_averaged) != "replicate")]

#resource_usage[which(resource_usage$ncores == 1 & resource_usage$Approach == "POMS" & resource_usage$num_mags == 50), ]

resource_usage_averaged$max_gbytes <- resource_usage_averaged$max_kbytes / 1e6

resource_usage_averaged$num_mags_factor <- factor(as.character(resource_usage_averaged$num_mags),
                                                  levels = sort(unique(resource_usage_averaged$num_mags)))

resource_usage_averaged$ncores_factor <- factor(as.character(resource_usage_averaged$ncores),
                                                levels = sort(unique(resource_usage_averaged$ncores)))

elapsed_seconds_heatmap <- ggplot(data = resource_usage_averaged, aes(x = num_mags_factor, y = ncores_factor, fill = elapsed_seconds)) +
                                  geom_tile() +
                                  geom_text(label = formatC(resource_usage_averaged$elapsed_seconds, format = 'f', digits = 0)) +
                                  xlab("Number of metagenome-assembled genomes ") +
                                  ylab("Number of cores") +
                                  facet_wrap(. ~ Approach, nrow = 2) +
                                  scale_fill_continuous(name = "Elapsed\nseconds", low = "light blue", high = "red", limits=c(0, 250)) +
                                  theme_bw()


memory_usage_heatmap <- ggplot(data = resource_usage_averaged, aes(x = num_mags_factor, y = ncores_factor, fill = max_gbytes)) +
                                  geom_tile() +
                                  geom_text(label = formatC(resource_usage_averaged$max_gbytes, format = 'f', digits = 2)) +
                                  xlab("Number of metagenome-assembled genomes ") +
                                  ylab("Number of cores") +
                                  facet_wrap(. ~ Approach, nrow = 2) +
                                  scale_fill_continuous(name = "Memory\nusage\n(GB)", low = "light blue", high = "red", limits=c(0, 1)) +
                                  theme_bw()

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_elapsed_seconds_resource_usage.png",
       plot = elapsed_seconds_heatmap,
       device = "png",
       width = 7,
       height = 6,
       dpi = 300)

