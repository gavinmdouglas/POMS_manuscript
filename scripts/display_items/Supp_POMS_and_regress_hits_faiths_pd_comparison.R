rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(ggbeeswarm)

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits")

almeida_output <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits/case.control_Almeida.2019.tsv",
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE)

almeida_output[which(almeida_output$Approach == "Phylo. regress. (specificity)"), "Approach"] <- "Phylo. regress.\n(specificity)"
almeida_output_edges1[which(almeida_output_edges1$Approach == "Phylo. regress. (specificity)"), "Approach"] <- "Phylo. regress.\n(specificity)"

almeida_output$Approach <- factor(almeida_output$Approach, levels = c("POMS", "Phylo. regress.\n(specificity)"))

almeida_output[which(almeida_output$Datatype == "kos"), "Datatype"] <- "KOs"
almeida_output[which(almeida_output$Datatype == "pathways"), "Datatype"] <- "Pathways"
almeida_output[which(almeida_output$Datatype == "modules"), "Datatype"] <- "Modules"
almeida_output$Datatype <- factor(almeida_output$Datatype, levels = c("KOs", "Pathways", "Modules"))

almeida_output <- almeida_output[which(almeida_output$Dataset != "Colorectal cancer"), ]
almeida_output <- almeida_output[which(almeida_output$Dataset != "Obesity 1" | almeida_output$Datatype != "Pathways"), ]

alemida_PD <- ggplot(data = almeida_output, aes(x = Datatype, y = PD, fill = Approach)) +
                      geom_quasirandom(colour = "grey", dodge.width=0.75) +
                      geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.5) +
                      theme_bw() +
                      facet_grid( ~ Dataset, drop = TRUE, scales = "free_x") +
                      ylab("Faith's Phylogenetic Diversity") +
                      xlab("") +
                      scale_fill_manual(values = c("#1f78b4", "blue"), drop = TRUE) +
                      ggtitle("Human case-control datasets") +
  theme(plot.title = element_text(hjust = 0.5))



TARA_output <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/data/results/Faiths_pd_sig_hits/Tara_Oceans.tsv",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)

TARA_output[which(TARA_output$Approach == "Phylo. regress. (specificity)"), "Approach"] <- "Phylo. regress.\n(specificity)"

TARA_output$Approach <- factor(TARA_output$Approach, levels = c("POMS", "Phylo. regress.\n(specificity)"))

TARA_output[which(TARA_output$Env_factor == "Mean_Salinity"), "Env_factor"] <- "Mean salinity"
TARA_output[which(TARA_output$Env_factor == "PO4"), "Env_factor"] <- "Phosphate"
TARA_output$Env_factor <- factor(TARA_output$Env_factor, levels = c("Mean salinity", "Phosphate"))

TARA_output[which(TARA_output$Datatype == "ko"), "Datatype"] <- "KOs"
TARA_output[which(TARA_output$Datatype == "pathway"), "Datatype"] <- "Pathways"
TARA_output[which(TARA_output$Datatype == "module"), "Datatype"] <- "Modules"
TARA_output$Datatype <- factor(TARA_output$Datatype, levels = c("KOs", "Pathways", "Modules"))

TARA_output <- TARA_output[which(TARA_output$Env_factor != "Mean salinity" | TARA_output$Datatype != "KOs"), ]
TARA_output <- TARA_output[which(TARA_output$Env_factor != "Phosphate" | TARA_output$Datatype != "Modules"), ]

TARA_PD <- ggplot(data = TARA_output, aes(x = Datatype, y = PD, fill = Approach)) +
  geom_quasirandom(colour = "grey", dodge.width=0.75) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.5) +
  theme_bw() +
  facet_grid( ~ Env_factor, drop = TRUE, scales = "free_x") +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("") +
  scale_fill_manual(values = c("#1f78b4", "blue"), drop = TRUE) +
  ggtitle("Tara Oceans") +
  theme(plot.title = element_text(hjust = 0.5))


Faiths_pd_boxplots <- plot_grid(TARA_PD, alemida_PD, labels = c('a', 'b'), nrow = 2)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_faiths_pd_sig_hits.pdf",
       plot = Faiths_pd_boxplots,
       device = "pdf",
       width = 6,
       height = 6,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_faiths_pd_sig_hits.png",
       plot = Faiths_pd_boxplots,
       device = "png",
       width = 6,
       height = 6,
       dpi = 300)


wilcox.test(almeida_output[which(almeida_output$Approach == "POMS"), "PD"], almeida_output[which(almeida_output$Approach != "POMS"), "PD"])

mean(almeida_output[which(almeida_output$Approach == "POMS" ), "PD"])
mean(almeida_output[which(almeida_output$Approach != "POMS"), "PD"])
sd(almeida_output[which(almeida_output$Approach == "POMS" ), "PD"])
sd(almeida_output[which(almeida_output$Approach != "POMS"), "PD"])


wilcox.test(TARA_output[which(TARA_output$Approach == "POMS"), "PD"], TARA_output[which(TARA_output$Approach != "POMS"), "PD"])

mean(TARA_output[which(TARA_output$Approach == "POMS" ), "PD"])
mean(TARA_output[which(TARA_output$Approach != "POMS"), "PD"])
sd(TARA_output[which(TARA_output$Approach == "POMS" ), "PD"])
sd(TARA_output[which(TARA_output$Approach != "POMS"), "PD"])
