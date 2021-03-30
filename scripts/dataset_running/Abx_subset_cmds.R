### Get SRA run ids of samples from timepoints 0 and 7 only - exclude 90 day samples to save disk space and time.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/Raymond_2016/")

abx_sra_table <- read.table("SraRunTable.txt", header=TRUE, sep=",",
                            stringsAsFactors = FALSE, comment.char = "")

abx_sra_table_no90 <- abx_sra_table[which(abx_sra_table$Timepoint != 90), ]

write.table(x = abx_sra_table_no90$Run, file="sra_runs_to_download.txt", col.names = FALSE, row.names = FALSE, quote=FALSE)
