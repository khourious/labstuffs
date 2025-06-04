# https://lpmor22.github.io/ | 2025-06-04

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load(ape, data.table, dplyr, rlang, this.path, treeio)

setwd(dirname(this.path()))

input <- "AKO_More80Cov.gisaid.aln.edited.nwk"
drop_tip <- c("EPI_ISL_11026046", "EPI_ISL_10101858", "EPI_ISL_9671457")

tree <- read.tree(input)
tree <- drop.tip(tree, drop_tip, trim.internal = TRUE)

write.tree(tree, file = sub("\\.nwk$", "-droped.nwk", input), append = FALSE)