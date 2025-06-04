# https://lpmor22.github.io/ | 2025-06-04

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load(cowplot, ggplot2, ggtree, phangorn, readxl, this.path, tidyr)

setwd(dirname(this.path()))

input1 <- "AKO_More80Cov.aln.edited.nwk"
input2 <- "AKO_More80Cov.aln.edited.xlsx"

tree <- read.tree(input1)
metadata <- read_excel(input2)

tree <- midpoint(tree)
tree <- reorder(tree)
tree <- ladderize(tree)
tree <- ggtree(tree, size = .4, color = "#8c8fae", ladderize = FALSE) +
  geom_rootedge(rootedge = .00001, size = .4)
tree <- tree %<+% metadata +
  geom_tippoint(aes(color = group2),
                shape = 19, size = 4, show.legend = TRUE) +
  scale_color_manual(na.translate = FALSE,
                       values = c(
                         "#660033",
                         "#9A6348",
                         "#D79B7D",
                         "#C0C741",
                         "#647D34",
                         "#E4943A",
                         "#CC3333",
                         "#D26471",
                         "#70377F",
                         "#7EC4C1",
                         "#34859D",
                         "#17434B",
                         "#441A3F",
                         "#584563",
                         "#C2C2DA"
                        )
  )

save_plot(sub("\\.nwk$", ".pdf", input1), tree, base_height = 4, base_width = 8)
save_plot(sub("\\.nwk$", ".png", input1), tree, base_height = 4, base_width = 8)