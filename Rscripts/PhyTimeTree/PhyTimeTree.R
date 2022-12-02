# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", dependencies = TRUE)
if (!requireNamespace("ggtree", quietly = TRUE)) BiocManager::install("ggtree", dependencies = TRUE)
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn", dependencies = TRUE)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("tidyverse", dependencies = TRUE)
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite", dependencies = TRUE)

library("ggplot2")
library("ggtree")
library("phangorn")
library("readxl")
library("svglite")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

tree <- read.tree("PhyTimeTree_Input1.nwk")

metadata <- read_excel('PhyTimeTree_Input2.xlsx')

tree <- midpoint(tree)
tree <- reorder(tree)
tree <- ladderize(tree)
tree <- ggtree(tree, size = 0.4, color = "#8c8fae", ladderize = FALSE) + geom_rootedge(rootedge = 0.00001, size = 0.4)
tree <- tree %<+% metadata + geom_tippoint(aes(color = group2), shape = 19, size = 6, show.legend = TRUE) + 
  scale_color_manual(values = c("#660033", "#9A6348", "#D79B7D", "#C0C741", "#647D34", "#E4943A", "#CC3333", "#D26471", "#70377F", "#7EC4C1", "#34859D", "#17434B", "#441A3F", "#584563", "#C2C2DA"), na.translate = FALSE)
tree
ggsave(filename = "PhyTimeTree_Output.pdf", dpi=600, width=8, height=4, units="in")