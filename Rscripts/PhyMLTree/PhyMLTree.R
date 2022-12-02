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

tree <- read.tree("PhyMLTree_Input1.nwk")

metadata <- read_excel("PhyMLTree_Input2.xlsx")

tree <- ggtree(tree, mrsd="2022-03-14", as.Date=FALSE, size=0.4, color="#8c8fae") + theme_tree2() + geom_rootedge(rootedge = 0.01, size=0.4)
tree <- tree %<+% metadata + geom_tippoint(aes(color=group2), shape=19, size=4, show.legend=TRUE) +
  scale_color_manual(values=c("#660033", "#9A6348", "#D79B7D", "#C0C741", "#647D34", "#E4943A", "#CC3333", "#D26471", "#70377F", "#7EC4C1", "#34859D", "#17434B", "#441A3F", "#584563", "#C2C2DA"), na.translate=FALSE) +
  scale_x_ggtree(breaks=c(2021.6657534246576, 2021.7479452054795, 2021.8328767123287, 2021.9150684931508, 2022.0, 2022.0849315068492, 2022.1616438356164), labels=c("Sep/2021", "Oct/2021", "Nov/2021", "Dez/2021", "Jan/2022", "Fev/2022", "Mar/2022"))
tree
ggsave(filename="PhyMLPlot_Output.pdf", dpi=600, width=8, height=12, units="in")