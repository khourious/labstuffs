# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", dependencies = TRUE)
if (!requireNamespace("ggtree", quietly = TRUE))
    BiocManager::install("ggtree", dependencies = TRUE)
if (!requireNamespace("phangorn", quietly = TRUE))
    install.packages("phangorn", dependencies = TRUE)
if (!requireNamespace("readxl", quietly = TRUE))
    install.packages("tidyverse", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE))
    install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("svglite", quietly = TRUE))
    install.packages("svglite", dependencies = TRUE)

library("ggplot2")
library("ggtree")
library("phangorn")
library("readxl")
library("rstudioapi")
library("svglite")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

input_file_1 <- "AKO_More80Cov.aln.edited.nwk"
input_file_2 <- "AKO_More80Cov.metadata.xlsx"
output_file <- "AKO_More80Cov.aln.edited.svg"

tree <- read.tree(input_file_1)
metadata <- read_excel(input_file_2)

tree <- midpoint(tree)
tree <- reorder(tree)
tree <- ladderize(tree)
tree <- ggtree(tree, size = .4, color = "#8c8fae", ladderize = FALSE) +
  geom_rootedge(rootedge = .00001, size = .4)
tree <- tree %<+% metadata + geom_tippoint(aes(color = group2), shape = 19, size = 6,
                                            show.legend = TRUE) +
    scale_color_manual(na.translate = FALSE, values = c(
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
        "#C2C2DA"))

ggsave(filename = output_file, dpi = 600, width = 8, height = 4, units = "in")
