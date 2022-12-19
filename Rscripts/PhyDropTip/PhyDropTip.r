# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("ape", quietly = TRUE))
    install.packages("ape", dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", dependencies = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", dependencies = TRUE)
if (!requireNamespace("rlang", quietly = TRUE))
    install.packages("rlang", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE))
    install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("treeio", quietly = TRUE))
    BiocManager::install("treeio", dependencies = TRUE)

library("ape")
library("data.table")
library("dplyr")
library("rlang")
library("rstudioapi")
library("treeio")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

tree <- read.tree("PhyDropTip_Input.nwk")

drop_tip <- c("EPI_ISL_11026046", "EPI_ISL_10101858", "EPI_ISL_9671457")

new_tree <- drop.tip(tree, drop_tip, trim.internal = TRUE)

write.tree(new_tree, file = "PhyDropTip_Ouput.nwk", append = FALSE)
