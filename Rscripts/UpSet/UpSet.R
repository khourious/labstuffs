# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", dependencies = TRUE)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR", dependencies = TRUE)

library("ComplexHeatmap")
library("ggplot2")
library("rstudioapi")
library("UpSetR")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

list <- read.csv("UpSet_Input.csv", check.names = FALSE)
list

matrix <- t(table(stack(list)))
matrix <- t(matrix[, -1])
matrix

write.csv(matrix, "UpSet_Output1.csv", row.names = TRUE)

pdf(file = "UpSet_Output2.pdf")
upset(fromList(list), 
      nintersects = NA, 
      # sets = c("SAMPLE_01", "SAMPLE_02", "SAMPLE_03", "SAMPLE_04", "SAMPLE_05"),
      keep.order = TRUE,
      # mainbar.y.label = "Intersection Size",
      mainbar.y.max = 50,
      # matrix.color = "#ff7b00",
      main.bar.color = "#e71837",
      sets.bar.color = "#27357e",
      # sets.x.label = "Set Size",
      point.size = 4.0,
      line.size = 1.0,
      # mb.ratio = c(0.6, 0.4),
      order.by = "freq",
      number.angles = 0,
      text.scale = c(1, 1, 1, 1, 1, 1))
dev.off()