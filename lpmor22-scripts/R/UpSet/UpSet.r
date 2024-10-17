# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = TRUE)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    BiocManager::install("ComplexHeatmap", dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE))
    install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("svglite", quietly = TRUE))
  install.packages("svglite", dependencies = TRUE)
if (!requireNamespace("UpSetR", quietly = TRUE))
    install.packages("UpSetR", dependencies = TRUE)

library("ComplexHeatmap")
library("ggplot2")
library("rstudioapi")
library("svglite")
library("UpSetR")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

input_file <- "Sars2Nextclade_Example.csv"
output_file_1 <- "Sars2Nextclade_Example_UpSetMatrix.csv"
output_file_2 <- "Sars2Nextclade_Example_UpSetPlot.svg"

df <- read.csv(input_file, check.names = FALSE)
df

matrix <- t(table(stack(df)))
matrix <- t(matrix[, -1])
matrix

write.csv(matrix, output_file_1, row.names = TRUE)

svg(output_file_2)
upset(fromList(df), nintersects = NA, keep.order = TRUE, point.size = 4, line.size = 1,
      # sets = c("SAMPLE_01", "SAMPLE_02", "SAMPLE_03", "SAMPLE_04", "SAMPLE_05"),
      # mainbar.y.label = "Intersection Size",
      mainbar.y.max = 50,
      # matrix.color = "#ff7b00",
      main.bar.color = "#e71837",
      sets.bar.color = "#27357e",
      # sets.x.label = "Set Size",
      # mb.ratio = c(0.6, 0.4),
      order.by = "freq", number.angles = 0, text.scale = c(1, 1, 1, 1, 1, 1))
dev.off()
