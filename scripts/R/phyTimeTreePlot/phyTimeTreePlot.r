# https://lpmor22.github.io/ | 2025-06-04

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load(cowplot, ggplot2, ggtree, phangorn, readxl, this.path, tidyr)

setwd(dirname(this.path()))

input1 <- "AKO_More80Cov.gisaid.aln.edited.timetree.nwk"
input2 <- "AKO_More80Cov.gisaid.aln.edited.timetree.xlsx"

tree <- read.tree(input1)
metadata <- read_excel(input2)

tree <- ggtree(tree, mrsd = "2022-03-14", as.Date = FALSE, size = .4, color = "#8c8fae") +
    theme_tree2() + geom_rootedge(rootedge = .01, linewidth = .4)
tree <- tree %<+% metadata + geom_tippoint(aes(color = group2), shape = 19, size = 4,
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
        "#C2C2DA")) +
    scale_x_ggtree(breaks = c(
        2021.6657534246576,
        2021.7479452054795,
        2021.8328767123287,
        2021.9150684931508,
        2022.0000000000000,
        2022.0849315068492,
        2022.1616438356164), labels = c(
            "Sep/2021",
            "Oct/2021",
            "Nov/2021",
            "Dez/2021",
            "Jan/2022",
            "Fev/2022",
            "Mar/2022"
    )
)

save_plot(sub("\\.nwk$", ".pdf", input1), tree, base_height = 12, base_width = 8)
save_plot(sub("\\.nwk$", ".png", input1), tree, base_height = 12, base_width = 8)