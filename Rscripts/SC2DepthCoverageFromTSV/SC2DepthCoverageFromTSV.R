# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("cowplot", quietly = TRUE))
    install.packages("cowplot", dependencies = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", dependencies = TRUE)
if (!requireNamespace("patchwork", quietly = TRUE))
    install.packages("patchwork", dependencies = TRUE)
if (!requireNamespace("plyr", quietly = TRUE))
    install.packages("plyr", dependencies = TRUE)
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("tidyverse", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("svglite", quietly = TRUE))
  install.packages("svglite", dependencies = TRUE)

library("cowplot")
library("dplyr")
library("ggplot2")
library("patchwork")
library("plyr")
library("readr")
library("rstudioapi")
library("svglite")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

# TSV file: samtools depth -a *.sorted.bam > *.tsv
input_file <- "103900.depth.tsv"
output_file <- "103900.coverage.svg"

df <- read.delim(input_file, sep = "\t", header = FALSE)
depth_cov_df <- data.frame(position = df$V2, depth = df$V3)

# https://www.ncbi.nlm.nih.gov/nuccore/1798174254
# https://doi.org/10.1038/s41586-020-2286-9
map_lv1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                "UTR", "5'UTR", 1, 265,
                "ORFs", "ORF1a", 266, 13468,
                "ORFs", "ORF1b", 13468, 21555,
                "Structural proteins", "S", 21563, 25384,
                "Structural proteins", "E", 26245, 26472,
                "Structural proteins", "M", 26523, 27191,
                "Structural proteins", "N", 28274, 29533,
                "UTR", "3'UTR", 29675, 29903)
map_lv2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                "NSPs", "nsp1", 266, 805,
                "NSPs", "nsp2", 806, 2719,
                "NSPs", "nsp3", 2720, 8554,
                "NSPs", "nsp4", 8555, 10054,
                "NSPs", "nsp5", 10055, 10972,
                "NSPs", "nsp6", 10973, 11842,
                "NSPs", "nsp8", 12092, 12685,
                "NSPs", "nsp10", 13025, 13441,
                "NSPs", "nsp12", 13442, 16236,
                "NSPs", "nsp13", 16237, 18039,
                "NSPs", "nsp14", 18040, 19620,
                "NSPs", "nsp15", 19621, 20658,
                "NSPs", "nsp16", 20659, 21552,
                "Accessory factors", "3ab", 25393, 26220,
                "Accessory factors", "6", 27202, 27387,
                "Accessory factors", "7a", 27394, 27759,
                "Accessory factors", "7b", 27756, 27887,
                "Accessory factors", "8", 27894, 28259,
                "Accessory factors", "9b", 28284, 28577,
                "Accessory factors", "10", 29558, 28674)
map_lv3 <- tribble(~"class", ~"gene", ~"start", ~"end",
                "NSPs", "nsp7", 11843, 12091,
                "NSPs", "nsp9", 12686, 13024,
                "NSPs", "nsp11", 13442, 13480)

depth_cov_plot <- ggplot() +
  geom_line(data = depth_cov_df, aes(x = position, y = depth),
            linewidth = .4, colour = "#000000") +
  labs(title = input_file,
       y = "Per base coverage (x)", x = NULL) +
  scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                     expand = expansion(0, 0), limits = c(0, 30000)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  theme_light(base_size = 10) +
  scale_y_log10(
#    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(hjust = 1, size = 8)) +
  geom_hline(yintercept = 10, linetype = "dotted", colour = "#5A5A5A")

map_lv1_plot <- map_lv1 %>% ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
            linewidth = .2, colour = "#000000", alpha = .3) +
  geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
  scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
  theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "UTR" = "#2F67CD",
    "ORFs" = "#FE0B12",
    "Structural proteins" = "#11961B"))

map_lv2_plot <- map_lv2 %>% ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
            linewidth = .2, colour = "#000000", alpha = .3) +
  geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
  scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
  theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "NSPs" = "#FF8A8E",
    "Accessory factors" = "#D860CF"))

map_lv3_plot <- map_lv3 %>% ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
            linewidth = .2, colour = "#000000", alpha = .3) +
  geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
  scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
  theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "NSPs" = "#FF8A8E"))

depth_cov_sars2 <- depth_cov_plot / map_lv1_plot / map_lv2_plot / map_lv3_plot +
  plot_layout(nrow = 4, heights = c(3, .4, .3, .3))
save_plot(output_file, depth_cov_sars2, base_height = 5, base_width = 16)
