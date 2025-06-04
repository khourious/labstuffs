# https://lpmor22.github.io/ | 2025-06-03

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load(cowplot, dplyr, ggplot2, patchwork, plyr, readr, this.path)

setwd(dirname(this.path()))

# samtools depth -a *.sorted.bam > *.tsv
input <- "103900.depth.tsv"

df <- read.delim(input, sep = "\t", header = FALSE)
depth_cov_df <- data.frame(position = df$V2, depth = df$V3)

depth_cov_plot <- ggplot() +
  geom_line(data = depth_cov_df, aes(x = position, y = depth),
            linewidth = .4, colour = "#000000") +
  labs(title = input,
       y = "Per base coverage (x)", x = NULL) +
  scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                     expand = expansion(0, 0), limits = c(0, 30000)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  theme_light(base_size = 10) +
  scale_y_log10(
    # labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  theme(plot.title = element_text(hjust = .5, size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(hjust = 1, size = 8),
        plot.margin = margin(10, 10, 10, 10)) +
  geom_hline(yintercept = 10, linetype = "dotted", colour = "#5A5A5A")

save_plot(sub("\\.tsv$", ".pdf", input), depth_cov_plot, base_height = 4, base_width = 16)
save_plot(sub("\\.tsv$", ".png", input), depth_cov_plot, base_height = 4, base_width = 16)