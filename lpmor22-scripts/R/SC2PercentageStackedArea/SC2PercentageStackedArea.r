# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("areaplot", quietly = TRUE))
    install.packages("areaplot", dependencies = TRUE)
if (!requireNamespace("cowplot", quietly = TRUE))
    install.packages("cowplot", dependencies = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", dependencies = TRUE)
if (!requireNamespace("lubridate", quietly = TRUE))
    install.packages("tidyverse", dependencies = TRUE)
if (!requireNamespace("plyr", quietly = TRUE))
    install.packages("plyr", dependencies = TRUE)
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("tidyverse", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE))
    install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("svglite", quietly = TRUE))
  install.packages("svglite", dependencies = TRUE)
if (!requireNamespace("tidyr", quietly = TRUE))
    install.packages("tidyr", dependencies = TRUE)

library("areaplot")
library("cowplot")
library("dplyr")
library("ggplot2")
library("lubridate")
library("plyr")
library("readr")
library("rstudioapi")
library("svglite")
library("tidyr")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

# https://www.epicov.org/epi3/frontend
# EpiCov -> Search
input_file <- list.files(pattern = "*.tar", full.names = TRUE)
ldply(.data = input_file, .fun = untar)

output_file <- "Sars2_PangoLineagesArea_Bahia_20221218.svg"

gisaid_metadata <- list.files(pattern = "*.metadata.tsv", full.names = TRUE)
gisaid_metadata_df <- ldply(gisaid_metadata, read_tsv, col_types = cols(.default = "c"))
gisaid_metadata_df$length <- as.numeric(gisaid_metadata_df$length)

gisaid_filt <- gisaid_metadata_df[gisaid_metadata_df$length > 25000, ]

gisaid_filt$epiym <- format(as.Date(gisaid_filt$date), "%Y-%m")

# https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
# 2020 Dec 15
sc2_variants <- list(
    "B.1.1" = "B.1.1",
    "B.1.1.28" = "B.1.1.28",
    "B.1.1.33" = "B.1.1.33",
    "B.1.1.7 (Alpha)" = c("B.1.1.7", unique(gisaid_filt$pangolin_lineage[grep("^Q\\.", gisaid_filt$pangolin_lineage)])),
    "B.1.351 (Beta)" = c("B.1.351", unique(gisaid_filt$pangolin_lineage[grep("^B.1.351\\.", gisaid_filt$pangolin_lineage)])),
    "P.1+P.1.* (Gamma)" = c("P.1", unique(gisaid_filt$pangolin_lineage[grep("^P\\.1\\.", gisaid_filt$pangolin_lineage)])),
    "B.1.617.2+AY.* (Delta)" = c("B.1.617.2", unique(gisaid_filt$pangolin_lineage[grep('^AY\\.', gisaid_filt$pangolin_lineage)])),
    "B.1.427+B.1.429 (Epsilon)" = c("B.1.427", "B.1.429"),
    "B.1.525 (Eta)" = "B.1.525",
    "B.1.526 (Iota)" = "B.1.526",
    "B.1.617.1 (Kappa)" = "B.1.617.1",
    "P.2 (Zeta)" = "P.2",
    "B.1.621+B.1.621.1 (Mu)" = c("B.1.621", "B.1.621.1"),
    "BA.1+BA.1.* (Omicron)" = c("B.1.1.529", "BA.1", unique(gisaid_filt$pangolin_lineage[grep("^BA\\.1\\.", gisaid_filt$pangolin_lineage)])),
    "BA.2+BA.2.* (Omicron)" = c("BA.2", unique(gisaid_filt$pangolin_lineage[grep("^BA\\.2\\.", gisaid_filt$pangolin_lineage)])),
    "BA.3+BA.3.* (Omicron)" = c("BA.3", unique(gisaid_filt$pangolin_lineage[grep("^BA\\.3\\.", gisaid_filt$pangolin_lineage)])),
    "BA.4+BA.4.* (Omicron)" = c("BA.4", unique(gisaid_filt$pangolin_lineage[grep("^BA\\.4\\.", gisaid_filt$pangolin_lineage)])),
    "BA.5+BA.5.* (Omicron)" = c("BA.5", unique(gisaid_filt$pangolin_lineage[grep("^BA\\.5\\.", gisaid_filt$pangolin_lineage)])),
    "BQ.1+BQ.1.* (Omicron)" = c("BQ.1", unique(gisaid_filt$pangolin_lineage[grep("^BQ\\.1\\.", gisaid_filt$pangolin_lineage)])),
    "BE.9 (Omicron)" = "BE.9")

sc2_variants_greek <- structure(c(rep(names(sc2_variants),
                                      sapply(sc2_variants, length))),
                                .Names = c(unlist(sc2_variants)))

gisaid_filt$greek <- ifelse(gisaid_filt$pangolin_lineage %in% names(sc2_variants_greek),
                     sc2_variants_greek[gisaid_filt$pangolin_lineage], "Others")

sc2_epiym <- ddply(gisaid_filt, .(gisaid_filt$epiym, gisaid_filt$greek), nrow, .drop = FALSE)
names(sc2_epiym) <- c("epiym", "variants", "n")
sc2_epiym <- ddply(sc2_epiym, .(epiym), transform, percentage = n/sum(n))

sc2_epiym_plot <- ggplot(sc2_epiym, aes(x = epiym, y = percentage)) +
  geom_area(aes(fill = variants, group = variants), position = position_fill(reverse = TRUE),
            linewidth = .4, colour = "#ffffff",  alpha = .6) +
  labs(x = NULL, y = "Frequency", fill = "Pangolin Lineages") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(0, .05)) +
  scale_fill_manual(values = c(
      "B.1.1" = "#FB6B5B",
      "B.1.1.28" = "#E49074",
      "B.1.1.33" = "#FD9EA2",
      "B.1.1.7 (Alpha)" = "#DF2A8E",
      "P.2 (Zeta)" = "#872B15",
      "P.1+P.1.* (Gamma)" = "#11961B",
      "B.1.617.2+AY.* (Delta)" = "#6F0F84",
      "BA.1+BA.1.* (Omicron)" = "#EECB3A",
      "BA.2+BA.2.* (Omicron)" = "#FE0B12",
      "BA.4+BA.4.* (Omicron)" = "#D860CF",
      "BA.5+BA.5.* (Omicron)" = "#2F67CD",
      "BQ.1+BQ.1.* (Omicron)" = "#FF7F07",
      "BE.9 (Omicron)" = "#73FBFD",
      "Others" = "#999999"), 
                    breaks = c(
                        "B.1.1",
                        "B.1.1.28",
                        "B.1.1.33",
                        "B.1.1.7 (Alpha)",
                        "P.2 (Zeta)",
                        "P.1+P.1.* (Gamma)",
                        "B.1.617.2+AY.* (Delta)",
                        "BA.1+BA.1.* (Omicron)",
                        "BA.2+BA.2.* (Omicron)",
                        "BA.4+BA.4.* (Omicron)",
                        "BA.5+BA.5.* (Omicron)",
                        "BQ.1+BQ.1.* (Omicron)",
                        "BE.9 (Omicron)",
                        "Others")) +
  theme_void() +
  theme(axis.title.x = element_text(size = 6), axis.title.y = element_text(angle = 90, size = 6),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = .8, size = 5),
        axis.text.y = element_text(hjust = 1.5, size = 5),
        legend.text = element_text(size = 5), legend.title = element_text(size = 6),
        legend.position = "top", legend.key.size = unit(.3, "cm"))

save_plot(output_file, sc2_epiym_plot, base_height = 3, base_width = 6)
