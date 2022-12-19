# @lpmor22 | https://lpmor22.github.io/

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
if (!requireNamespace("tidyr", quietly = TRUE))
    install.packages("tidyr", dependencies = TRUE)

library("cowplot")
library("dplyr")
library("ggplot2")
library("lubridate")
library("plyr")
library("readr")
library("rstudioapi")
library("tidyr")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

# https://www.epicov.org/epi3/frontend
# EpiCov -> Search
gisaid_tar <- list.files(pattern = "*.tar", full.names = TRUE)
ldply(.data = gisaid_tar, .fun = untar)

gisaid_metadata <- list.files(pattern = "*.metadata.tsv", full.names = TRUE)
input = ldply(gisaid_metadata, read_tsv, col_types = cols(.default = "c"))
input$length <- as.numeric(input$length)
filt <- input[input$length > 25000, ]

filt$epiym <- format(as.Date(filt$date), "%Y-%m")

# https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
# 2022 Dec 15
sc2_variants <- list(
    "B.1.1" = "B.1.1",
    "B.1.1.28" = "B.1.1.28",
    "B.1.1.33" = "B.1.1.33",
    "B.1.1.7 (Alpha)" = c("B.1.1.7", unique(input$pangolin_lineage[grep("^Q\\.", input$pangolin_lineage)])),
    "B.1.351 (Beta)" = c("B.1.351", unique(input$pangolin_lineage[grep("^B.1.351\\.", input$pangolin_lineage)])),
    "P.1+P.1.* (Gamma)" = c("P.1", unique(input$pangolin_lineage[grep("^P\\.1\\.", input$pangolin_lineage)])),
    "B.1.617.2+AY.* (Delta)" = c("B.1.617.2", unique(input$pangolin_lineage[grep('^AY\\.',input$pangolin_lineage)])),
    "B.1.427+B.1.429 (Epsilon)" = c("B.1.427", "B.1.429"),
    "B.1.525 (Eta)" = "B.1.525",
    "B.1.526 (Iota)" = "B.1.526",
    "B.1.617.1 (Kappa)" = "B.1.617.1",
    "P.2 (Zeta)" = "P.2",
    "B.1.621+B.1.621.1 (Mu)" = c("B.1.621", "B.1.621.1"),
    "BA.1+BA.1.* (Omicron)" = c("B.1.1.529", "BA.1", unique(input$pangolin_lineage[grep("^BA\\.1\\.", input$pangolin_lineage)])),
    "BA.2+BA.2.* (Omicron)" = c("BA.2", unique(input$pangolin_lineage[grep("^BA\\.2\\.", input$pangolin_lineage)])),
    "BA.3+BA.3.* (Omicron)" = c("BA.3", unique(input$pangolin_lineage[grep("^BA\\.3\\.", input$pangolin_lineage)])),
    "BA.4+BA.4.* (Omicron)" = c("BA.4", unique(input$pangolin_lineage[grep("^BA\\.4\\.", input$pangolin_lineage)])),
    "BA.5+BA.5.* (Omicron)" = c("BA.5", unique(input$pangolin_lineage[grep("^BA\\.5\\.", input$pangolin_lineage)])),
    "BQ.1+BQ.1.* (Omicron)" = c("BQ.1", unique(input$pangolin_lineage[grep("^BQ\\.1\\.", input$pangolin_lineage)])),
    "BE.9 (Omicron)" = "BE.9"
)

sc2_variants_greek <- structure(c(rep(names(sc2_variants),
                                      sapply(sc2_variants, length))),
                                .Names = c(unlist(sc2_variants)))

filt$greek <- ifelse(filt$pangolin_lineage %in% names(sc2_variants_greek),
                     sc2_variants_greek[filt$pangolin_lineage], "Others")

sc2_epiym <- ddply(filt, .(filt$epiym, filt$greek), nrow)
names(sc2_epiym) <- c("epiym", "variants", "n")

ggp <- ggplot(sc2_epiym, aes(x = epiym, y = n, fill = variants)) +
  geom_col(position = position_fill(reverse = TRUE), width = .95, alpha = .6) +
  labs(x = NULL, y = "Frequency", fill = "Pangolin Lineages") +
  scale_y_continuous(labels = scales::percent, expand = expansion(0, .05)) +
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
        axis.text.y = element_text(hjust = 1, size = 5),
        legend.text = element_text(size = 5), legend.title = element_text(size = 6),
        legend.position = "top", legend.key.size = unit(.3, "cm"))

save_plot("SC2ProportionalStackedBar_Output.svg", ggp, base_height = 2.5, base_width = 6)