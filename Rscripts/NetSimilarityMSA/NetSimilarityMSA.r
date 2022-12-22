# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = TRUE)
if (!requireNamespace("DECIPHER", quietly = TRUE))
    BiocManager::install("DECIPHER", dependencies = TRUE)
if (!requireNamespace("igraph", quietly = TRUE))
    install.packages("igraph", dependencies = TRUE)
if (!requireNamespace("Matrix", quietly = TRUE))
    install.packages("Matrix", dependencies = TRUE)
if (!requireNamespace("rlang", quietly = TRUE))
    install.packages("rlang", dependencies = TRUE)
if (!requireNamespace("rstudioapi", quietly = TRUE))
    install.packages("rstudioapi", dependencies = TRUE)
if (!requireNamespace("smacof", quietly = TRUE))
    install.packages("smacof", dependencies = TRUE)
if (!requireNamespace("svglite", quietly = TRUE))
  install.packages("svglite", dependencies = TRUE)

library("DECIPHER")
library("igraph")
library("Matrix")
library("rlang")
library("rstudioapi")
library("smacof")
library("svglite")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

input_file <- "AKO_More80Cov.aln.edited.sorted.fasta"
output_file_1 <- "AKO_More80Cov_NetSim.svg"
output_file_2 <- "AKO_More80Cov_MatrixSim.csv"

fasta <- input_file

dna_string <- readDNAStringSet(fasta)
dna_string

dist <- DistanceMatrix(dna_string, type = "dist", includeTerminalGaps = FALSE,
                       penalizeGapLetterMatches = TRUE, correction = "none",
                       processors = 12, verbose = TRUE)

# sim <- sim2diss(dist, method = "corr", to.dist = FALSE)
# sim <- sim2diss(dist, method = "reverse", to.dist = FALSE)
# sim <- sim2diss(dist, method = "reciprocal", to.dist = FALSE)
# sim <- sim2diss(dist, method = "ranks", to.dist = FALSE)
sim <- sim2diss(dist, method = "exp", to.dist = FALSE)
# sim <- sim2diss(dist, method = "Gaussian", to.dist = FALSE)
# sim <- sim2diss(dist, method = "cooccurrence", to.dist = FALSE)
# sim <- sim2diss(dist, method = "gravity", to.dist = FALSE)
# sim <- sim2diss(dist, method = "confusion", to.dist = FALSE)
# sim <- sim2diss(dist, method = "transition", to.dist = FALSE)
# sim <- sim2diss(dist, method = "membership", to.dist = FALSE)
# sim <- sim2diss(dist, method = "probability", to.dist = FALSE)

net <- graph_from_adjacency_matrix(sim, mode = "directed", weighted = TRUE, diag = TRUE)
net <- simplify(net)

snet <- subgraph.edges(net, E(net)[E(net)$weight > 2], del = FALSE)

V(snet)$color <- "#C2C2DA"
V(snet)["101150_03"]$color <- "#660033"
V(snet)["101150_05"]$color <- "#660033"
V(snet)["102050_01"]$color <- "#9A6348"
V(snet)["102050_02"]$color <- "#9A6348"
V(snet)["103010_01"]$color <- "#D79B7D"
V(snet)["103010_04"]$color <- "#D79B7D"
V(snet)["103020_01"]$color <- "#C0C741"
V(snet)["103020_02"]$color <- "#C0C741"
V(snet)["103020_03"]$color <- "#C0C741"
V(snet)["103020_04"]$color <- "#C0C741"
V(snet)["103020_05"]$color <- "#C0C741"
V(snet)["109110_01"]$color <- "#647D34"
V(snet)["109110_04"]$color <- "#647D34"
V(snet)["127140_01"]$color <- "#E4943A"
V(snet)["127140_02"]$color <- "#E4943A"
V(snet)["207090_02"]$color <- "#CC3333"
V(snet)["207090_03"]$color <- "#CC3333"
V(snet)["213080_01"]$color <- "#D26471"
V(snet)["213080_05"]$color <- "#D26471"
V(snet)["219020_01"]$color <- "#70377F"
V(snet)["219020_02"]$color <- "#70377F"
V(snet)["220070_01"]$color <- "#7EC4C1"
V(snet)["220070_02"]$color <- "#7EC4C1"
V(snet)["220070_03"]$color <- "#7EC4C1"
V(snet)["220070_04"]$color <- "#7EC4C1"
V(snet)["229110_01"]$color <- "#34859D"
V(snet)["229110_02"]$color <- "#34859D"
V(snet)["412230_02"]$color <- "#17434B"
V(snet)["412230_04"]$color <- "#17434B"
V(snet)["426300_02"]$color <- "#441A3F"
V(snet)["426300_13"]$color <- "#441A3F"
V(snet)["430160_01"]$color <- "#584563"
V(snet)["430160_03"]$color <- "#584563"
# layout <- layout_in_circle
# layout <- layout_as_star
layout <- layout_with_kk
# layout <- layout_with_fr

svg(output_file_1)
plot(snet, layout = layout,
     edge.arrow.size = 0,
     edge.color = "#8C8FAE",
     edge.width = E(net)$weight/10,
     # edge.width = (edge_attr(snet)$weight)/100,
     # vertex.color = "white",
     vertex.label = NA,
     # vertex.label.cex = .5,
     # vertex.label.color = "#000000",
     # vertex.label.dist = 1,
     vertex.frame.color = "#000000",
     vertex.frame.width = .5,
     # vertex.size = igraph::degree(snet, mode = "all", loops = TRUE),
     vertex.size = 5,
     vertex.shape = "circle")
dev.off()

sim[upper.tri(sim, diag = FALSE)] <- ""
write.csv(sim, output_file_2)
