#!/usr/bin/env Rscript
# copap_modules.R -- find modules of co-evolving gene families from CoPAP / GLOOME output
#
# Idea (as in the CoPAP paper): keep the significant co-evolving pairs, build a network
# (families = nodes, significant pairs = edges) and cluster it with Transitivity Clustering,
# i.e. weighted cluster editing: choose clusters so that the total weight of edges you must
# add or delete to make every cluster a perfect clique is minimal.
#
# Usage:
#   Rscript copap_modules.R Coevolution.tsv annotation.tsv outdir [q_cutoff] [direction]
#     Coevolution.tsv : GLOOME/CoPAP table with columns POS1 POS2 r rate-bin p-value q-value direction
#     annotation.tsv  : tab-separated, first column = family position (as in POS1/POS2),
#                       further columns free (e.g. orthogroup, annotation, KEGG)
#     q_cutoff        : default 0.05  (use e.g. 0.1 or p-value for tiny toy data, see PVALUE below)
#     direction       : "co-occurrence" (default), "avoidance" or "both"
#   Optional environment variables:  EDGE_COL=q-value|p-value  SOLVER=heuristic|exact  MIN_SIZE=3
#
# Needs: R, igraph.  Optional: anticlust + Rglpk (exact solver, SOLVER=exact).
suppressPackageStartupMessages(library(igraph))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript copap_modules.R Coevolution.tsv annotation.tsv outdir [q_cutoff] [direction]")
coev_file <- args[1]; ann_file <- args[2]; outdir <- args[3]
cutoff <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
direction <- if (length(args) >= 5) args[5] else "co-occurrence"
edge_col <- Sys.getenv("EDGE_COL", "q-value")
solver <- Sys.getenv("SOLVER", "heuristic")
min_size <- as.integer(Sys.getenv("MIN_SIZE", "3"))   # smallest cluster reported as a module
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

coev <- read.delim(coev_file, check.names = FALSE, stringsAsFactors = FALSE)
ann <- read.delim(ann_file, check.names = FALSE, stringsAsFactors = FALSE)
names(ann)[1] <- "pos"
stopifnot(all(c("POS1", "POS2", "r", edge_col, "direction") %in% names(coev)))

# ---- 1. significant pairs = edges
keep <- coev[[edge_col]] < cutoff
if (direction != "both") keep <- keep & coev$direction == direction
sig <- coev[keep, ]
cat(sprintf("%d of %d pairs pass %s < %g (%s)\n", nrow(sig), nrow(coev), edge_col, cutoff, direction))
if (nrow(sig) == 0) stop("No significant pairs. Try a looser cutoff or EDGE_COL=p-value.")

nodes <- sort(unique(c(sig$POS1, sig$POS2)))
n <- length(nodes)

# ---- 2. similarity matrix for cluster editing
# Weight = |r| on significant pairs; non-significant pairs get a negative cost (-mean |r| of edges),
# so adding an edge that was NOT significant is penalised, deleting a significant one is penalised too.
S <- matrix(-mean(abs(sig$r)), n, n, dimnames = list(nodes, nodes))
for (k in seq_len(nrow(sig))) {
  i <- match(sig$POS1[k], nodes); j <- match(sig$POS2[k], nodes)
  S[i, j] <- S[j, i] <- abs(sig$r[k])
}
diag(S) <- 0

# ---- 3. weighted cluster editing
cost <- function(cl) {  # cost of making clusters perfect cliques
  same <- outer(cl, cl, "==")
  sum(abs(S)[upper.tri(S)] * ifelse(S[upper.tri(S)] > 0, !same[upper.tri(S)], same[upper.tri(S)]))
}
heuristic_cluster <- function(S, restarts = 50) {
  # start from connected components, then greedy single-node moves/merges until no improvement
  g <- graph_from_adjacency_matrix((S > 0) * 1, mode = "undirected", diag = FALSE)
  best <- NULL; bestc <- Inf
  for (r in seq_len(restarts)) {
    cl <- if (r == 1) as.integer(components(g)$membership) else sample(seq_len(nrow(S)), nrow(S), replace = TRUE)
    improved <- TRUE
    while (improved) {
      improved <- FALSE
      for (v in sample(seq_len(nrow(S)))) {
        labs <- c(unique(cl), max(cl) + 1L)       # existing clusters or a new singleton
        gain <- sapply(labs, function(L) sum(S[v, cl == L & seq_along(cl) != v]))
        cur <- sum(S[v, cl == cl[v] & seq_along(cl) != v])
        if (max(gain) > cur + 1e-12) { cl[v] <- labs[which.max(gain)]; improved <- TRUE }
      }
    }
    cl <- as.integer(factor(cl)); cc <- cost(cl)
    if (cc < bestc) { best <- cl; bestc <- cc }
  }
  best
}
set.seed(1)
if (solver == "exact" && requireNamespace("anticlust", quietly = TRUE)) {
  cl <- as.integer(anticlust::wce(S)$cluster)    # exact weighted cluster editing (needs Rglpk)
} else {
  if (solver == "exact") message("anticlust not installed: using the heuristic.")
  cl <- heuristic_cluster(S)
}
names(cl) <- nodes
sizes <- table(cl)
# stable naming: largest module first
ord <- names(sort(sizes, decreasing = TRUE)); ord <- ord[sizes[ord] >= min_size]
mod_name <- setNames(paste0("M", seq_along(ord)), ord)
module <- mod_name[as.character(cl)]; names(module) <- nodes

# ---- 4. outputs
res <- data.frame(pos = nodes, module = unname(module), stringsAsFactors = FALSE)
res$degree <- sapply(nodes, function(v) sum(S[as.character(v), ] > 0))
res <- merge(res, ann, by = "pos", all.x = TRUE)
res <- res[order(is.na(res$module), res$module, res$pos), ]
write.table(res, file.path(outdir, "modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# per-module summary
mods <- sort(unique(na.omit(res$module)))
summ <- do.call(rbind, lapply(mods, function(m) {
  v <- nodes[which(module == m)]
  ss <- S[as.character(v), as.character(v)]; pos <- ss[upper.tri(ss)]
  data.frame(module = m, n_families = length(v), families = paste(v, collapse = ","),
             edge_density = round(mean(pos > 0), 2), mean_abs_r = round(mean(pos[pos > 0]), 2))
}))
write.table(summ, file.path(outdir, "module_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nModules:\n"); print(summ, row.names = FALSE)
cat(sprintf("\nSingletons (significant but not in a module): %d\n", sum(is.na(module))))

# network plot
g <- graph_from_data_frame(data.frame(from = as.character(sig$POS1), to = as.character(sig$POS2), r = abs(sig$r)),
                           directed = FALSE, vertices = data.frame(name = as.character(nodes)))
pal <- c("#17695D", "#D9822B", "#2B5C9E", "#8E4A8A", "#B5A642", "#C0504D", "#4F9D69", "#5B6E91")
col <- ifelse(is.na(module), "#BBBBBB", pal[(match(module, mods) - 1) %% length(pal) + 1])
lab <- ann[match(nodes, ann$pos), 2]
png(file.path(outdir, "modules_network.png"), width = 1400, height = 1000, res = 130)
par(mar = c(1, 1, 2, 1))
set.seed(3)
plot(g, vertex.color = col, vertex.label = paste0(nodes), vertex.size = 14, vertex.label.color = "white",
     vertex.label.cex = 0.8, vertex.frame.color = NA, edge.width = 1 + 4 * E(g)$r, edge.color = "#55625D88",
     layout = layout_with_fr(g, weights = E(g)$r), main = sprintf("Co-evolving families (%s < %g): colour = module", edge_col, cutoff))
if (length(mods)) legend("bottomleft", legend = mods, fill = pal[(seq_along(mods) - 1) %% length(pal) + 1], bty = "n", cex = 0.9)
invisible(dev.off())
cat("\nWrote modules.tsv, module_summary.tsv, modules_network.png to", outdir, "\n")
