# Real treeWAS (Collins & Didelot 2018, PLoS Comput Biol) on the toy example
# install.packages("remotes"); remotes::install_git("https://github.com/caitiecollins/treeWAS.git", build_vignettes=FALSE)
library(treeWAS); library(ape)
tree <- read.tree("toy2_tree.nwk")
pa   <- read.delim("toy2_PA_matrix.tsv", row.names=1, check.names=FALSE)
ph   <- read.delim("toy2_trait.tsv", row.names=1)
pa   <- as.matrix(pa[tree$tip.label,])
y    <- setNames(ph[tree$tip.label,"phenotype"], tree$tip.label)
out  <- treeWAS(snps = pa, phen = y, tree = tree,
                n.subs = NULL, seed = 1, plot.tree = FALSE,
                snps.reconstruction = "parsimony", phen.reconstruction = "parsimony",
                p.value = 0.01, p.value.correct = "bonf", p.value.by = "count")
print(out$terminal$sig.snps); print(out$simultaneous$sig.snps); print(out$subsequent$sig.snps)
# truth: M1 genes are in toy2_truth.tsv (role == module_M1); decoys = decoy_cladeA / decoy_cladesCD
