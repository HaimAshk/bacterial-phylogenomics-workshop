## Two conda environments

| Environment | Used for | Create |
|---|---|---|
| `phylo-genomes` | steps 1-4: genome download, skani, PopPUNK, dRep, plots | `conda env create -f envs/phylo-genomes.yml` |
| `phylo-gwas` | steps 5-6: GLOOME / CoPAP, treeWAS, pyseer | `conda env create -f envs/phylo-gwas.yml` then install treeWAS (below) |

```bash
conda activate phylo-gwas
Rscript -e 'devtools::install_github("caitiecollins/treeWAS", build_vignettes=FALSE)'
gloome --help
```

The tools are split because GLOOME pins `numpy<2`, which can conflict with the packages in the genome environment.
Not tested together on a clean machine yet: run the `conda env create` lines before the workshop.
