## Two conda environments

| Environment | Used for | Create |
|---|---|---|
| `phylo-genomes` | steps 1-4: genome download, skani, PopPUNK, dRep, plots | `conda env create -f envs/phylo-genomes.yml -p ./conda_envs/phylo-genomes` |
| `phylo-gwas` | steps 5-6: GLOOME / CoPAP, treeWAS, pyseer | `conda env create -f envs/phylo-gwas.yml -p ./conda_envs/phylo-gwas` then install treeWAS (below) |

```bash
# run from inside the repository folder (bacterial-phylogenomics-workshop/)
conda activate ./conda_envs/phylo-gwas
Rscript -e 'remotes::install_github("caitiecollins/treeWAS", build_vignettes=FALSE)'
gloome --help
```

Activate by path. Steps 1-4: `conda activate ./conda_envs/phylo-genomes`; steps 5-6: `conda activate ./conda_envs/phylo-gwas` (from inside `work/`, use `../conda_envs/...`).
`-p` puts the environments in a `conda_envs/` folder next to the repository files (useful when your home directory is small). Without `-p`, `conda env create -f envs/phylo-gwas.yml` creates a named environment instead, and `conda activate phylo-gwas` works.

The tools are split because GLOOME pins `numpy<2`, which can conflict with the packages in the genome environment.
Not tested together on a clean machine yet: run the `conda env create` lines before the workshop.
