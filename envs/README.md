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

If your home directory is small or you prefer a specific location, create the environments in a folder of your choice with `-p` (prefix) instead of `-n`/the name in the file, and activate them by path:

```bash
conda env create -f envs/phylo-genomes.yml -p /path/to/envs/phylo-genomes
conda env create -f envs/phylo-gwas.yml    -p /path/to/envs/phylo-gwas

conda activate /path/to/envs/phylo-genomes
conda activate /path/to/envs/phylo-gwas
```

The tools are split because GLOOME pins `numpy<2`, which can conflict with the packages in the genome environment.
Not tested together on a clean machine yet: run the `conda env create` lines before the workshop.
