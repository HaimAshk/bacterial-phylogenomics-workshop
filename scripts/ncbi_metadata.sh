#!/bin/bash
# Metadata (host, isolation source, pathovar/strain, country, year) for the 63 genomes via NCBI Datasets
# install: conda install -c conda-forge ncbi-datasets-cli   (provides `datasets` and `dataformat`)
ACC=accessions_63.txt

# 1) assembly-level report (organism, strain, isolate, BioSample accession, assembly level, CheckM...)
datasets summary genome accession --inputfile $ACC --as-json-lines > genomes.jsonl
dataformat tsv genome --inputfile genomes.jsonl \
  --fields accession,organism-name,organism-infraspecific-names-strain,organism-infraspecific-names-isolate,assminfo-biosample-accession,assminfo-level,checkm-completeness,assmstats-total-sequence-len,assminfo-biosample-host,assminfo-biosample-isolation-source,assminfo-biosample-geo-loc-name,assminfo-biosample-collection-date \
  > pse_assembly_metadata.tsv

# 2) if host/isolation_source columns are empty, pull full BioSample attributes (edirect: conda install -c bioconda entrez-direct)
cut -f5 pse_assembly_metadata.tsv | tail -n +2 | sort -u > biosamples.txt
while read b; do
  efetch -db biosample -id $b -format docsum -mode xml | \
  xtract -pattern DocumentSummary -element Accession -block Attribute -if Attribute@attribute_name -equals host -element Attribute \
         -block Attribute -if Attribute@attribute_name -equals isolation_source -element Attribute
done < biosamples.txt > pse_biosample_host.tsv
