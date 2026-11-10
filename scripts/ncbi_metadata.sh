#!/bin/bash
# Metadata (host, isolation source, strain, country, date, assembly quality) for a list of genomes via NCBI Datasets
# usage:  bash ncbi_metadata.sh accessions_63.txt        # writes metadata.tsv
#         bash ncbi_metadata.sh accessions_63.txt --biosample   # optional: also fetch full BioSample host / isolation_source (needs edirect and access to NCBI E-utilities)
# install: conda install -c conda-forge ncbi-datasets-cli   (provides `datasets` and `dataformat`)
ACC=${1:-accessions_63.txt}

# 1) assembly-level report: organism, strain, isolate, BioSample accession, assembly level, CheckM, length, host, isolation source, country, date
datasets summary genome accession --inputfile "$ACC" --as-json-lines > genomes.jsonl
dataformat tsv genome --inputfile genomes.jsonl \
  --fields accession,organism-name,organism-infraspecific-strain,organism-infraspecific-isolate,assminfo-biosample-accession,assminfo-level,checkm-completeness,assmstats-total-sequence-len,assminfo-biosample-host,assminfo-biosample-isolation-source,assminfo-biosample-geo-loc-name,assminfo-biosample-collection-date \
  > metadata.tsv
echo "wrote metadata.tsv ($(($(wc -l < metadata.tsv) - 1)) genomes)"

# 2) optional: if host / isolation_source columns are empty, pull the full BioSample attributes
#    (edirect: conda install -c bioconda entrez-direct; some networks block E-utilities, then skip this step)
if [ "$2" = "--biosample" ]; then
  cut -f5 metadata.tsv | tail -n +2 | sort -u > biosamples.txt
  while read b; do
    efetch -db biosample -id $b -format docsum -mode xml | \
    xtract -pattern DocumentSummary -element Accession -block Attribute -if Attribute@attribute_name -equals host -element Attribute \
           -block Attribute -if Attribute@attribute_name -equals isolation_source -element Attribute
  done < biosamples.txt > biosample_host.tsv
fi
