#!/usr/bin/env bash
# Build a Microbializer input zip (one FASTA per genome) from a list of NCBI accessions.
# Usage: bash make_microbializer_zip.sh accessions.txt [output_name]
# Requires: ncbi-datasets-cli (conda install -c conda-forge ncbi-datasets-cli), unzip, zip
set -euo pipefail

ACC_FILE="${1:?Usage: bash make_microbializer_zip.sh accessions.txt [output_name]}"
OUT="${2:-microbializer_input}"

for tool in datasets unzip zip; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done

WORK="${OUT}_work"
rm -rf "$WORK" "$OUT" "$OUT.zip"
mkdir -p "$WORK" "$OUT"

# Drop blank lines and stray whitespace
grep -v '^[[:space:]]*$' "$ACC_FILE" | tr -d '\r' | sed 's/[[:space:]]*$//' > "$WORK/accessions_clean.txt"
N_EXPECTED=$(wc -l < "$WORK/accessions_clean.txt")
echo "Downloading $N_EXPECTED genomes..."

datasets download genome accession --inputfile "$WORK/accessions_clean.txt" \
  --include genome --filename "$WORK/genomes.zip"
unzip -q "$WORK/genomes.zip" -d "$WORK/unzipped"

# One FASTA per genome. All sequences of an assembly (chromosome + plasmids) go into one file.
# Dots in accessions become underscores (GCF_000007805.1 -> GCF_000007805_1.fna) so the
# genome name has no extra dot.
for d in "$WORK"/unzipped/ncbi_dataset/data/GC*; do
  [ -d "$d" ] || continue
  acc=$(basename "$d")
  safe=${acc//./_}
  cat "$d"/*.fna > "$OUT/${safe}.fna"
done

N_FOUND=$(ls "$OUT"/*.fna | wc -l)
echo "Wrote $N_FOUND FASTA files to $OUT/ (expected $N_EXPECTED)"
if [ "$N_FOUND" -ne "$N_EXPECTED" ]; then
  echo "WARNING: count mismatch. Accessions with no file:" >&2
  while read -r acc; do
    [ -f "$OUT/${acc//./_}.fna" ] || echo "  $acc" >&2
  done < "$WORK/accessions_clean.txt"
fi

# Sanity check: every file should start with a FASTA header and not be tiny
for f in "$OUT"/*.fna; do
  head -c1 "$f" | grep -q '>' || echo "WARNING: $f does not start with '>'" >&2
  [ "$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f")" -gt 1000000 ] || echo "WARNING: $f is under 1 MB" >&2
done

zip -qr "$OUT.zip" "$OUT"
echo "Created $OUT.zip ($(du -h "$OUT.zip" | cut -f1)). Upload this to https://microbializer.tau.ac.il/"
echo "Intermediate files are in $WORK/ and can be deleted."
