#!/usr/bin/env python3
"""Thin out a phyletic-pattern FASTA (0/1 strings, one record per genome) before GLOOME / CoPAP.

CoPAP tests every pair of families, so run time and memory grow with the square of the number of
families. Families that are the same in (almost) all genomes cannot show co-evolution anyway.

Filters, applied in this order:
  1. invariant families (all 0 or all 1) are always dropped
  2. --min-minor N : keep a family only if the rarer state occurs in at least N genomes (default 3)
  3. --parsimony F : optional GLOOME ParsimonyAndHomoplasyScores.tsv from a first run on the same
                     pattern; with --min-steps K keep only families with parsimony score S >= K
                     (default 2: families that changed more than once on the tree)
  4. --dedupe      : identical columns are collapsed to one representative (members are listed in the map)
  5. --max-families M : if still more than M families, keep the M with the most balanced presence
                     (closest to 50% of genomes), a simple proxy for "most informative"

Outputs (prefix given by --out):
  PREFIX.fas          filtered phyletic pattern (same genome names, same order)
  PREFIX.map.tsv      new_pos (0-based, as GLOOME numbers them) -> orig_pos (0-based), n_present, merged_members
Join your annotation table to the map on orig_pos to annotate GLOOME / CoPAP output.

Example:
  python3 filter_phyletic_pattern.py pattern.fas --out filtered --min-minor 5 --dedupe --max-families 1000
Standard library only.
"""
import argparse
import csv
import sys


def read_fasta(path):
    names, seqs, cur = [], [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                names.append(line[1:])
                seqs.append([])
                cur = seqs[-1]
            else:
                cur.append(line)
    seqs = ["".join(s) for s in seqs]
    if not seqs:
        sys.exit("no sequences found")
    n = len(seqs[0])
    for nm, s in zip(names, seqs):
        if len(s) != n:
            sys.exit(f"sequence {nm} has length {len(s)}, expected {n}")
        if set(s) - set("01"):
            sys.exit(f"sequence {nm} contains characters other than 0 and 1")
    return names, seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-minor", type=int, default=3)
    ap.add_argument("--parsimony", help="GLOOME ParsimonyAndHomoplasyScores.tsv from a first run")
    ap.add_argument("--min-steps", type=int, default=2)
    ap.add_argument("--dedupe", action="store_true")
    ap.add_argument("--max-families", type=int)
    a = ap.parse_args()

    names, seqs = read_fasta(a.fasta)
    g, n = len(seqs), len(seqs[0])
    cols = ["".join(s[i] for s in seqs) for i in range(n)]
    present = [c.count("1") for c in cols]
    keep = [i for i in range(n) if 0 < present[i] < g]
    print(f"{n} families, {g} genomes; {n - len(keep)} invariant dropped", file=sys.stderr)

    keep = [i for i in keep if min(present[i], g - present[i]) >= a.min_minor]
    print(f"{len(keep)} left after --min-minor {a.min_minor}", file=sys.stderr)

    if a.parsimony:
        steps = {}
        with open(a.parsimony) as fh:
            rd = csv.reader(fh, delimiter="\t")
            head = next(rd)
            si = [k for k, h in enumerate(head) if h.startswith("Parsimony Score")][0]
            for row in rd:
                if row[0] != "TOTAL":
                    steps[int(row[0])] = int(float(row[si]))
        keep = [i for i in keep if steps.get(i, 0) >= a.min_steps]
        print(f"{len(keep)} left after parsimony score >= {a.min_steps}", file=sys.stderr)

    members = {i: [] for i in keep}
    if a.dedupe:
        seen, uniq = {}, []
        for i in keep:
            if cols[i] in seen:
                members[seen[cols[i]]].append(i)
            else:
                seen[cols[i]] = i
                uniq.append(i)
        keep = uniq
        print(f"{len(keep)} left after collapsing identical patterns", file=sys.stderr)

    if a.max_families and len(keep) > a.max_families:
        keep = sorted(sorted(keep, key=lambda i: abs(present[i] / g - 0.5))[: a.max_families])
        print(f"{len(keep)} kept (most balanced presence)", file=sys.stderr)

    if not keep:
        sys.exit("no families left: loosen the filters")

    with open(a.out + ".fas", "w") as fh:
        for nm, s in zip(names, seqs):
            fh.write(f">{nm}\n{''.join(s[i] for i in keep)}\n")
    with open(a.out + ".map.tsv", "w") as fh:
        fh.write("new_pos\torig_pos\tn_present\tmerged_members\n")
        for j, i in enumerate(keep):
            fh.write(f"{j}\t{i}\t{present[i]}\t{','.join(map(str, members.get(i, [])))}\n")
    print(f"wrote {a.out}.fas ({len(keep)} families) and {a.out}.map.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
