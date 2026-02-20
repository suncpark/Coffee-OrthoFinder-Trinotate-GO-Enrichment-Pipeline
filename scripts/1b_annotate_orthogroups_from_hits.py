#!/usr/bin/env python3
"""
annotate_orthogroups_from_hits.py
------------------------------------------------------------
Parses:
  1) OrthoFinder Orthogroups.txt (OG0000000: gene1 gene2 ...)
  2) DIAMOND blastp outfmt6 (qseqid sseqid evalue bitscore length)
  3) hmmsearch --domtblout output (Pfam or other HMM db)

Produces:
  A) Genefamily_UniProtAnnotation_DiamondBlastP.tab
     famID  size  HitFreq  UniProHit  Description
     - UniProHit chosen as most frequent sseqid among family members,
       counting only hits with evalue <= --blast-cutoff (default 1e-10).

  B) Genefamily_HmmAnnotation_Hmmsearch.tab
     famID  fam.size  HmmHit
     - HmmHit is a '^' joined list of PFAM target_name;accession pairs
       with (freq >= ceil(size * --hmm-min-frac)) and evalue <= --hmm-cutoff

  C) Genefamily_UniProtAnnotation_CoffeaCount.tab  (optional)
     famID size UniProHit hmm_hit rob eug arabica
     - requires --gene-count-tsv with OrthoFinder Orthogroups.GeneCount.tsv
       and column names for coffee species (defaults provided).

Optional:
  --uniprot-desc-tsv : tab-separated mapping {id <tab> description}
     - used to populate Description in output A
     - if not provided, Description="NA"
------------------------------------------------------------
"""

from __future__ import annotations
import argparse
import math
import sys
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional


def read_orthogroups_txt(path: Path) -> Dict[str, List[str]]:
    og = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or ":" not in line:
                continue
            fam, rest = line.split(":", 1)
            fam = fam.strip()
            genes = rest.strip().split()
            og[fam] = genes
    return og


def read_uniprot_desc_tsv(path: Path) -> Dict[str, str]:
    d = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            uid = parts[0].strip()
            desc = parts[1].strip()
            d[uid] = desc
    return d


def parse_diamond_outfmt6(path: Path, evalue_cutoff: float) -> Dict[str, List[Tuple[str, float]]]:
    """
    Returns:
      hits[qseqid] = [(sseqid, evalue), ...]  (only rows <= evalue_cutoff)
    """
    hits = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            q, s, ev = parts[0], parts[1], parts[2]
            try:
                e = float(ev)
            except ValueError:
                continue
            if e <= evalue_cutoff:
                hits[q].append((s, e))
    return hits


def parse_hmm_domtblout(path: Path, evalue_cutoff: float) -> Dict[str, List[Tuple[str, str, float]]]:
    """
    hmmsearch --domtblout format:
      target name, target acc, query name, ... , i-Evalue (column 13) ...
    We'll parse minimally:
      - target_name (parts[0])
      - target_acc  (parts[1])
      - query_name  (parts[3])
      - i_evalue    (parts[12])  (1-based: 13th column)

    Returns:
      hits[query_name] = [(target_name, target_acc, i_evalue), ...]
    """
    hits = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for ln, line in enumerate(f, 1):
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 13:
                continue
            target_name = parts[0]
            target_acc = parts[1]
            query_name = parts[3]
            try:
                i_eval = float(parts[12])
            except ValueError:
                continue
            if i_eval <= evalue_cutoff:
                hits[query_name].append((target_name, target_acc, i_eval))
    return hits


def read_gene_count_tsv(path: Path) -> Tuple[List[str], Dict[str, Dict[str, int]]]:
    """
    OrthoFinder Orthogroups.GeneCount.tsv is typically:
      Orthogroup <tab> Species1 <tab> Species2 ...
    Returns:
      (header_species, counts[orthogroup][species] = int)
    """
    counts = {}
    header_species: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 2:
            raise ValueError("GeneCount TSV header looks too short.")
        # header[0] is 'Orthogroup' (or similar)
        header_species = header[1:]
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            og = parts[0]
            row = {}
            for sp, val in zip(header_species, parts[1:]):
                try:
                    row[sp] = int(val)
                except ValueError:
                    row[sp] = 0
            counts[og] = row
    return header_species, counts


def most_frequent_uniprot_hit(
    fam_genes: List[str],
    blast_hits: Dict[str, List[Tuple[str, float]]],
) -> Tuple[str, str]:
    """
    Returns (hit_id, freq_as_str)
    If none found => ("NA","NA")
    Frequency is count of occurrences across all qualifying hits for all genes in family.
    """
    c = Counter()
    for g in fam_genes:
        for (hit, _e) in blast_hits.get(g, []):
            c[hit] += 1
    if not c:
        return "NA", "NA"
    hit, freq = c.most_common(1)[0]
    return hit, str(freq)


def hmm_hits_above_fraction(
    fam_genes: List[str],
    hmm_hits: Dict[str, List[Tuple[str, str, float]]],
    min_frac: float,
) -> str:
    """
    Returns '^'-joined list of "name;acc" for hits appearing in >= ceil(size*min_frac) members
    (counting all qualifying hits across genes).
    If none => "NA"
    """
    size = len(fam_genes)
    if size == 0:
        return "NA"
    threshold = max(1, math.ceil(size * min_frac))
    c = Counter()
    for g in fam_genes:
        for (name, acc, _e) in hmm_hits.get(g, []):
            c[f"{name};{acc}"] += 1
    if not c:
        return "NA"
    kept = [k for k, v in c.most_common() if v >= threshold]
    return "^".join(kept) if kept else "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthogroups-txt", required=True, type=Path,
                    help="OrthoFinder Orthogroups.txt")
    ap.add_argument("--diamond-outfmt6", required=True, type=Path,
                    help="DIAMOND blastp outfmt6 (qseqid sseqid evalue bitscore length)")
    ap.add_argument("--hmm-domtblout", required=True, type=Path,
                    help="hmmsearch --domtblout output")
    ap.add_argument("--blast-cutoff", type=float, default=1e-10,
                    help="E-value cutoff for counting DIAMOND hits (default: 1e-10)")
    ap.add_argument("--hmm-cutoff", type=float, default=1e-5,
                    help="E-value cutoff for counting HMM hits (default: 1e-5)")
    ap.add_argument("--hmm-min-frac", type=float, default=0.05,
                    help="Minimum fraction of family members supporting an HMM hit (default: 0.05 = 5%%)")
    ap.add_argument("--uniprot-desc-tsv", type=Path, default=None,
                    help="Optional TSV: UniProtID<TAB>Description (used for Description column)")
    ap.add_argument("--gene-count-tsv", type=Path, default=None,
                    help="Optional OrthoFinder Orthogroups.GeneCount.tsv to add coffee counts")
    ap.add_argument("--robusta-col", default="CoffeaRobusta",
                    help="Species column name for robusta in GeneCount TSV")
    ap.add_argument("--eugenioides-col", default="CoffeaEugenioides",
                    help="Species column name for eugenioides in GeneCount TSV")
    ap.add_argument("--arabica-col", default="CoffeaArabica",
                    help="Species column name for arabica in GeneCount TSV")

    args = ap.parse_args()

    og = read_orthogroups_txt(args.orthogroups_txt)
    blast = parse_diamond_outfmt6(args.diamond_outfmt6, args.blast_cutoff)
    hmm = parse_hmm_domtblout(args.hmm_domtblout, args.hmm_cutoff)

    uniprot_desc = {}
    if args.uniprot_desc_tsv:
        uniprot_desc = read_uniprot_desc_tsv(args.uniprot_desc_tsv)

    # Output A: UniProt top hit per family
    outA = Path("Genefamily_UniProtAnnotation_DiamondBlastP.tab")
    with outA.open("w", encoding="utf-8") as w:
        w.write("\t".join(["famID", "size", "HitFreq", "UniProHit", "Description"]) + "\n")
        for fam in sorted(og.keys()):
            genes = og[fam]
            hit, freq = most_frequent_uniprot_hit(genes, blast)
            desc = uniprot_desc.get(hit, "NA")
            desc = desc.replace(",", " ")  # mimic your older cleanup
            w.write("\t".join([fam, str(len(genes)), freq, hit, desc]) + "\n")

    # Output B: HMM hits >=5% frequency (default)
    outB = Path("Genefamily_HmmAnnotation_Hmmsearch.tab")
    with outB.open("w", encoding="utf-8") as w:
        w.write("\t".join(["famID", "fam.size", "HmmHit"]) + "\n")
        for fam in sorted(og.keys()):
            genes = og[fam]
            hits = hmm_hits_above_fraction(genes, hmm, args.hmm_min_frac)
            w.write("\t".join([fam, str(len(genes)), hits]) + "\n")

    # Output C: Coffee count summary (optional)
    if args.gene_count_tsv:
        header_species, counts = read_gene_count_tsv(args.gene_count_tsv)
        outC = Path("Genefamily_UniProtAnnotation_CoffeaCount.tab")
        with outC.open("w", encoding="utf-8") as w:
            w.write("\t".join(["famID", "size", "UniProHit", "hmm_hit", "rob", "eug", "arabica"]) + "\n")
            for fam in sorted(og.keys()):
                genes = og[fam]
                hit, _freq = most_frequent_uniprot_hit(genes, blast)
                hmmhit = hmm_hits_above_fraction(genes, hmm, args.hmm_min_frac)
                row = counts.get(fam, {})
                rob = row.get(args.robusta_col, 0)
                eug = row.get(args.eugenioides_col, 0)
                ara = row.get(args.arabica_col, 0)
                w.write("\t".join([fam, str(len(genes)), hit, hmmhit, str(rob), str(eug), str(ara)]) + "\n")

    print("[DONE] Wrote:", outA, outB, ("Genefamily_UniProtAnnotation_CoffeaCount.tab" if args.gene_count_tsv else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
