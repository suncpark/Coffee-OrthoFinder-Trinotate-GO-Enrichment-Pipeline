# 🌱 Coffee OrthoFinder + Trinotate + GO Enrichment Pipeline

## Overview

This repository contains a reproducible comparative genomics pipeline for:

- Protein functional annotation using DIAMOND and HMMER (Pfam)
- Protein-only annotation using Trinotate
- OrthoFinder orthogroup functional summarization
- Identification of coffee-specific orthogroups and unassigned genes
- GO enrichment analysis using hypergeometric testing (R + BH FDR)

Originally developed for comparative analysis of Coffea species 
(C. arabica, C. canephora, C. eugenioides), the pipeline is fully generalizable to any species.

---

## Repository Structure

scripts/
│
├── 0_run_trinotate.sh  
├── 1a_run_family_annotation.sh  
├── 1b_annotate_orthogroups_from_hits.py  
├── 2_coffee_unique_sets.pl  
├── 3_GO_enrichment_assay.pl  

---

## Workflow Summary

### 1️⃣ Protein Family Annotation (DIAMOND + HMMER)

Script: scripts/1a_run_family_annotation.sh

Example:

```bash
bash scripts/1a_run_family_annotation.sh   --faa-dir ./all_cleaned_proteins   --combined combined.fasta   --dmnd uniprot_sprot.dmnd   --diamond-out all_species.outfmt6   --pfam-hmm Pfam-A.hmm   --hmm-out all_species.PFAMout   --cpus 40
```

Output:
- Combined FASTA
- DIAMOND outfmt6 file
- PFAM domtblout file
- HMMER log

---

### 2️⃣ Orthogroup Functional Annotation

Script: scripts/1b_annotate_orthogroups_from_hits.py

Example:

```bash
python scripts/1b_annotate_orthogroups_from_hits.py   --orthogroups-txt Orthogroups.txt   --diamond-outfmt6 all_species.outfmt6   --hmm-domtblout all_species.PFAMout   --gene-count-tsv Orthogroups.GeneCount.tsv
```

Output:
- Genefamily_UniProtAnnotation_DiamondBlastP.tab
- Genefamily_HmmAnnotation_Hmmsearch.tab
- Genefamily_UniProtAnnotation_CoffeaCount.tab (optional)

---

### 3️⃣ Coffee-Specific Orthogroups & Unassigned Genes

Script: scripts/2_coffee_unique_sets.pl

Example:

```bash
perl scripts/2_coffee_unique_sets.pl   --familywise 0_each_ortho_familywiseCount.tsv   --unassigned Orthogroups_UnassignedGenes.tsv   --outdir unique_sets
```

---

### 4️⃣ Trinotate Protein-Only Annotation

Script: scripts/0_run_trinotate.sh

Example:

```bash
bash scripts/0_run_trinotate.sh
```

Output:
- arabica.sqlite
- arabica_report.tab
- arabica_GO_trinotate.txt
- arabica_GO_trinotate_ancestral.txt

---

### 5️⃣ GO Enrichment Analysis

Script: scripts/3_GO_enrichment_assay.pl

Example:

```bash
perl scripts/3_GO_enrichment_assay.pl   --map mldbm_paths.tsv   --targets targets.tsv   --gene-dir ./gene_lists   --pvalue 0.01   --aspect P
```

Output:
- CSV file with enriched GO terms (includes p-value and FDR)

---

## Requirements

Core Tools:
- DIAMOND
- HMMER (hmmsearch)
- Trinotate
- OrthoFinder

Languages:
- Bash
- Python ≥ 3.8
- Perl ≥ 5.26

Perl Modules:
- Getopt::Long
- MLDBM
- DB_File
- Storable
- Statistics::R

R:
- Required for hypergeometric testing (phyper)

---

## Recommended Execution Order

1. 1a_run_family_annotation.sh
2. 1b_annotate_orthogroups_from_hits.py
3. 2_coffee_unique_sets.pl
4. 0_run_trinotate.sh
5. 3_GO_enrichment_assay.pl

---

## License

Add your preferred license (MIT recommended).
