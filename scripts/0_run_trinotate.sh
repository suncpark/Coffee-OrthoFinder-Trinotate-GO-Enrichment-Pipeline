#!/usr/bin/env bash
# ============================================================
# Trinotate (protein-only) workflow
# ============================================================
# Purpose:
#   Build a Trinotate SQLite DB for *Coffea arabica* proteins,
#   load DIAMOND BLASTP + PFAM hmmsearch results, generate a report,
#   and extract GO terms (with/without ancestral terms).
#
# Assumptions:
#   - You already ran DIAMOND blastp and hmmsearch on Ubuntu (or elsewhere)
#     and have:
#       1) C.arabica_diamondblastp_uniprot.outfmt6
#       2) arabica.PFAMout   (hmmsearch --domtblout output)
#   - You have a working Trinotate installation in conda env "TN"
#   - TRINOTATE_DATA_DIR points to Trinotate resource DB directory
#
# Run:
#   bash trinotate_arabica_onepage.sh
#
# ============================================================
set -euo pipefail

# ----------------------------
# Step 0. Load environment
# ----------------------------
# Activate Trinotate conda environment (must contain Trinotate + helper scripts)
conda activate TN

# Tell Trinotate where its resource databases live (SwissProt, Pfam mapping, etc.)
export TRINOTATE_DATA_DIR=/mnt/d/ps3_GenomeWideGeneFamilyAnalysis/lettuce_trinotate/db

# Move to a working directory where you want the sqlite and outputs to live
cd /mnt/d/plantGenomes/Genitianales_coffea/all_cleaned_proteins_coffea/OrthoFinder/0_trinotate

# ----------------------------
# Step 1. Define inputs/outputs
# ----------------------------
# Protein FASTA (your Arabica protein set)
PEP_FASTA="/mnt/d/plantGenomes/Genitianales_coffea/all_cleaned_proteins_coffea/CoffeaArabica_LP_NoPartial.fas"

# DIAMOND blastp outfmt6 against UniProt/SwissProt (qseqid sseqid evalue bitscore length ...)
BLASTP_OUT="/mnt/d/plantGenomes/Genitianales_coffea/all_cleaned_proteins_coffea/OrthoFinder/1_annotate_family/C.arabica_diamondblastp_uniprot.outfmt6"

# PFAM hmmsearch output (recommended: hmmsearch --domtblout)
PFAM_OUT="/mnt/d/plantGenomes/Genitianales_coffea/all_cleaned_proteins_coffea/OrthoFinder/1_annotate_family/arabica.PFAMout"

# Trinotate SQLite DB to create
SQLITE_DB="arabica.sqlite"

# Symlink name (Trinotate likes a local input file; link keeps path simple)
PEP_LINK="arabica.pep"

# Report and GO outputs
REPORT_TAB="arabica_report.tab"
GO_TXT="arabica_GO_trinotate.txt"
GO_ANC_TXT="arabica_GO_trinotate_ancestral.txt"

# ----------------------------
# Step 2. Create Trinotate SQLite DB
# ----------------------------
# Creates an empty Trinotate database schema in arabica.sqlite
Trinotate --db "${SQLITE_DB}" --create

# ----------------------------
# Step 3. Load protein sequences into the DB
# ----------------------------
# Create/update symlink to your protein FASTA
ln -sf "${PEP_FASTA}" "${PEP_LINK}"

# Initialize DB using protein sequences only (no transcript/CDS needed)
# This loads the protein IDs into the Trinotate schema.
Trinotate --db "${SQLITE_DB}" --init --transdecoder_pep "${PEP_LINK}"

# If the above is slow or you prefer explicit loader:
# TrinotateSeqLoader.pl --sqlite "${SQLITE_DB}" --transdecoder_pep "${PEP_LINK}" --bulk_load

# ----------------------------
# Step 4. Load DIAMOND BLASTP results
# ----------------------------
# Loads your blastp hits into the sqlite.
Trinotate --db "${SQLITE_DB}" --LOAD_swissprot_blastp "${BLASTP_OUT}"

# If the above command fails for any reason, use the direct loader:
# Trinotate_BLAST_loader.pl --sqlite "${SQLITE_DB}" --outfmt6 "${BLASTP_OUT}"

# ----------------------------
# Step 5. Load PFAM hmmsearch results
# ----------------------------
# Loads PFAM domain hits into the sqlite.
Trinotate --db "${SQLITE_DB}" --LOAD_pfam "${PFAM_OUT}"

# If the above fails, use the direct loader:
# Trinotate_PFAM_loader.pl --sqlite "${SQLITE_DB}" --pfam "${PFAM_OUT}"

# ----------------------------
# Step 6. Generate Trinotate report
# ----------------------------
# Produces the main report table (can take a long time depending on DB size).
Trinotate --db "${SQLITE_DB}" --report > "${REPORT_TAB}"

# ----------------------------
# Step 7. Extract GO assignments
# ----------------------------
# Extract GO terms (gene-level) WITHOUT ancestral/parent GO terms
extract_GO_assignments_from_Trinotate_xls.pl \
  --Trinotate_xls "${REPORT_TAB}" \
  -G > "${GO_TXT}"

# Extract GO terms (gene-level) WITH ancestral/parent GO terms
extract_GO_assignments_from_Trinotate_xls.pl \
  --Trinotate_xls "${REPORT_TAB}" \
  -G -I > "${GO_ANC_TXT}"


echo "[DONE]"
echo "SQLite DB   : ${SQLITE_DB}"
echo "Report      : ${REPORT_TAB}"
echo "GO          : ${GO_TXT}"
echo "GO ancestral: ${GO_ANC_TXT}"
