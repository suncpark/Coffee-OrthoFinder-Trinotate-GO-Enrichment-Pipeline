#!/usr/bin/env bash
set -euo pipefail

# run_family_annotation.sh
# ------------------------------------------------------------
# 1) Combine protein FASTAs in a directory
# 2) Run DIAMOND blastp against a UniProt .dmnd database
# 3) Run hmmsearch against Pfam-A.hmm (or any HMM database)
#
# Example:
#   bash run_family_annotation.sh \
#     --faa-dir ./all_cleaned_proteins_coffea \
#     --combined _all_coffea_genomes.fas \
#     --dmnd /path/to/uniprot_sprot.dmnd \
#     --diamond-out all_coffea_diamond.outfmt6 \
#     --pfam-hmm /path/to/Pfam-A.hmm \
#     --hmm-out all_coffea.PFAMout \
#     --cpus 40
#
# Notes:
# - The DIAMOND output is outfmt6:
#   qseqid sseqid evalue bitscore length
# - hmmsearch output is --domtblout format
# ------------------------------------------------------------

usage() {
  cat <<'EOF'
Usage:
  run_family_annotation.sh --faa-dir DIR --combined FILE --dmnd FILE --diamond-out FILE --pfam-hmm FILE --hmm-out FILE [--cpus N]

Required:
  --faa-dir        Directory containing *.fa/*.faa/*.fasta/*.fas protein FASTA files
  --combined       Output combined FASTA filename
  --dmnd           DIAMOND database (.dmnd), e.g., uniprot_sprot.dmnd
  --diamond-out    DIAMOND output file (outfmt6)
  --pfam-hmm       HMM database, e.g., Pfam-A.hmm
  --hmm-out        hmmsearch domtblout output filename

Optional:
  --cpus           Threads (default: 20)
EOF
}

CPUS=20
FAA_DIR=""
COMBINED=""
DMND=""
DIAMOND_OUT=""
PFAM_HMM=""
HMM_OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --faa-dir) FAA_DIR="$2"; shift 2;;
    --combined) COMBINED="$2"; shift 2;;
    --dmnd) DMND="$2"; shift 2;;
    --diamond-out) DIAMOND_OUT="$2"; shift 2;;
    --pfam-hmm) PFAM_HMM="$2"; shift 2;;
    --hmm-out) HMM_OUT="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1;;
  esac
done

[[ -z "$FAA_DIR" || -z "$COMBINED" || -z "$DMND" || -z "$DIAMOND_OUT" || -z "$PFAM_HMM" || -z "$HMM_OUT" ]] && { usage; exit 1; }
[[ ! -d "$FAA_DIR" ]] && { echo "ERROR: --faa-dir not found: $FAA_DIR" >&2; exit 1; }
[[ ! -f "$DMND" ]] && { echo "ERROR: --dmnd not found: $DMND" >&2; exit 1; }
[[ ! -f "$PFAM_HMM" ]] && { echo "ERROR: --pfam-hmm not found: $PFAM_HMM" >&2; exit 1; }

echo "[INFO] CPUs = $CPUS"
echo "[INFO] Combining FASTA files in: $FAA_DIR"
echo "[INFO] Output combined FASTA: $COMBINED"

# Combine FASTAs robustly:
# - Writes each record as:
#   >ID
#   SEQUENCE
python3 - <<'PY' "$FAA_DIR" "$COMBINED"
import sys, os, glob

faa_dir = sys.argv[1]
out_fa  = sys.argv[2]

patterns = ["*.fa", "*.faa", "*.fasta", "*.fas"]
files = []
for p in patterns:
    files.extend(glob.glob(os.path.join(faa_dir, p)))
files = sorted(set(files))

if not files:
    raise SystemExit(f"ERROR: No FASTA files found in {faa_dir} with patterns {patterns}")

def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]  # keep only first token as ID
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)

with open(out_fa, "w", encoding="utf-8") as out:
    for fp in files:
        for hid, seq in fasta_iter(fp):
            out.write(f">{hid}\n{seq}\n")

print(f"[INFO] Combined {len(files)} files into {out_fa}", file=sys.stderr)
PY

echo "[INFO] Running DIAMOND blastp..."
echo "  diamond blastp --query $COMBINED --db $DMND --outfmt 6 qseqid sseqid evalue bitscore length --evalue 1e-5 --max-hsps 1 --threads $CPUS --out $DIAMOND_OUT"
diamond blastp \
  --query "$COMBINED" \
  --db "$DMND" \
  --outfmt 6 qseqid sseqid evalue bitscore length \
  --evalue 1e-5 \
  --max-hsps 1 \
  --threads "$CPUS" \
  --out "$DIAMOND_OUT"

echo "[INFO] Running hmmsearch..."
echo "  hmmsearch --cpu $CPUS --noali --domtblout $HMM_OUT $PFAM_HMM $COMBINED"
hmmsearch \
  --cpu "$CPUS" \
  --noali \
  --domtblout "$HMM_OUT" \
  "$PFAM_HMM" \
  "$COMBINED" \
  > "${HMM_OUT}.log"

echo "[DONE] Outputs:"
echo "  Combined FASTA : $COMBINED"
echo "  DIAMOND outfmt6: $DIAMOND_OUT"
echo "  HMM domtblout  : $HMM_OUT"
echo "  HMM log        : ${HMM_OUT}.log"
