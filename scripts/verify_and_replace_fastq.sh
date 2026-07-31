#!/usr/bin/env bash
# Verify a replacement FASTQ is a valid gzip, then install it into the local run copy
# and (optionally) the network share.
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <new_fastq.fq.gz> <dest_fastq.fq.gz> [also_copy_to_network_path]"
  echo "Example:"
  echo "  $0 ~/Downloads/SNU719_Rta-Zta_cntl1_1.fq.gz \\"
  echo "     /Volumes/TUNGSACore3/rnaseq_runs/SNU719_Rta-Zta_2025-04-10/cntl/SNU719_Rta-Zta_cntl1_1.fq.gz \\"
  echo "     /Volumes/FlemingtonLabMain1/2b_Flemington_Lab_Experiments/SNU719_Rta-Zta_2025-04-10/cntl/SNU719_Rta-Zta_cntl1_1.fq.gz"
  exit 1
fi
NEW=$1
DEST=$2
NET=${3:-}
echo "gzip -t $NEW ..."
gzip -t "$NEW"
echo "OK. Installing -> $DEST"
mkdir -p "$(dirname "$DEST")"
cp -f "$NEW" "$DEST"
gzip -t "$DEST"
echo "Installed and re-verified."
if [[ -n "$NET" ]]; then
  echo "Also installing -> $NET"
  mkdir -p "$(dirname "$NET")"
  cp -f "$NEW" "$NET"
  gzip -t "$NET"
fi
echo "Done."
