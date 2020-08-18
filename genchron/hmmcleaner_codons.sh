#!/bin/bash

set -euo pipefail

help="USAGE: hmmcleaner_codons.sh <input cDNA fasta>"

inputal="${1:-}"

[[ -n "$1" ]] || { echo "$help" >&2 && exit 2; }


# Remove extension
if [[ "$inputal" = *.* ]]; then
    inputbase="${1%.*}"
else
    inputbase="$inputal"
fi


SCRIPTS="$HOME/scripts"

[[ -d "$SCRIPTS" ]] || { echo "Script directory does not exist." >&2 && exit 1; }

# Ungap and Convert to protein alignment
#$SCRIPTS/seqtools/ungap.py "$inputal" |\
$SCRIPTS/seqtools/fasta_translate.py "${inputal}" "${inputbase}_prot.fa"

# Run HmmCleaner on the proteins
HmmCleaner.pl "${inputbase}_prot.fa" >"${inputbase}_prot_cleanedcounts.txt"

# Find out which of the codon position correspond to the trimmed aa.
$SCRIPTS/seqtools/fillpositions.py -c -o "${inputbase}_hmm.fa" \
         "${inputbase}_prot_hmm.log" "${inputal}"

# Back-translate
#treebest backtrans "${inputal}_prot_hmm.fa" "${}" > ""
