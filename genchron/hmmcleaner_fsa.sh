#!/bin/bash

set -euo pipefail

help="USAGE: hmmcleaner_fsa.sh <input cDNA fasta>"
# Modified from hmmcleaner_codons to reuse the protein alignment already existing.

inputal="${1:-}"

[[ -n "$1" ]] || { echo "$help" >&2 && exit 2; }


# Remove extension
if [[ "$inputal" = *.* ]]; then
    inputroot="${inputal%.*}"
else
    inputroot="$inputal"
fi
# Find the dirname and basename
if [[ "$inputroot" = */* ]]; then
    inputdir="${inputroot%/*}"
    inputbase="${inputroot##*/}"
else
    inputdir='.'
    inputbase="$inputroot"
fi

# Output directory is the parent folder of "realign/"
outputroot="${inputdir%/realign}/${inputbase}"


SCRIPTS="$HOME/scripts"

[[ -d "$SCRIPTS" ]] || { echo "Script directory does not exist." >&2 && exit 1; }

# Ungap and Convert to protein alignment
#$SCRIPTS/seqtools/ungap.py "$inputal" |\
#$SCRIPTS/seqtools/fasta_translate.py "${inputal}" "${inputbase}_prot.fa"

# Run HmmCleaner on the proteins
HmmCleaner.pl "${inputal}" > "${inputdir}/${inputbase}_cleanedcounts.txt"

# Find out which of the codon position correspond to the trimmed aa.
$SCRIPTS/seqtools/fillpositions.py -c -o "${outputroot%_protfsa}_fsa_hmm.fa" \
         "${inputdir}/${inputbase}_hmm.log" "${outputroot%_protfsa}_fsa.fa"

# Back-translate
#treebest backtrans "${inputal}_prot_hmm.fa" "${}" > ""
