#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

help="
For each input fasta (.bz2) file, output one line:

<file name><tab><list of sequence ids>
"


sed -n 'v' <<< ""  || { echo "Error: GNU sed not installed.">&2 && exit 1; }

# One-liner for already decompressed files
#sed -snr '/^>/{s/^>(\S+).*$/\1/;H};${x;F;s/\n/ /gp}' $@ | sed 'N;s/\n /\t/'

rebzip() {
    bzip2 "$fastafile"
}

trap 'rebzip' TERM INT ERR


while (( $# )); do
    fastafile="${1}"
    shift

    bz2=0
    if [[ "$fastafile" =~ \.bz2$ ]]; then
        bunzip2 "$fastafile"
        fastafile=${fastafile%.bz2}
        bz2=1
    fi
    echo -ne "${fastafile}\t"
    sed -nr '/^>/{s/^>(\S+).*$/\1/;1h;1!H};${x;s/\n/ /gp}' "$fastafile"
    (( bz2 )) && bzip2 "$fastafile" &
done
