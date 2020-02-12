#!/bin/bash
set -euo pipefail

[[ "$2" =~ [0-9]+ ]] || { echo "USAGE: ./generax_family1by1.sh <familyfile> <family number (integer)>">&2 ; exit 2 ; }

report_mem_usage() {
    echo 'Current Memory usage:'
    # Apparently, ps is not the right tool reporting memory usage.
    ps -axh -q "$$" -o 'vsz'
}

maxmem="${3:-}"  # In MB
if [[ -n "$maxmem" ]]; then
    # Kill the process if exceeds given virtual memory
    ulimit -a
    memKB=$((maxmem*1024))
    echo "Setting Hard memory limit: ${maxmem} MB ($memKB)."
    ulimit -v $memKB

    trap report_mem_usage ERR SIGINT SIGTERM SIGKILL
fi

familychunk=$(/bin/sed -n '1p;'"$(($2 * 4 + 2))"',+3p' "$1")
echo -e "--- Select family ${2}:\n${familychunk}\n---"

printf -v nb '%05d' "$2"
mkdir -p "generax_family-1by1/${nb::2}"

#spfile="generax_family-1by1/${2::2}/${2:2}/labelled_species_tree.newick"
## Check that $spfile is not currently being written to.
## -n makes lsof faster.
## -F will output the fields:
## f: file descriptor (expect an integer);
## a: access mode (w for write, u for read and write);
## l: lock type: space for none, r and R for reading, w|W|u|U and also N|x|X
#if [ -f "$spfile" ]; then
#while lsof -n -Ffal $spfile | egrep -xq 'a[wWuU]|l[wWuUNxX]'; do
#    echo "File ${spfile} is busy (written too); retrying in 2 sec."
#    sleep 2;
#done;
#fi
echo '--- Status of the output directory (lsof):'
lsof -n "generax_family-1by1/${nb::2}" || echo 'Unused.'
echo '---'

generax -f <(echo "${familychunk}") -s ../PhylTree.TimeTree201901.Ensembl93-like.goodQual.binary.shortnames.nwk -p "generax_family-1by1/${nb::2}/${nb:2}"

