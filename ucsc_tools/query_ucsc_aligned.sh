#!/bin/bash

set -euo pipefail

help="USAGE: ./query_ucsc_aligned.sh <db> <table> <position>

<table> ideally a 'net' table
"

if [ $# -ne 3 ]; then
	echo "$help"
	exit 1
else
	db=$1
	table=$2
	position=$3
fi

#OPTIND=1         # Reset in case getopts has been used previously in the shell.
#
## Initialize our own variables:
#while getopts "h?vbq" opt; do
#	case "$opt" in
#	h|\?)
#		echo "$help"
#		exit 0
#		;;
#	v)  valgrind=1
#		;;
#	b)  buffer=1
#		;;
#	q)  quiet=1
#		;;
#	esac
#done
#
#shift $((OPTIND-1))
#
#[ "$1" = "--" ] && shift
#
#if [ -z "$1" ]; then
#	echo "$help"
#	exit 1
#fi

#db='mm10'
##table='chainProCap1'
#table='netProCap1'
#position='chr17:71,344,493-71,475,343'
position="${position//,/}"
tname="${position%%:*}"
range="${position#*:}"
tstart="${range%%-*}"
tend="${range#*-}"

#columns="level,
#         tName,
#         tStart,
#         tEnd,
#         strand,
#         qName,
#         qStart,
#         qEnd,
#         qSize,
#         ali,
#         type,
#         score"

#columns="tName, tStart, tEnd, qStrand, qName, qStart, qEnd, qSize, score"
columns="level,type,tName,tStart,tEnd,strand,qName,qStart,qEnd,chainId,score"


# Level of alignment
# Target chromosome
# Start on target
# End on target
# Orientation of query. + or -
# Query chromosome
# Start on query
# End on query
# Bases in gap-free alignments
# Score - a number proportional to 100x matching bases

# The output is already ordered by tName, tStart

query="SELECT $columns FROM $table
        WHERE ( type != 'gap' AND tName = "'"'"$tname"'"'" AND
                (   $tstart <= tStart AND tStart < $tend
                 OR $tstart <= tEnd   AND tEnd   <= $tend))
        ORDER BY tStart;"

echo "$query" >&2

echo "$query" | mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A $db
# The -A flag is optional but is recommended for speed.
