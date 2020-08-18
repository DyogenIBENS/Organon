#!/bin/bash

set -e
set -u
IFS=$'\n\t'
trap "exit 1" INT

echo -e 'name\tseq_bytes\tseq_nb\tseq_len\ttime\tmemory'

# Number of leaves
for ctlfile in *.ctl; do
	base=`basename "$ctlfile" .ctl`
	outfile="$base.condor.stdout"
	logfile="$base.condor.log"
	if [ -f "$outfile" -a -f "$logfile" ]; then
		#treefile=`sed -rn 's/^\s*treefile\s*=\s*(\S+)[\s*]?.*$/\1/p' "$ctlfile"`
		seqfile=`sed -rn 's/^\s*seqfile\s*=\s*(\S+)[\s*]?.*$/\1/p' "$ctlfile"`
		if [ -z "$seqfile" ]; then
			echo "Error: 'seqfile' not found in $ctlfile" >&2
			exit 1
		fi

		# Get: bytesize of alignment, nb of sequences, length of sequences
		echo -ne "$base\t"
		echo -ne "`wc -c $seqfile | cut -d' ' -f1`\t"
		echo -ne "`head -1 $seqfile | cut -f-2`\t"
		
		if grep -q 'Job terminated' "$logfile"; then
			
			# Time used:
			grep '^Time used:' "$outfile" | \
				sed 's/^Time used:\s\+//' | \
				tr '\n' '\t'

			# Memory used (MB)
			# FIXME: the log file (from condor) is appended at every run.
			#        therefore you need to grep the last occurrence of 'Memory'
			grep '^\s\+Memory' "$logfile" | \
				sed -r 's/^\s+Memory \(MB\)\s+:\s+([0-9]+)\s+[0-9]+\s+[0-9]+$/\1/'
		else
			echo -e '\t'
		fi
	else
		echo "Codeml Job not started for $ctlfile" >&2
	fi
done

