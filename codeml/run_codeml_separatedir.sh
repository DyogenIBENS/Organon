#!/bin/bash

# Tell the script to quit at the first command that fails
set -e
# For the same behavior with every command in a pipe (not wanted here):
#set -o pipefail
# Error when trying to use unset variable
#set -u

help='
DESCRIPTION

Run codeml in a separate directory, then move files in the current dir
(output files with standard names (rub, rst, 2NG.dS, etc) will be prefixed
by the basename of the control file)

USAGE

To run the control file `ENSGT00790000122969.ctl`, do:

    ./run_codeml_separatedir.sh ENSGT00790000122969.ctl

If the control file is in a directory (`subdir/ENSGT00790000122969.ctl`), do:

    ./run_codeml_separatedir.sh subdir/ENSGT00790000122969.ctl

results will then appear in `subdir`.

OPTIONS
    -h  Print help and exit.
    -v  Use valgrind --tool=massif to measure memory consumption
        (slows down the program). Output to `$genetree.massif.out`
    -b  Disable the use of `unbuffer` in front of the command line
        (will delay stdout messages).
    -q  Quiet mode.
    -s  Keep the output in its separate subdirectory.
'
### PARSE COMMAND LINE ###

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
valgrind=0
buffer=0
quiet=0
subdir=0

while getopts "h?vbqs" opt; do
	case "$opt" in
	h|\?)
		echo "$help"
		exit 0
		;;
	v)  valgrind=1
		;;
	b)  buffer=1
		;;
	q)  quiet=1
		;;
	s)  subdir=1
		;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [ -z "$1" ]; then
	echo "$help"
	exit 1
fi

### Set up command executed in case of error
function verbose_exit {
	msg="Exit information:
    current directory: $(pwd -P)
    current ctl file: ${genetree:-}.${ctlext:-}
    seqfile: ${seqfile:-}
    treefile: ${treefile:-}
    outfile: ${outfile:-}"
	echo "$msg" >&2
}

trap verbose_exit EXIT

### SET UP FILENAMES AND TEMPORARY DIRECTORY ###

# remove extension
ctlfile=$1
genetree=${1%.*}
ctlext=${1##*.}

# parse seqfile, treefile and outfile from the control file:
seqfile=`sed -rn 's/^\s*seqfile\s*=\s*(\S+)[\s*]?.*$/\1/p'   $ctlfile`
treefile=`sed -rn 's/^\s*treefile\s*=\s*(\S+)[\s*]?.*$/\1/p' $ctlfile`
outfile=`sed -rn 's/^\s*outfile\s*=\s*(\S+)[\s*]?.*$/\1/p'   $ctlfile`

genetreedir=$(dirname $genetree)
genetree=$(basename $genetree)

# Go to the directory where the ctl file is
cd $genetreedir

# Create a temporary directory for output files.
if [ ! -d "$genetree" ]; then
	mkdir "$genetree"
fi

# convert to absolute path. Needed in the temporary ctl file.
# Function from: http://stackoverflow.com/a/13087801/4614641
# Works only for existing files.
function abspath {
    if [[ -d "$1" ]]
    then
        pushd "$1" >/dev/null
        pwd
        popd >/dev/null
    elif [[ -e $1 ]]
    then
        pushd "$(dirname "$1")" >/dev/null
        echo "$(pwd)/$(basename "$1")"
        popd >/dev/null
    else
        echo "$1" does not exist! >&2
        return 127
    fi
}

seqfile=$(abspath $seqfile)
treefile=$(abspath $treefile)
outbase=$(basename $outfile)
outdir=$(abspath $(dirname $outfile))
outfile="$outdir/$outbase"

# To avoid bad surprises, link the input files in the current dir. That way, 
# the filename length cannot exceed 128 characters (max size allowed by codeml)
#echo "Linking seqfile and treefile in temporary directory." >&2
ln -fs -T "$seqfile" "$genetree/seqfile.phy"
ln -fs -T "$treefile" "$genetree/treefile.nwk"

files="seqfile  = seqfile.phy
treefile = treefile.nwk
outfile  = outfile.mlc"
#echo "$files" >&2

# Create a temporary control file in the temporary directory, containing the 
# absolute paths to the input and output files
tmpctl="$genetree/$genetree.$ctlext"
echo "$files" > "$tmpctl"
sed -rn '/^\s*(seqfile|treefile|outfile)\s*=/!p' "$genetree.$ctlext" >> "$tmpctl"

cd "$genetree"

### Execute CODEML ###
if [ "$valgrind" -eq 1 ]; then
	valgrind --tool=massif --massif-out-file=../"$genetree.massif.out" \
		codeml "$genetree.$ctlext"
elif [ "$buffer" -eq 1 ]; then
	codeml "$genetree.$ctlext" | tee ../"$genetree.log"
elif [ "$quiet" -eq 1 ]; then
	codeml "$genetree.$ctlext" > ../"$genetree.log"
else
	unbuffer codeml "$genetree.$ctlext" | tee ../"$genetree.log"
fi
# TODO: build command as string, then eval.


### Rename, move and clean up ###
if [ ${PIPESTATUS[0]} -eq 0 ]; then
	if [ $subdir -eq 0 ]; then
		for f in "2NG.dS" "2NG.dN" "2NG.t" "rub" "rst" "rst1" "lnf" "4fold.nuc"; do
			mv $f ../${genetree}_$f
		done
	
		# remove the temporary control file
		rm "$genetree.$ctlext"
		# remove symbolic links (to seqfile and treefile)
		#find -type l -delete
		rm seqfile.phy treefile.nwk
		
		# move the result file if present in the output directory
		#maybeoutfile=$(basename $outfile)
		#if [ -f $maybeoutfile ]; then
		#	# check if need overwriting an existing file
		#	#if [ -f ../$maybeoutfile ]; then ;fi
		#	mv $maybeoutfile ../$maybeoutfile
		#fi
		mv outfile.mlc "$outfile"

		cd ..
		rmdir $genetree/ #|| echo >&2 "$!"
	fi
else
	echo "codeml failed!" >&2
	exit 1
fi

# Unset the trap
trap - EXIT
