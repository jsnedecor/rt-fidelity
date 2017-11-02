#!/bin/bash

PROJDIR=

samples="$PROJDIR/samples.csv"
sampledir="$PROJDIR/samples"
scriptdir="$PROJDIR/reports/bin"
resultdir="$PROJDIR/reports/summary"

mkdir -p "$resultdir"

cd "$resultdir"

for mode in fse sse
do
    for x in $(ls "$sampledir")
    do
	### combine summary.csv.*
	"$scriptdir"/merge.pl --mode $mode "$sampledir"/$x/summary/summary.csv.* | cut -d, -f1-26 >> "$resultdir"/summary-norm.tmp
    done

    ### reformat summary
    ( head --lines 1 "$resultdir"/summary-norm.tmp && egrep -v "^NP" "$resultdir"/summary-norm.tmp ) > "$resultdir"/summary-norm.2.tmp

    ### reformat samples
    cut -d, -f1,2,3,4 "$samples" > "$resultdir"/samples.tmp

    ### combine pieces
    paste -d, "$resultdir"/{samples.tmp,summary-norm.2.tmp} > "$resultdir"/summary.$mode.csv

    ### cleanup
    rm -f "$resultdir"/*.tmp

    ### generate tables
    Rscript "$scriptdir"/table.R "$resultdir"/summary.$mode.csv "$resultdir"/table $mode
done
