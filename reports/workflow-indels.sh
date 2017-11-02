#!/bin/bash

PROJDIR=

samples="$PROJDIR/samples.csv"
sampledir="$PROJDIR/samples"
scriptdir="$PROJDIR/reports/bin"
resultdir="$PROJDIR/reports/indels"

enzymes=("ProtoScript II Reverse Transcriptase" "M-MuLV Reverse Transcriptase" "AMV Reverse Transcriptase" "Bst 2.0 DNA Polymerase" "Bst 3.0 DNA Polymerase")

mkdir -p "$resultdir"

cd "$resultdir"

(echo "Enzyme,Template,Amplicon,ST,S1,S2,S3,S4,SS,DT,D1,D2,D3,D4,DS,IT,I1,I2,I3,I4,IS"
    for enzyme in "${enzymes[@]}"
    do
	for template in "Regular" "N6mA" "P-U" "5mC" "5mU" "5hmU"
	do
    	    head --lines 1 "$samples" > samples.tmp
    	    cat "$samples" | grep "$enzyme" | grep "$template" >> samples.tmp
    	    echo -n "$enzyme,$template,NA,"
    	    $scriptdir/mutation-map.pl --mode FSE samples.tmp "$sampledir"
	done
    done
) > table-indels-fse.csv

(echo "Enzyme,Template,Amplicon,ST,S1,S2,S3,S4,SS,DT,D1,D2,D3,D4,DS,IT,I1,I2,I3,I4,IS"
    for enzyme in "${enzymes[@]}"
    do
	for template in "Regular"
	do
    	    head --lines 1 "$samples" > samples.tmp
    	    cat "$samples" | grep "$enzyme" | grep "$template" >> samples.tmp
    	    echo -n "$enzyme,$template,NA,"
    	    $scriptdir/mutation-map.pl --mode SSE samples.tmp "$sampledir"
	done
    done
) > table-indels-sse.csv

rm -f "$resultdir"/*.tmp
