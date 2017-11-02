#!/bin/bash

samples=$1
command=$2

known_commands=(
    view
    consensus
    chimeric
    mutation
    summary
)

isknown=0
for cmd in ${known_commands[@]}
do
    if [ "$command" == "$cmd" ] ; then
	isknown=1
	break
    fi
done

if [ "$isknown" == 0 ] ; then
    for cmd in ${known_commands[@]}
    do
	echo "usage: $0 samples.csv $cmd"
    done
    
    exit 1
fi

### define top project directory
PROJDIR=

cd $PROJDIR

while read line
do
    if ( [[ ! $line =~ "SampleID" ]] && [[ ! $line =~ "#" ]] )
    then
	### extract SampleID, Amplicon, and data path for each sample
	sampleId=`echo $line | cut -d, -f1`
	amplicon=`echo $line | cut -d, -f4`
	collectionPathUri=`echo $line | cut -d, -f5`
	
	if [ "$command" == "view" ]
	then
	    ### preview run info
	    echo ""
	    echo "sampleId=$sampleId"
	    echo "amplicon=$amplicon"
	    echo "collectionPathUri=$collectionPathUri"
	else
	    logdir=`printf "%s/samples/%05i" $PROJDIR $sampleId`
	    mkdir -p "$logdir"
    	    qsub -v root="$PROJDIR",amplicon="$amplicon",sampleId="$sampleId",collectionPathUri="$PROJDIR/$collectionPathUri" -N "ccs$sampleId" -o $logdir/$command.log -j yes $PROJDIR/scripts/ccs2-$command.sh
	fi
    fi
done < $samples
