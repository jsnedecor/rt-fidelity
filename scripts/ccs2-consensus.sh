#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 8

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load smrtlink-3.1.1
module load samtools-1.3.1

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

### create run directory
rundir=`printf "%s/samples/%05i" $root $sampleId`
echo $rundir
mkdir -p $rundir

### switch to working directory
cd $rundir

### look up sequencing data
echo ""
find "$collectionPathUri" -name "*.bax.h5" | sort -u > input.fofn
echo "Task 1 completed at $(date)"

### convert
echo ""

### check if FOFN contains filenames with spaces
num=`cat input.fofn | grep ' ' | wc -l`

if [ $num -gt 0 ]
then
    ### there are spaces, use workaround
    IFS=$'\r\n' GLOBIGNORE='*' command eval 'hdf=($(<input.fofn))'
    bax2bam "${hdf[@]}" -o movie
else
    ### no spaces, proceed as usual
    bax2bam --fofn=input.fofn -o movie
fi

echo "Task 2 completed at $(date)"

### map reads
echo ""
blasr movie.subreads.bam "$reference" -bam -out aligned_reads.bam -useQuality -nproc 8 -bestn 1
echo "Task 3 completed at $(date)"

### extract mapping direction
echo ""
$root/bin/ccs2-map.pl aligned_reads.bam | bzip2 - > clusters.csv.bz2
echo "Task 4 completed at $(date)"

### split forward and reverse reads
echo ""
$root/bin/ccs2-split.pl movie.subreads.bam clusters.csv.bz2
echo "Task 5 completed at $(date)"

### strand-specific CCS reads
for strand in fwd rev
do
    ### convert to bam
    echo ""
    samtools view -Sb subreads.${strand}.sam > subreads.${strand}.bam
    echo "Task 6 ($strand) completed at $(date)"
    
    ### build ccs
    echo ""
    ccs --reportFile=subreads_ccs.${strand}.csv --logFile=subreads_ccs.${strand}.log --numThreads=7 --minPasses=1 subreads.${strand}.bam subreads_ccs.${strand}.bam
    echo "Task 7 ($strand) completed at $(date)"
done

### cleanup
echo ""
rm -f movie.*
rm -f subreads.*
rm -f aligned_reads.*
rm -f clusters.csv.bz2
echo "Task 8 completed at $(date)"

echo ""
echo "Finished on $(date)"
