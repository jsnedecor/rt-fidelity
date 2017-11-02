#!/bin/bash

#$ -S /bin/bash

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load p7zip

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

### run directory
rundir=`printf "%s/samples/%05i/summary" $root $sampleId`
mkdir -p $rundir
cd $rundir

### tally up mutations
echo ""
sampledir=`printf "%s/samples/%05i" $root $sampleId`

n=9     # number of chunks
c=20000 # number of Zmws per chunk

for x in `seq 1 $n`
do
    ### define range
    b=`echo "scale=0; ($x-1)*$c+1" | bc`
    e=`echo "scale=0; $x*$c"       | bc`
    range="$b-$e"
    
    pbsfile="$rundir/summary.pbs.$x"
    logfile="$rundir/summary.log.$x"

    rm -f $logfile

    cat <<EOF >$pbsfile
#!/bin/bash

#\$ -S /bin/bash
#\$ -P longrun
#\$ -o $logfile
#\$ -j yes
#\$ -pe smp 2
#\$ -cwd

echo ""
echo "Started on \$(date) on \$(uname -n)"

module load p7zip

$root/bin/ccs2-summary.pl \\
    --range $range \\
    --strand fwd,rev \\
    --mutation-dir "mutation" \\
    --np 15 \\
    --qv 93 \\
    --mapq 60 \\
    --lb 40 \\
    --ub 40 \\
    --rlen-min-a 50 \\
    --rlen-max-a 50 \\
    --ctxfile context.csv.$x \\
    $reference $sampledir >summary.csv.$x

echo ""
echo "Finished on \$(date)"

EOF

    echo $pbsfile
    qsub $pbsfile
done

echo ""
echo "Finished on $(date)"
