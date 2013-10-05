#!/bin/sh
#PBS -l nodes=1:ppn=4
#PBS -l mem=4gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

perl SmallRNA.pl --read SRR771499_1.fastq.gz --step 1234 --cpu 4

echo "qsub Done"

