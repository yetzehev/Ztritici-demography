#!/bin/bash
#SBATCH --job-name=call-mask #Give your job a name.
#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=12G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.
#This script generates VCF and mask files from individual diploid BAM files.

data="./../data/haploid"
out="./../data/haploid2x/"
sample="ZT549"
msmctools="../../packages/msmc-tools/"
refGenome="../data/ncbi-genomes-2022-09-28/GCF_000219625.1_MYCGR_v2.0_genomic.fna"
#estimating the average sequencing depth using all site on chromosome 2 = NC_018217.1

DEPTH=$(samtools depth -r NC_018218.1 $data/$sample.ipo323.rm.dp.rg.bam | awk '{sum += $3} END {print sum / NR}')
echo $DEPTH > $out/tmp
for CHR in NC_018218.1 NC_018217.1 NC_018216.1 NC_018215.1 NC_018214.1 NC_018213.1 NC_018212.1 NC_018211.1 NC_018210.1 NC_018209.1 NC_018208.1 NC_018207.1 NC_018206.1 ; do
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $CHR -f $refGenome $data/$sample.ipo323.rm.dp.rg.bam | bcftools call -c -V indels | $msmctools/bamCaller.py $DEPTH $out/$sample.mask.chr$CHR.bed.gz | gzip -c > $out/$sample.chr$CHR.vcf.gz;
done

