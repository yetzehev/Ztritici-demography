#!/bin/bash
#SBATCH --job-name=input-msmc #Give your job a name.
#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 #Multithreading.
#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=5G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.

out="./../data/pseudodiploid_vcf/"
msmcinput="./../data/msmc_16hap"

SAMPLE1="ZT442.ZT431" #aeg-pop1
SAMPLE2="ZT454.ZT429" #aeg-pop1
SAMPLE3="ZT460.ZT528" #aeg-pop1
SAMPLE4="ZT441.ZT498" #aeg-pop1
SAMPLE5="ZT617.ZT576" #wht-pop2
SAMPLE6="ZT567.ZT649" #wht-pop2
SAMPLE7="ZT611.ZT638" #wht-pop2
SAMPLE8="ZT676.ZT583" #wht-pop2

#This script generates the masking files for 8 pseudodiploid bamfiles created on step 1; four for Aegilops and 4 for Triticum
for CHR in NC_018218.1 NC_018217.1 NC_018216.1 NC_018215.1 NC_018214.1 NC_018213.1 NC_018212.1 NC_018211.1 NC_018210.1 NC_018209.1 NC_018208.1 NC_018207.1 NC_018206.1; do #NC_018205.1 NC_018204.1 NC_018203.1 NC_018202.1 NC_018201.1 NC_018200.1 NC_018199.1 NC_018198.1; do

../../packages/msmc-tools/generate_multihetsep.py \
--mask=$out/$SAMPLE1.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE2.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE3.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE4.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE5.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE6.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE7.mask.chr$CHR.bed.gz \
--mask=$out/$SAMPLE8.mask.chr$CHR.bed.gz \
$out/$SAMPLE1.chr$CHR.vcf.gz $out/$SAMPLE2.chr$CHR.vcf.gz $out/$SAMPLE3.chr$CHR.vcf.gz $out/$SAMPLE4.chr$CHR.vcf.gz $out/$SAMPLE5.chr$CHR.vcf.gz $out/$SAMPLE6.chr$CHR.vcf.gz $out/$SAMPLE7.chr$CHR.vcf.gz $out/$SAMPLE8.chr$CHR.vcf.gz> $msmcinput/aeg-wht.16hapl.comb1.chr.$CHR.msmc.input;done
