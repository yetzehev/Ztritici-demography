#!/bin/bash
 
#SBATCH --job-name=cross-pop2 #Give your job a name.
#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 #Multithreading.
#SBATCH --time=48:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=24G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.

#This script uses as input two pseudo-diploid files to runs msmc for the 13 chromosomes of Z. tritici

for i in 1 2 3 4 5 6 7 8 9 10; do

out="./../data/msmc_16hap"

msmc2 -I 0,1,2,3,4,5,6,7 -o $out/aeg.16hapl.allChr.within-pop1.comb$i.output $out/aeg-wht.16hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018218.1.msmc.input

msmc2 -I 8,9,10,11,12,13,14,15 -o $out/wht.16hapl.allChr.within-pop2.comb$i.output $out/aeg-wht.16hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018218.1.msmc.input

msmc2 -I 0-8,0-9,0-10,0-11,0-12,0-13,0-14,0-15,1-8,1-9,1-10,1-11,1-12,1-13,1-14,1-15,2-8,2-9,2-10,2-11,2-12,2-13,2-14,2-15,3-8,3-9,3-10,3-11,3-12,3-13,3-14,3-15,4-8,4-9,4-10,4-11,4-12,4-13,4-14,4-15,5-8,5-9,5-10,5-11,5-12,5-13,5-14,5-15,6-8,6-9,6-10,6-11,6-12,6-13,6-14,6-15,7-8,7-9,7-10,7-11,7-12,7-13,7-14,7-15 -o $out/aeg-wht.16hapl.allChr.across-pop12.comb$i.output $out/aeg-wht.16hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018218.1.msmc.input ;
rm $out/aeg-wht.16hapl.comb$i.chr.NC_*.1.msmc.input;
done




