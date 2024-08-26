#!/bin/bash
 
#SBATCH --job-name=combine-coal #Give your job a name. 
#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if comb to 2, could run two softwares in the script at the same time.
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

msmctools="/data/biosoftware/msmc-tools/msmc-tools"
data="./../data"
out="./../data/msmc_haploid2x_08hap/" 

## Compute with input -r 0.15  rho/theta
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do # 1 2 3 4 5 6 7 8 9 10
$msmctools/combineCrossCoal.py $out/aeg-wht.08hapl.allChr.across-pop12.comb$i.output.final.txt $out/aeg.08hapl.allChr.within-pop1.comb$i.output.final.txt $out/wht.08hapl.allChr.within-pop2.comb$i.output.final.txt > $out/combined_wht-aeg_08hapl_msmc.comb$i.r015.final.txt;
done
