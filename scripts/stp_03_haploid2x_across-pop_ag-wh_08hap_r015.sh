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
#SBATCH --partition=global #Request a specific partition for the resource allocation.

## This script uses as input two homozygous-diploid files to runs msmc for the 13 chromosomes of Z. tritici
## since files are homodiploid I'm only comparing only one genotype per file

for i in 1 2 3 4 5 6 7 8 9 10  11 12 13 14 15 16 17 18 19 20; do
out="./../data/msmc_haploid2x_08hap/"

## Note: Only one genome index was used per input file. For example, in the first file, genome-index 0 was included while genome-index 1 was not. In the second file, genome-index 2 was included, and the following index was excluded, and so on.

msmc2 -I 0,2,4,6 -r 0.15 -o $out/aeg.08hapl.allChr.within-pop1.comb$i.output $out/aeg-wht.08hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018218.1.msmc.input

msmc2 -I 8,10,12,14 -r 0.15 -o $out/wht.08hapl.allChr.within-pop2.comb$i.output $out/aeg-wht.08hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018218.1.msmc.input

msmc2 -I 0-8,0-10,0-12,0-14,2-8,2-10,2-12,2-14,4-8,4-10,4-12,4-14,6-8,6-10,6-12,6-14 -r 0.15 -o $out/aeg-wht.08hapl.allChr.across-pop12.comb$i.output $out/aeg-wht.08hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.08hapl.comb$i.chr.NC_018218.1.msmc.input ;

rm $out/aeg-wht.08hapl.comb$i.chr.NC_*.1.msmc.input;

done




