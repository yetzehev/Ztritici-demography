#!/bin/bash
 
#SBATCH --job-name=make-stp_01 #Give your job a name.
#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 #Multithreading.
#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=8G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.

#Based on the post "How to merge two haploid samples into a pseudo-diploid" available at "https://www.biostars.org/p/341696/"

#Note: This script generates the bash files to make the pseudodiploid bam file for Triticum aestivum isolates
#Ztritici.wht.txt | while read s; do for d in `cat Ztritici.wht.txt`; do if [ $s != $d ]; then echo $s,$d; fi; done; done > Ztritici.wht.diploid.comb
#Generate a list of 100 random combinations with the comand below, those combination will be used for the SNP calling
#shuf -n 100 Ztritici.wht.diploid.comb > Ztritici.wht.diploid.100comb; tr "\n" " " < Ztritici.wht.diploid.100comb


#Instructions
#1) Execute this file with ./make_stp_01_pseudoDipl_TrtAes.sh to generate the bash files to submit 
#2) Submit each task on wallace as an individual task

#=================================
# Triticum
#=================================

for i in ZT617,ZT576 ZT567,ZT649 ZT611,ZT638 ZT676,ZT583 ZT705,ZT576 ZT640,ZT671 ZT661,ZT635 ZT611,ZT650 ZT643,ZT662 ZT640,ZT559 ZT663,ZT611 ZT704,ZT611 ZT662,ZT642 ZT565,ZT655 ZT573,ZT705 ZT710,ZT709 ZT704,ZT643 ZT671,ZT661 ZT565,ZT709 ZT611,ZT710 ZT709,ZT644 ZT721,ZT681 ZT640,ZT676 ZT587,ZT710 ZT555,ZT668 ZT678,ZT643 ZT567,ZT555 ZT678,ZT699 ZT559,ZT570 ZT642,ZT699 ZT705,ZT699 ZT655,ZT635 ZT661,ZT649 ZT549,ZT651 ZT549,ZT559 ZT555,ZT570 ZT559,ZT616 ZT583,ZT573 ZT555,ZT638 ZT681,ZT634 ZT567,ZT640 ZT682,ZT655 ZT620,ZT721 ZT716,ZT709 ZT704,ZT682 ZT634,ZT721 ZT676,ZT638 ZT635,ZT642 ZT611,ZT549 ZT583,ZT651 ZT709,ZT721 ZT555,ZT682 ZT721,ZT555 ZT649,ZT635 ZT559,ZT681 ZT627,ZT635 ZT650,ZT627 ZT620,ZT649 ZT617,ZT559 ZT640,ZT649 ZT583,ZT644 ZT663,ZT567 ZT567,ZT663 ZT640,ZT587 ZT671,ZT682 ZT643,ZT611 ZT644,ZT642 ZT634,ZT671 ZT655,ZT612 ZT650,ZT655 ZT655,ZT616 ZT634,ZT576 ZT576,ZT672 ZT662,ZT611 ZT650,ZT634 ZT671,ZT676 ZT576,ZT634 ZT649,ZT620 ZT555,ZT661 ZT555,ZT705 ZT705,ZT612 ZT627,ZT573 ZT709,ZT620 ZT705,ZT642 ZT617,ZT573 ZT644,ZT671 ZT651,ZT721 ZT668,ZT573 ZT616,ZT668 ZT638,ZT611 ZT650,ZT616 ZT672,ZT617 ZT649,ZT678 ZT668,ZT662 ZT650,ZT555 ZT617,ZT627 ZT643,ZT661 ZT655,ZT643 ZT663,ZT676 ZT612,ZT640 ; do    KEY=${i%,*};   VAL=${i#*,};  
echo -e '#!/bin/bash'"\n#SBATCH --job-name=psuedoDipl #Give your job a name.\n#SBATCH --nodes=1 #Only increase for openmpi jobs.\n#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1 #Multithreading.\n#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.\n#SBATCH --mem=8G #Total Memory per node to use for the job\n#SBATCH --error=job.%J.err #Std Error write standard error to this file\n#SBATCH --output=job.%J.out #Std Out write standard output to this file\n#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)\n#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line\n#SBATCH --partition=standard #Request a specific partition for the resource allocation.\n#This script takes merge two haploid files and update the index header to createa new pseudosiploid file\n\ndata=\"./../data/haploid\"\nout=\"./../data/pseudodiploid\"""\nSAMPLE1=\""$KEY"\"\nSAMPLE2=\""$VAL"\"\nsamtools merge \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam \$data/\$SAMPLE1.ipo323.rm.dp.rg.bam \$data/\$SAMPLE2.ipo323.rm.dp.rg.bam\njava -jar picard.jar AddOrReplaceReadGroups \\\\\nI=\$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam \\\\\nO=\$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.readGroup.bam \\\\\nRGID=\$SAMPLE1.\$SAMPLE2 \\\\\nRGLB=lib1 \\\\\nRGPL=illumina \\\\\nRGPU=unit1 \\\\\nRGSM=\$SAMPLE1.\$SAMPLE2 \\\\\n\nsamtools index -b \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.readGroup.bam \\\\\n\nrm \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam" > stp_01_pseudoDipl_$KEY-$VAL- wht.sh;done


