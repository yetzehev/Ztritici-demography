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

#This script generates a executing file to create pseudodiploid genotypes. It is based on the post "How to merge two haploid samples into a pseudo-diploid" available at "https://www.biostars.org/p/341696/"

#1) Combine the 45 genotypes in pairs separated by a comma
#cat Ztritici.aeg.txt | while read s; do for d in `cat Ztritici.aeg.txt`; do if [ $s != $d ]; then echo $s,$d; fi; done; done > Ztritici.aeg.diploid.comb

#Then make a list with 100 random combinations 
#shuf -n 100 Ztritici.aeg.diploid.comb > Ztritici.aeg.diploid.100comb; tr "\n" " " < Ztritici.aeg.diploid.100comb


#Instructions
#1) Execute this file with ./make_stp_01_pseudoDipl_aegCyl.sh to generate the bash files to submit 
#2) Submit each task on wallace as an individual task

#=================================
# Aegilops
#=================================

for i in ZT442,ZT431 ZT500,ZT427 ZT500,ZT440 ZT431,ZT460 ZT460,ZT508 ZT442,ZT476 ZT453,ZT442 ZT536,ZT475 ZT441,ZT498 ZT431,ZT536 ZT442,ZT429 ZT495,ZT517 ZT485,ZT504 ZT476,ZT479 ZT481,ZT508 ZT460,ZT528 ZT517,ZT537 ZT481,ZT534 ZT454,ZT429 ZT536,ZT485 ZT534,ZT485 ZT479,ZT440 ZT460,ZT537 ZT448,ZT427 ZT536,ZT440 ZT485,ZT537 ZT431,ZT485 ZT508,ZT441 ZT481,ZT429 ZT508,ZT498 ZT485,ZT469 ZT476,ZT429 ZT508,ZT442 ZT536,ZT453 ZT454,ZT485 ZT469,ZT481 ZT508,ZT460 ZT534,ZT475 ZT469,ZT427 ZT436,ZT454 ZT485,ZT427 ZT534,ZT442 ZT537,ZT534 ZT453,ZT508 ZT536,ZT454 ZT427,ZT515 ZT508,ZT489 ZT440,ZT448 ZT495,ZT476 ZT517,ZT489 ZT454,ZT453 ZT440,ZT500 ZT495,ZT498 ZT440,ZT534 ZT448,ZT431 ZT517,ZT442 ZT498,ZT453 ZT528,ZT500 ZT427,ZT517 ZT537,ZT479 ZT534,ZT500 ZT508,ZT528 ZT508,ZT504 ZT460,ZT441 ZT427,ZT475 ZT427,ZT489 ZT431,ZT454 ZT515,ZT517 ZT537,ZT448 ZT436,ZT460 ZT427,ZT500 ZT481,ZT537 ZT440,ZT495 ZT515,ZT453 ZT504,ZT460 ZT528,ZT481 ZT489,ZT481 ZT537,ZT476 ZT508,ZT453 ZT517,ZT528 ZT489,ZT436 ZT429,ZT536 ZT489,ZT460 ZT431,ZT436 ZT498,ZT515 ZT528,ZT475 ZT454,ZT528 ZT476,ZT442 ZT429,ZT448 ZT481,ZT515 ZT537,ZT515 ZT427,ZT537 ZT517,ZT427 ZT534,ZT489 ZT479,ZT442 ZT498,ZT469 ZT495,ZT429 ZT479,ZT453 ZT469,ZT440 ZT440,ZT453 ; do    KEY=${i%,*};   VAL=${i#*,};  
echo -e '#!/bin/bash'"\n#SBATCH --job-name=psuedoDipl #Give your job a name.\n#SBATCH --nodes=1 #Only increase for openmpi jobs.\n#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1 #Multithreading.\n#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.\n#SBATCH --mem=8G #Total Memory per node to use for the job\n#SBATCH --error=job.%J.err #Std Error write standard error to this file\n#SBATCH --output=job.%J.out #Std Out write standard output to this file\n#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)\n#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line\n#SBATCH --partition=standard #Request a specific partition for the resource allocation.\n#This script takes merge two haploid files and update the index header to createa new pseudosiploid file\n\ndata=\"./../data/haploid\"\nout=\"./../data/pseudodiploid\"""\nSAMPLE1=\""$KEY"\"\nSAMPLE2=\""$VAL"\"\nsamtools merge \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam \$data/\$SAMPLE1.ipo323.rm.dp.rg.bam \$data/\$SAMPLE2.ipo323.rm.dp.rg.bam\njava -jar picard.jar AddOrReplaceReadGroups \\\\\nI=\$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam \\\\\nO=\$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.readGroup.bam \\\\\nRGID=\$SAMPLE1.\$SAMPLE2 \\\\\nRGLB=lib1 \\\\\nRGPL=illumina \\\\\nRGPU=unit1 \\\\\nRGSM=\$SAMPLE1.\$SAMPLE2 \\\\\n\nsamtools index -b \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.readGroup.bam \\\\\n\nrm \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam" > stp_01_pseudoDipl_$KEY-$VAL-aeg.sh;done


