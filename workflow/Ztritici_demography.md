# Demography of two divergent populations of *Triticum*- and *Aegilops*-infecting  *Z. tritici* 

# Pipeline description

###  Aim
To assess the  demographic history of *Aegilops*- and *Triticum*-infecting *Z. tritici* populations

### Population divergence scenarios of *Aegilops*- and *Triticum*-infecting Z. tritici

A. Compute population divergence with different number of time segments: 32, 16 and 12

B. Compute divergence with diferent genomes subsets

## Data sets

A total number of 74 genotypes:
- *Triticum*-infecting  infecting Z.tritici genotypes (n = 45)
- *Aegilops cylindrica*-infecting Z. tritici genotypes (n= 29)

Note: Outiers isolates that infect *A. tauschii* and high IBS *Triticum*-infecting isolates were not included

### Step 1: Pseudidiploid bam files

B.Define the 100 pseudodiploid genotypes to perform the analysis. I randomly combined the IDs of 29 *Aegilops*, or the IDs of the 45 *Triticum*-infecting isolates.
	
    #Combine the IDs in pairs
    cat Ztritici.aeg.txt | while read s; do for d in `cat Ztritici.aeg.txt`; do if [ $s != $d ]; then echo $s,$d; fi; done; done > Ztritici.aeg.diploid.comb
	#Generate a list with 100 random combinations the comand 
	shuf -n 100 Ztritici.aeg.diploid.comb > Ztritici.aeg.diploid.100comb; tr "\n" " " < Ztritici.aeg.diploid.100comb
    
C. The list of pseudidiploid combination was used to generate pseudodiploid bam files. Execute the script below with the command 

	./make_stp_01_pseudoDipl_aegCyl.sh 


To display the content of the script with the command
	
    cat ./make_stp_01_pseudoDipl_aegCyl.sh
    

This command will show:
```
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


    for i in ZT442,ZT431 ZT500,ZT427 ZT500,ZT440 ZT431,ZT460 ZT460,ZT508 ZT442,ZT476 ZT453,ZT442 ZT536,ZT475 ZT441,ZT498 ZT431,ZT536 ZT442,ZT429 ZT495,ZT517 ZT485,ZT504 ZT476,ZT479 ZT481,ZT508 ZT460,ZT528 ZT517,ZT537 ZT481,ZT534 ZT454,ZT429 ZT536,ZT485 ZT534,ZT485 ZT479,ZT440 ZT460,ZT537 ZT448,ZT427 ZT536,ZT440 ZT485,ZT537 ZT431,ZT485 ZT508,ZT441 ZT481,ZT429 ZT508,ZT498 ZT485,ZT469 ZT476,ZT429 ZT508,ZT442 ZT536,ZT453 ZT454,ZT485 ZT469,ZT481 ZT508,ZT460 ZT534,ZT475 ZT469,ZT427 ZT436,ZT454 ZT485,ZT427 ZT534,ZT442 ZT537,ZT534 ZT453,ZT508 ZT536,ZT454 ZT427,ZT515 ZT508,ZT489 ZT440,ZT448 ZT495,ZT476 ZT517,ZT489 ZT454,ZT453 ZT440,ZT500 ZT495,ZT498 ZT440,ZT534 ZT448,ZT431 ZT517,ZT442 ZT498,ZT453 ZT528,ZT500 ZT427,ZT517 ZT537,ZT479 ZT534,ZT500 ZT508,ZT528 ZT508,ZT504 ZT460,ZT441 ZT427,ZT475 ZT427,ZT489 ZT431,ZT454 ZT515,ZT517 ZT537,ZT448 ZT436,ZT460 ZT427,ZT500 ZT481,ZT537 ZT440,ZT495 ZT515,ZT453 ZT504,ZT460 ZT528,ZT481 ZT489,ZT481 ZT537,ZT476 ZT508,ZT453 ZT517,ZT528 ZT489,ZT436 ZT429,ZT536 ZT489,ZT460 ZT431,ZT436 ZT498,ZT515 ZT528,ZT475 ZT454,ZT528 ZT476,ZT442 ZT429,ZT448 ZT481,ZT515 ZT537,ZT515 ZT427,ZT537 ZT517,ZT427 ZT534,ZT489 ZT479,ZT442 ZT498,ZT469 ZT495,ZT429 ZT479,ZT453 ZT469,ZT440 ZT440,ZT453 ; do    KEY=${i%,*};   VAL=${i#*,};
    echo -e '#!/bin/bash'"\n#SBATCH --job-name=psuedoDipl #Give your job a name.\n#SBATCH --nodes=1 #Only increase for openmpi jobs.\n#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1 #Multithreading.\n#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.\n#SBATCH --mem=8G #Total Memory per node to use for the job\n#SBATCH --error=job.%J.err #Std Error write standard error to this file\n#SBATCH --output=job.%J.out #Std Out write standard output to this file\n#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)\n#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line\n#SBATCH --partition=standard #Request a specific partition for the resource allocation.\n#This script takes merge two haploid files and update the index header to createa new pseudosiploid file\n\ndata=\"./../data/haploid\"\nout=\"./../data/pseudodiploid\"""\nSAMPLE1=\""$KEY"\"\nSAMPLE2=\""$VAL"\"\nsamtools merge \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam \$data/\$SAMPLE1.ipo323.rm.dp.rg.bam \$data/\$SAMPLE2.ipo323.rm.dp.rg.bam\njava -jar picard.jar AddOrReplaceReadGroups \\\\\nI=\$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam \\\\\nO=\$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.readGroup.bam \\\\\nRGID=\$SAMPLE1.\$SAMPLE2 \\\\\nRGLB=lib1 \\\\\nRGPL=illumina \\\\\nRGPU=unit1 \\\\\nRGSM=\$SAMPLE1.\$SAMPLE2 \\\\\n\nsamtools index -b \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.readGroup.bam \\\\\n\nrm \$out/\$SAMPLE1.\$SAMPLE2.pseudoDiploid.bam" > stp_01_pseudoDipl_$KEY-$VAL-aeg.sh;done
    
```

D. Submit the individual files  generated as independent jobs.

	for i in *aeg.sh; do sbatch $i;done
    

### Step2: Perform the SNP calling with bcftools

A. After generating the pseudodiploid bamfiles, use the as input files for the SNPcallin with bcftools. This scrpit will generate an individual file for each pseudiploid genotype.

B. Execute with the command

	02_make_stp_02_call_mask_wht.sh
    
C. To visualize the content of the script execute:
	
    cat 02_make_stp_02_call_mask_wht.sh

This command will show:



	 
    #!/bin/bash
    #SBATCH --job-name=make-stp_02 #Give your job a name.
    #SBATCH --nodes=1 #Only increase for openmpi jobs.
    #SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1 #Multithreading.
    #SBATCH --time=48:00:00 #Time for a job to run given as hh:mm:ss.
    #SBATCH --mem=8G #Total Memory per node to use for the job
    #SBATCH --error=job.%J.err #Std Error write standard error to this file
    #SBATCH --output=job.%J.out #Std Out write standard output to this file
    #SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
    #SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
    #SBATCH --partition=standard #Request a specific partition for the resource allocation.

    #This script generates a bash script to perform the SNPcalling and masking for the pseudodiploid files.

    for i in ZT442,ZT431 ZT500,ZT427 ZT500,ZT440 ZT431,ZT460 ZT460,ZT508 ZT442,ZT476 ZT453,ZT442 ZT536,ZT475 ZT441,ZT498 ZT431,ZT536 ZT442,ZT429 ZT495,ZT517 ZT485,ZT504 ZT476,ZT479 ZT481,ZT508 ZT460,ZT528 ZT517,ZT537 ZT481,ZT534 ZT454,ZT429 ZT536,ZT485 ZT534,ZT485 ZT479,ZT440 ZT460,ZT537 ZT448,ZT427 ZT536,ZT440 ZT485,ZT537 ZT431,ZT485 ZT508,ZT441 ZT481,ZT429 ZT508,ZT498 ZT485,ZT469 ZT476,ZT429 ZT508,ZT442 ZT536,ZT453 ZT454,ZT485 ZT469,ZT481 ZT508,ZT460 ZT534,ZT475 ZT469,ZT427 ZT436,ZT454 ZT485,ZT427 ZT534,ZT442 ZT537,ZT534 ZT453,ZT508 ZT536,ZT454 ZT427,ZT515 ZT508,ZT489 ZT440,ZT448 ZT495,ZT476 ZT517,ZT489 ZT454,ZT453 ZT440,ZT500 ZT495,ZT498 ZT440,ZT534 ZT448,ZT431 ZT517,ZT442 ZT498,ZT453 ZT528,ZT500 ZT427,ZT517 ZT537,ZT479 ZT534,ZT500 ZT508,ZT528 ZT508,ZT504 ZT460,ZT441 ZT427,ZT475 ZT427,ZT489 ZT431,ZT454 ZT515,ZT517 ZT537,ZT448 ZT436,ZT460 ZT427,ZT500 ZT481,ZT537 ZT440,ZT495 ZT515,ZT453 ZT504,ZT460 ZT528,ZT481 ZT489,ZT481 ZT537,ZT476 ZT508,ZT453 ZT517,ZT528 ZT489,ZT436 ZT429,ZT536 ZT489,ZT460 ZT431,ZT436 ZT498,ZT515 ZT528,ZT475 ZT454,ZT528 ZT476,ZT442 ZT429,ZT448 ZT481,ZT515 ZT537,ZT515 ZT427,ZT537 ZT517,ZT427 ZT534,ZT489 ZT479,ZT442 ZT498,ZT469 ZT495,ZT429 ZT479,ZT453 ZT469,ZT440 ZT440,ZT453 ; do echo -e '#!/bin/bash'"\n#SBATCH --job-name=call-mask #Give your job a name.\n#SBATCH --nodes=1 #Only increase for openmpi jobs.\n#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1 #Multithreading.\n#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.\n#SBATCH --mem=12G #Total Memory per node to use for the job\n#SBATCH --error=job.%J.err #Std Error write standard error to this file\n#SBATCH --output=job.%J.out #Std Out write standard output to this file\n#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)\n#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line\n#SBATCH --partition=standard #Request a specific partition for the resource allocation.\n#This script generates VCF and mask files from individual diploid BAM files.\n\ndata=\"./../data/pseudodiploid\"\nout=\"./../data/pseudodiploid_vcf\"""\nsample=\"$i\"\nmsmctools=\"../../packages/msmc-tools/\"
    refGenome=\"../data/GCF_000219625.1_MYCGR_v2.0_genomic.fna\"\n#estimating the average sequencing depth using all site on chromosome 2 = NC_018217.1\n\nDEPTH=\$(samtools depth -r NC_018218.1 \$data/\$sample.pseudoDiploid.readGroup.bam | awk '{sum += \$3} END {print sum / NR}')\necho "\$"DEPTH"" > \$out/tmp\nfor CHR in NC_018218.1 NC_018217.1 NC_018216.1 NC_018215.1 NC_018214.1 NC_018213.1 NC_018212.1 NC_018211.1 NC_018210.1 NC_018209.1 NC_018208.1 NC_018207.1 NC_018206.1 NC_018205.1 NC_018204.1 NC_018203.1 NC_018202.1 NC_018201.1 NC_018200.1 NC_018199.1 NC_018198.1; do\nsamtools mpileup -B -q 20 -Q 20 -C 50 -g -r \$CHR -f \$refGenome \$data/\$sample.pseudoDiploid.readGroup.bam | bcftools call -c -V indels | \$msmctools/bamCaller.py \$DEPTH \$out/\$sample.mask.chr\$CHR.bed.gz | gzip -c > \$out/\$sample.chr\$CHR.vcf.gz;\ndone
    " > stp_02_call-mask-aeg_$i.sh;done


D. The last step will generate invidual scripts to peforM the SNVcalling and masking at once: 
script name: stp_02_call-mask-wht_$i.sh, where $i is each of the pseudodiploid genotypes.


E. Display the content with the command
	stp_02_call-mask-wht_$i.sh`

```
	#!/bin/bash
 
	#SBATCH --job-name=call-mask
	#SBATCH --nodes=1 
	#SBATCH --ntasks=1 
	#SBATCH --ntasks-per-node=1
	#SBATCH --cpus-per-task=1 
	#SBATCH --time=24:00:00 
	#SBATCH --mem=64G 
	#SBATCH --error=job.%J.err 
	#SBATCH --output=job.%J.out
	#SBATCH --mail-type=FAIL 
	#SBATCH --mail-user=rojas@evolbio.mpg.de 
	#SBATCH --partition=standard 

	#This script genarates VCF and mask files from individual diploid BAM files.

	data="./../data"
	out="./../out"
	bin="./../bin"
	sample="ZT559.ZT561"
	refGenome="GCF_000219625.1_MYCGR_v2.0_genomic.fna"

	#estimating the average sequencing depth using all site on chromosome 1 = NC_018218.1

	DEPTH=$(samtools depth -r NC_018218.1 $data/$sample.pseudoDiploid.readGroup.bam | awk '{sum += $3} END {print sum / NR}')

	echo "$DEPTH" > $out/tmp

	for CHR in NC_018218.1 NC_018217.1 NC_018216.1 NC_018215.1 NC_018214.1 NC_018213.1 NC_018212.1 NC_018211.1 NC_018210.1 NC_018209.1 NC_018208.1 NC_018207.1 NC_018206.1 NC_018205.1 NC_018204.1 NC_018203.1 NC_018202.1 NC_018201.1 NC_018200.1 NC_018199.1 NC_018198.1; do
	samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $CHR -f $data/$refGenome $data/$sample.pseudoDiploid.readGroup.bam | bcftools call -c -V indels | ./msmc-tools/bamCaller.py $DEPTH $out/$sample.mask.chr$CHR.bed.gz | gzip -c > $out/$sample.chr$CHR.vcf.gz
	done
```

F. Submit the individual files  generated as independent jobs.

	for i in *aeg.sh; do sbatch $i;done
    

## Step 3. Select combinations for msmc2 input file

msmc2 requires an input file tha includes genotypes from both divergent populations, hereinafter pop1 and pop2. Each individual input files can consist of 2 haplotypes (1 from each pop), 4 haplotypes (2 from each pop), or 16 haploptypes (8 from each pop).

For *Z. tritici* I only considered the 13 chromosomes from the core genome, and then created 10 input files based on the permutaions in table 1.


**Table 1. List of genotype permutations to generate input files for msmc2**

Permutation|	Aegilops cylindrica (pop1)                   |	Triticum (pop2)                                      |
-----------|--------------------------------------------|--------------------------------------------------|
1	|ZT442.ZT431, ZT454.ZT429, ZT460.ZT528, ZT441.ZT498	|ZT617.ZT576, ZT567.ZT649, ZT611.ZT638, ZT676.ZT583|
2	|ZT442.ZT429, ZT431.ZT536, ZT460.ZT537, ZT479.ZT440	|ZT709.ZT620, ZT705.ZT642, ZT617.ZT573, ZT644.ZT671|
3	|ZT442.ZT476, ZT448.ZT427, ZT536.ZT475, ZT460.ZT508	|ZT662.ZT642, ZT565.ZT655, ZT573.ZT705, ZT710.ZT709|
4	|ZT536.ZT485, ZT431.ZT460, ZT476.ZT479, ZT495.ZT517	|ZT611.ZT549, ZT583.ZT651, ZT709.ZT721, ZT555.ZT682|
5	|ZT536.ZT440, ZT485.ZT504, ZT500.ZT427, ZT517.ZT537	|ZT704.ZT643, ZT671.ZT661, ZT565.ZT709, ZT611.ZT710|
6	|ZT481.ZT508, ZT448.ZT427, ZT489.ZT436, ZT431.ZT460	|ZT709.ZT644, ZT721.ZT681, ZT640.ZT676, ZT587.ZT710|
7	|ZT481.ZT534, ZT500.ZT440, ZT454.ZT429, ZT441.ZT498	|ZT555.ZT668, ZT678.ZT643, ZT567.ZT555, ZT678.ZT699|
8	|ZT441.ZT498, ZT536.ZT475, ZT431.ZT460, ZT517.ZT537	|ZT559.ZT570, ZT642.ZT699, ZT576.ZT634, ZT655.ZT635|
9	|ZT481.ZT508, ZT534.ZT485, ZT476.ZT479, ZT500.ZT440	|ZT661.ZT649, ZT549.ZT651, ZT655.ZT643, ZT555.ZT570|
10	|ZT453.ZT442, ZT495.ZT517, ZT448.ZT427, ZT460.ZT537	|ZT559.ZT616, ZT583.ZT573, ZT555.ZT638, ZT681.ZT634|

The script below will generate an input file for ach of the 13 chromosomes of *Z. tritici*. I run each combination on an independet script.

A. To display the content of the script use

	cat 03_make_stp_03_input-msmc_16hap_aeg-trt_nonMskd.comb1.sh
This will show
```
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
SAMPLE5="ZT617.ZT576" #trt-pop2
SAMPLE6="ZT567.ZT649" #trt-pop2
SAMPLE7="ZT611.ZT638" #trt-pop2
SAMPLE8="ZT676.ZT583" #trt-pop2

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
```

## Step 4: Population size estimation
msmc 2 will generate one input file per chromosome, those can be used to compute three coalescence rate functions. Two within  Aegilops and Triticum infecting-populations, these are used to estimate effective population size for each pop;  and and a third one that is across populations and will be used in the next step to compute Cross Coalescence Rate (CCR).

Using a higher number of haplotypes (16 haps) allows to get a higher resolution when estimating Ne for recent times.


In the  script: stp_04_across-pop_ag-wh_16hap.sh I looped through the 10 haplotype combinations to compute the within and between coalescence rates.

To display the script execute the command:

	cat ./stp_04_across-pop_ag-wh_16hap.sh

This will shows:
```
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

#This script uses as input two pseudo-diploid files to runs msmc for the 21 chromosomes of Z. tritici

for i in 1 2 3 4 5 6 7 8 9 10; do

out="./../data/msmc_16hap"

msmc2 -I 0,1,2,3,4,5,6,7 -o $out/aeg.16hapl.allChr.within-pop1.comb$i.output $out/aeg-wht.16hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018218.1.msmc.input

msmc2 -I 8,9,10,11,12,13,14,15 -o $out/wht.16hapl.allChr.within-pop2.comb$i.output $out/aeg-wht.16hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018218.1.msmc.input

msmc2 -I 0-8,0-9,0-10,0-11,0-12,0-13,0-14,0-15,1-8,1-9,1-10,1-11,1-12,1-13,1-14,1-15,2-8,2-9,2-10,2-11,2-12,2-13,2-14,2-15,3-8,3-9,3-10,3-11,3-12,3-13,3-14,3-15,4-8,4-9,4-10,4-11,4-12,4-13,4-14,4-15,5-8,5-9,5-10,5-11,5-12,5-13,5-14,5-15,6-8,6-9,6-10,6-11,6-12,6-13,6-14,6-15,7-8,7-9,7-10,7-11,7-12,7-13,7-14,7-15 -o $out/aeg-wht.16hapl.allChr.across-pop12.comb$i.output $out/aeg-wht.16hapl.comb$i.chr.NC_018206.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018207.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018208.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018209.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018210.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018211.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018212.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018213.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018214.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018215.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018216.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018217.1.msmc.input $out/aeg-wht.16hapl.comb$i.chr.NC_018218.1.msmc.input ;
rm $out/aeg-wht.16hapl.comb$i.chr.NC_*.1.msmc.input;
done
```


## Step 5. Compute CCR

msmc2 provides the python script: combineCrossCoal.py to  compute  RCCR between pop1 and po2. I computed RCCCR for each of the 10 combinaations.
```
#!/bin/bash
 
#SBATCH --job-name=combine-coal #Give your job a name.
#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if comb to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 #Multithreading.
#SBATCH --time=48:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=64G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.

#This script uses as input two pseudo-diploid files to runs msmc for the 13 chromosomes of Z. tritici

msmctools="/data/biosoftware/msmc-tools/msmc-tools"
data="./../data"
out="./../data/msmc_16hap"

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb1.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb1.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb1.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb1.final.txt

$msmctools/combineCrossCoal.py $out/wht-aeg.16hapl.comb2.12seg.allChr.across-pop12.output.final.txt $out/wheat.16hapl.12seg.comb2.allChr.msmc.output.final.txt $out/aeg.16hapl.12seg.comb2.allChr.msmc.output.final.txt > $out/combined_wht-aeg_16hap_msmc.12seg.comb2.final.txt

$msmctools/combineCrossCoal.py $out/wht-aeg.16hapl.comb3.12seg.allChr.across-pop12.output.final.txt $out/wheat.16hapl.12seg.comb3.allChr.msmc.output.final.txt $out/aeg.16hapl.12seg.comb3.allChr.msmc.output.final.txt > $out/combined_wht-aeg_16hap_msmc.12seg.comb3.final.txt

$msmctools/combineCrossCoal.py $out/wht-aeg.16hapl.comb4.12seg.allChr.across-pop12.output.final.txt $out/wheat.16hapl.12seg.comb4.allChr.msmc.output.final.txt $out/aeg.16hapl.12seg.comb4.allChr.msmc.output.final.txt > $out/combined_wht-aeg_16hap_msmc.12seg.comb4.final.txt

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb5.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb5.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb5.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb5.final.txt 

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb6.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb6.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb6.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb6.final.txt 

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb7.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb7.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb7.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb7.final.txt

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb8.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb8.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb8.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb8.final.txt .

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb9.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb9.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb9.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb9.final.txt

$msmctools/combineCrossCoal.py $out/aeg-wht.16hapl.allChr.across-pop12.comb10.output.final.txt $out/aeg.16hapl.allChr.within-pop1.comb10.output.final.txt $out/wht.16hapl.allChr.within-pop2.comb10.output.final.txt > $out/combined_wht-aeg_16hap_msmc.comb10.final.txt
```

## Step 6. Plot RCCR and Ne

After computing Ne and RCCR these can be plot on R. To scale Ne and RCCR to years or number of  generations, it must be provided a mutation rate (u) per cell per generation and number of generations per year, if one generation per year g=1, if two generations per year g=0.5, etc.

The code to plot Ne and RCCR is provided in the R markdown file **stp_06_Ne_CCR.Rmd**



