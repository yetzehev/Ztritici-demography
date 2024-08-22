## Demography of two divergent populations of *Triticum*- and *Aegilops*-infecting  *Z. tritici* 

## Pipeline description

###  Aim
To assess the  demographic history of *Aegilops*- and *Triticum*-infecting *Z. tritici* populations

### Data sets

A total number of 74 genotypes:
- *Triticum*-infecting  infecting Z.tritici genotypes (n = 45)
- *Aegilops cylindrica*-infecting Z. tritici genotypes (n= 28)

Note: Outiers isolates that infect *A. tauschii* and high IBS *Triticum*-infecting isolates were not included.
    

### Step 1: Perform the SNP calling and depth-masking with bcftools

Although *Zymoseptoria tritici* has a haploid genome, recombination still occurs during its sexual reproduction. However, the msmcm2 pipeline requires a diploid, phased VCF file as input. To address this, SNP calling was performed using the msmc2 pipeline with haploid data, but with the ploidy flag set to 2. This approach generated a diploid, homozygous VCF file, which was then used to create the MSMC2 input file. For the cross-coalescence computation, only one genome index was considered (see `stp_03_haploid2x_across-pop_ag-wh_08hap_r015.sh`).

Note: Accessory chromosomes were excluded from msmc2 demographic analysis.

Individual scripts were generated for each genotype and SNP calling and depth-masking were run in parallel. A sample script is showed bellow

	stp_01_call-mask-trt_ZT549.haploid2x.sh

E. Display the content with the command
``cat stp_01_call-mask-trt_ZT549.haploid2x.sh``

```
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
```

### Step 2.  Genotype combinations for msmc2 input file

msmc2 requires an input file tha includes genotypes from both divergent populations, hereinafter pop1 and pop2. Each individual input files can consist of 2 haplotypes (1 from each pop), 4 haplotypes (2 from each pop), or 16 haploptypes (8 from each pop).

**Table 1. List of genotype permutations to generate input files for msmc2**

Permutation|	Aegilops cylindrica (pop1)                   |	Triticum (pop2)                                      
-----------|--------------------------------------------|--------------------------------------------------
1	| ZT495, ZT454, ZT448, ZT440	| ZT617, ZT567, ZT611, ZT678 
2	| ZT460, ZT528, ZT441, ZT498	| ZT611, ZT638, ZT676, ZT583
3	| ZT442, ZT429, ZT431, ZT536	| ZT617, ZT573, ZT644, ZT671
4	| ZT460, ZT537, ZT479, ZT440	| ZT709, ZT620, ZT705, ZT642
5	| ZT442, ZT476, ZT448, ZT427	| ZT573, ZT705, ZT710, ZT709
6	| ZT536, ZT475, ZT460, ZT508	| ZT662, ZT642, ZT565, ZT655
7	| ZT536, ZT485, ZT431, ZT460	| ZT709, ZT721, ZT555, ZT682
8	| ZT476, ZT479, ZT495, ZT517	| ZT611, ZT549, ZT583, ZT651
9	| ZT536, ZT440, ZT485, ZT504	| ZT565, ZT709, ZT611, ZT710
10	| ZT500, ZT427, ZT517, ZT537	| ZT704, ZT643, ZT671, ZT661
11	| ZT481, ZT508, ZT448, ZT427	| ZT640, ZT676, ZT587 ZT710
12	| ZT489, ZT436, ZT431, ZT460	| ZT709, ZT644, ZT721, ZT681
13	| ZT481, ZT534, ZT500, ZT440	| ZT567, ZT555, ZT678, ZT699
14	| ZT454, ZT429, ZT441, ZT498	| ZT555, ZT668, ZT678, ZT643
15	| ZT441, ZT498, ZT536, ZT475	| ZT576, ZT634, ZT655, ZT635
16	| ZT431, ZT460, ZT517, ZT537	| ZT559, ZT570, ZT642, ZT699
17	| ZT481, ZT508, ZT534, ZT485	| ZT655, ZT643, ZT555, ZT570
18	| ZT476, ZT479, ZT500, ZT440	| ZT661, ZT649, ZT549, ZT651
19	| ZT453, ZT442, ZT495, ZT517	| ZT555, ZT638, ZT681, ZT634
20	| ZT448, ZT427, ZT460, ZT537	| ZT559, ZT616, ZT583, ZT573


The script below generates an msmc2 input file for each of the 13 chromosomes in the *Z. tritici* core genome. Each combination from Table 1 was run in a separate script.

A. To display the content of the script use

	cat 02_haploid2x_make_stp_02_input-msmc_08hap_aeg-trt_nonMskd.comb1.sh
This will show:
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

out="./../data/haploid2x"
msmcinput="./../data/msmc_haploid2x_08hap"

SAMPLE1=ZT495 #aeg-pop1
SAMPLE2=ZT454 #aeg-pop1
SAMPLE3=ZT448 #aeg-pop1
SAMPLE4=ZT440 #aeg-pop1
SAMPLE5=ZT617 #trt-pop2
SAMPLE6=ZT567 #trt-pop2
SAMPLE7=ZT611 #trt-pop2
SAMPLE8=ZT678 #trt-pop2
##N. combination
combN=comb1

#This script generates the masking files for 8 haploid2x bamfiles created on step 1; four for Aegilops and 4 for Triticum
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
$out/$SAMPLE1.chr$CHR.vcf.gz $out/$SAMPLE2.chr$CHR.vcf.gz $out/$SAMPLE3.chr$CHR.vcf.gz $out/$SAMPLE4.chr$CHR.vcf.gz $out/$SAMPLE5.chr$CHR.vcf.gz $out/$SAMPLE6.chr$CHR.vcf.gz $out/$SAMPLE7.chr$CHR.vcf.gz $out/$SAMPLE8.chr$CHR.vcf.gz> $msmcinput/aeg-wht.08hapl.$combN.chr.$CHR.msmc.input;done

```

### Step 3: Combine all the chromosomes per population and across populations
After creating one input file per chromosome, these will be combined  for computing the effective population size  per population and the  Relative Cross Coalescence Rate (CCR)  between populations.

The following script generated the  chromosome combination for the  20 haplotype combinations  described  in Table 1. **Important:** only one genome index was used per homozygous-diploid file. 

	stp_03_haploid2x_across-pop_ag-wh_08hap_r015.sh

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
```

## Step 4. Compute RCCR

msmc2 provides the python script: combineCrossCoal.py to  compute the relative cross coalescence rate (RCCR) between pop1 and pop2.  RCCR was computed  for each of the 20  haplotype combinations.

``cat stp_04_haploid2x_comb-cross_ag-wh_08hap_comb1-20.sh``



```
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
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do 
$msmctools/combineCrossCoal.py $out/aeg-wht.08hapl.allChr.across-pop12.comb$i.output.final.txt $out/aeg.08hapl.allChr.within-pop1.comb$i.output.final.txt $out/wht.08hapl.allChr.within-pop2.comb$i.output.final.txt > $out/combined_wht-aeg_08hapl_msmc.comb$i.r015.final.txt;
done
```

## Step 6. Plot RCCR and Ne

After computing Ne and RCCR these can be plot on R. To scale Ne and RCCR to years or number of  generations, it must be provided a mutation rate (u) per cell per generation and number of generations per year, if one generation per year g=1, if two generations per year g=0.5, etc.

The code to plot Ne and RCCR is provided in the R markdown file **stp_06_Ne_CCR.Rmd**