### PacBio HiFi reads pre-processing ###

## 1. bam to fq 
#a. (using bedtools v.2.29.1) ##
module purge
module load lang/Anaconda3/2020.11
source activate bedtools-env2

bedtools bamtofastq -i /path/filename.ccs.bam -fq  /path/reads/filename_ccs.fq

#b. bam to .fq and to .fasta (bam2fastx using SMRTtool) ##
#1. index
module purge
module load lang/Anaconda3/2020.11
source activate pbbam-env

pbindex /path/m64046_210816_140027.ccs.bam

#2. convert bam to fq
module load lang/Anaconda3/2020.11
source activate bam2fastx-env

bam2fasta -u -o /path/bam2fastx/ /path/m64046_210816_140027.ccs.bam

bam2fastq -u -o /path/bam2fastx/ /path/m64046_210816_140027.ccs.bam

conda deactivate

## 2. Trimming PacBio SMRT adapters from reads
cutadapt -j 10 -o /path/reads_noAdapt.fastq /path/reads.fastq -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT --discard-trimmed

#4.5%  of reads removed.

### Reads Statistics ###

#1. Fastqc on trimmed reads

fastqc -t 6 /path/reads_noAdapt.fastq -o /path/fastqc_TrimmedReads

#2. Plot reads length
python plot_fasta_length_v2.py /path/reads/reads_noAdapt.fasta /path/reads_statistics/reads_noAdapt.fasta.length.png

#3. Reads statistics
#a. Count reads on .fasta file:
grep ">" /path_to_reads/reads_noAdapt.fasta | wc -l
#N = 2 489 222

#b. Count reads on .fastq file:
expr $(cat /path_to_reads/reads_noAdapt.fastq | wc -l) / 4
#N = 2 489 222

#c. Using asmstats script:

module purge
module load lang/Anaconda3/2020.11
source activate seqtk-env

asmstats /path_to_reads/reads_noAdapt.fastq  > /path/reads_statistics/reads_noAdapt.fastq.stats

#4. Kmer analysis
# Assuming kmer size 21
module purge
module load lang/Anaconda3/2020.11
source activate jellyfish-env

#jellyfish count
jellyfish count -C -m 21 -s 1000000 -t 20 -o /path/reads_statistics/YelCar_reads_trim.k21.jf /path/reads/reads_noAdapt.fastq

#jellyfish histo
jellyfish histo /path/reads_statistics/YelCar_reads_trim.k21.jf > /path/reads_statistics/YelCar_reads_trim.k21.histo

#Assumin kmer size 31
module purge
module load lang/Anaconda3/2020.11
source activate jellyfish-env

#jellyfish count
jellyfish count -C -m 31 -s 1000000 -t 20 -o /path/reads_statistics/YelCar_reads_trim.k31.jf /path/reads/reads_noAdapt.fastq

#jellyfish histo
jellyfish histo /path/reads_statistics/YelCar_reads_trim.k31.jf > /path/reads_statistics/YelCar_reads_trim.k31.histo

#open the .histo file created in genomescope2 website: http://qb.cshl.edu/genomescope/genomescope2.0/

### GENOME ASSEMBLY ###
# 1a. Hifiasm with different parameters (l0,l1,l2 and nok, k21, k31)
# https://github.com/chhylp123/hifiasm

#scripts for assemblies (change parameters -l and -k accordingly)

module purge
module load lang/Anaconda3/2020.11
source activate hifiasm-env

hifiasm -t 28 -l2 -- primary -o /path/asm/hifiasm.trim.l2.k31.asm /path/reads/reads_noAdapt.fastq -k 31

# it took Real time: 27160.206 sec; CPU: 922008.523 sec; Peak RSS: 54.377 GB requesting 38 cpus-per-task with 16G mem-per-cpu.
# .gfa file can be viewed with Bandage

# 1b. script for statistics for the asm:

/path/gfa2fa /path/hifiasm.trim.l2.k31.asm.bp.p_ctg.gfa > /path/YelCar_hifiasm_trimmed_l2_k31.asm.p_ctg.fa

module purge
module load lang/Anaconda3/2020.11
source activate seqtk-env

export PATH=$PATH:/work/dominguez/YelCar_GenomeAssembly_MPI/

asmstats /path/YelCar_hifiasm_trimmed_l2_k31.asm.p_ctg.fa > /path/hifiasm.trim.l2.k31.asm.p_ctg.fa.stats

# 1c. First Round of purging (purge PRI p1 to obtain c1):
module purge
module load lang/Anaconda3/2020.11
source activate purge_dups_env

## Change this: ##########################################################
OUT_DIR="/path/hifiasm/purge_dups"
NAME="c1.l2.k31"
ASM="/path/YelCar_hifiasm_trimmed_l2_k31.asm.p_ctg.fa"
ASM_NAME=$(basename $ASM .fa)
HIFI_DIR="/path_to_reads"
TRANS=21
MAX_DEPTH=62
CUTOFFS=$TRANS\_$MAX_DEPTH
###############################################################################

if [ ! -d $OUT_DIR/$NAME ]
then
  mkdir -p $OUT_DIR/$NAME
fi

##
# Align HiFi reads and generate BAM files, then calculate read depth histogram and base-level read depth
cd $OUT_DIR/$NAME

for i in $(ls $HIFI_DIR/*_noAdapt.fastq.gz)
do
	READ_NAME=$(basename $i .fastq.gz)
	minimap2 -xasm20 -I 4G -t 24 $ASM $i | gzip -c - > $OUT_DIR/$NAME/$READ_NAME\.paf.gz
done
## -xmap-hifi is not in latest release. See https://github.com/lh3/minimap2/issues/739
##  we will use the Reads mapped to contigs pairwise mapping format (PAF) file for calculating some statistics required in a later stage
## purge_dups initially produces a read-depth histogram from base-level coverages. This information is used for estimating the coverage cutoffs, taking into account that collapsed haplotype contigs will lead to reads from both alleles mapping to those contigs, whereas if the alleles have assembled as separate contigs, then the reads will be split over the two contigs, resulting in half the read-depth (Roach et al. 2018).

pbcstat *.paf.gz 
# this produces PB.base.cov and PB.stat files
calcuts PB.stat > cutoffs_default 2>calcults_default.log
python3 /path/hist_plot.py -c cutoffs_default PB.stat PB_default.png

### VGP suggested cutoffs:
calcuts -m$TRANS -u$MAX_DEPTH PB.stat > cutoffs_$CUTOFFS
python3 /path/hist_plot.py -c cutoffs_$CUTOFFS PB.stat PB_$CUTOFFS\.png

##
# Split an assembly and do a self-self alignment
cd $OUT_DIR/$NAME
split_fa $ASM > $ASM_NAME\.split
minimap2 -t 18 -xasm5 -DP $ASM_NAME\.split $ASM_NAME\.split | pigz -p 6 -c - > $ASM_NAME\.split.self.paf.gz

##### 

# Purge haplotigs and overlaps:
cd $OUT_DIR/$NAME
purge_dups -2 -T cutoffs_default -c PB.base.cov $ASM_NAME\.split.self.paf.gz > dups_default.bed 2> purge_dups_default.log
purge_dups -2 -T cutoffs_$CUTOFFS -c PB.base.cov $ASM_NAME\.split.self.paf.gz > dups_$CUTOFFS\.bed 2> purge_dups_$CUTOFFS\.log

##
# Get purged primary and haplotig sequences from draft assembly:
cd $OUT_DIR/$NAME
get_seqs -e dups_default.bed $ASM
mv purged.fa purged_default.fa
mv hap.fa hap_default.fa
get_seqs -e dups_$CUTOFFS\.bed $ASM
mv purged.fa purged_$CUTOFFS\.fa
mv hap.fa hap_$CUTOFFS\.fa

############################################################################
module purge
module load lang/Anaconda3/2020.11
source activate seqtk-env

export PATH=$PATH:/work/dominguez/YelCar_GenomeAssembly_MPI/

asmstats /path/hifiasm/purge_dups/c1.l2.k31/hap_21_62.fa > /path/hifiasm/purge_dups/c1.l2.k31/k31_hap_21_62.fa.stats
asmstats /path/hifiasm/purge_dups/c1.l2.k31/purged_21_62.fa > /path/hifiasm/purge_dups/c1.l2.k31/k31_purged_21_62.fa.stats
asmstats /path/hifiasm/purge_dups/c1.l2.k31/hap_default.fa > /path/hifiasm/purge_dups/c1.l2.k31/k31_hap_default.fa.stats
asmstats /path/hifiasm/purge_dups/c1.l2.k31/purged_default.fa > /path/hifiasm/purge_dups/c1.l2.k31/k31_purged_default.fa.stats
########################################################################

# "cutoffs" file has six numbers: the first is low coverage, any contig with an average coverge less than the number is a junk contig, the forth number is diploid (middle) coverage, the last number is high coverage for repeats, contigs with an average coverage over this number are usually repeats 

# Next step: merge hap.fa and alt_asm and redo the above steps to get the haplotig set.

#job was ran with --cpus-per-task=24, #SBATCH --mem=64G. It took 15 min.


# 1d. Second Round of purging (purge ALT c2p2 to obtain q2):
#using VGP suggested cutoffs:

module purge
module load lang/Anaconda3/2020.11
source activate purge_dups_env

## Change this: ##########################################################
OUT_DIR="/path/hifiasm/purge_dups"
NAME="c1.l2.k31"
ASM="/path/YelCar_hifiasm_trimmed_l2_k31.asm.p_ctg.fa"
ASM_NAME=$(basename $ASM .fa)
ALT="/path/YelCar_hifiasm_trimmed_l2_k31.asm.a_ctg.fa"
HIFI_DIR="/path_to_reads"
TRANS=21
MAX_DEPTH=62
CUTOFFS=$TRANS\_$MAX_DEPTH
###############################################################################

if [ ! -d $OUT_DIR/$NAME/q2_VGPcutoffs ]
then
  mkdir -p $OUT_DIR/$NAME/q2_VGPcutoffs
fi

#Merge hap_.fa and ALT 
##c2 is before purging so it is in another folder
cd $OUT_DIR/$NAME/q2_VGPcutoffs
cat $ALT >> c2p2.fa; cat /path/hifiasm/purge_dups/c1.l2.k31/hap_21_62.fa >> c2p2.fa

#align HiFi reads to c2p2 and generate BAM files, then calculate read depth histogram and base-level read depth
cd $OUT_DIR/$NAME/q2_VGPcutoffs
for i in $(ls $HIFI_DIR/*_noAdapt.fastq.gz)
do
        READ_NAME=$(basename $i .fastq.gz)
        minimap2 -xasm20 -I 4G -t 24 c2p2.fa $i | gzip -c - > $OUT_DIR/$NAME/q2_VGPcutoffs/$READ_NAME\.paf.gz
done
pbcstat *.paf.gz

### VGP suggested cutoffs:
calcuts -m$TRANS -u$MAX_DEPTH PB.stat > cutoffs_$CUTOFFS
python3 /path/hist_plot.py -c cutoffs_$CUTOFFS PB.stat PB_$CUTOFFS\.png
##

# Split an assembly and do a self-self alignment
cd $OUT_DIR/$NAME/q2_VGPcutoffs
split_fa c2p2.fa > c2p2.split
minimap2 -t 18 -xasm5 -DP c2p2.split c2p2.split | pigz -p 6 -c - > c2p2.split.self.paf.gz

#############
# Purge haplotigs and overlaps
cd $OUT_DIR/$NAME/q2_VGPcutoffs
purge_dups -2 -T cutoffs_$CUTOFFS -c PB.base.cov c2p2.split.self.paf.gz > dups_$CUTOFFS\.bed 2> purge_dups_$CUTOFFS\.log


# Get purged and haplotig sequences from draft assembly
cd $OUT_DIR/$NAME/q2_VGPcutoffs
get_seqs -e dups_$CUTOFFS\.bed c2p2.fa
mv purged.fa purged_$CUTOFFS\.fa
mv hap.fa hap_$CUTOFFS\.fa

#########################################################################
module purge
module load lang/Anaconda3/2020.11
source activate seqtk-env

export PATH=$PATH:/work/dominguez/YelCar_GenomeAssembly_MPI/

asmstats $OUT_DIR/$NAME/q2_VGPcutoffs/purged_$CUTOFFS\.fa > $OUT_DIR/$NAME/q2_VGPcutoffs/purged_$CUTOFFS\.fa.stats
asmstats $OUT_DIR/$NAME/q2_VGPcutoffs/hap_$CUTOFFS\.fa > $OUT_DIR/$NAME/q2_VGPcutoffs/hap_$CUTOFFS\.fa.stats

## requested --cpus-per-task=24, --mem=64G, and it took 15 min to run

#using default values:

module purge
module load lang/Anaconda3/2020.11
source activate purge_dups_env

## change this: ##########################################################
OUT_DIR="/path/hifiasm/purge_dups"
NAME="c1.l2.k31"
ASM="/path/YelCar_hifiasm_trimmed_l2_k31.asm.p_ctg.fa"
ASM_NAME=$(basename $ASM .fa)
ALT="/path/YelCar_hifiasm_trimmed_l2_k31.asm.a_ctg.fa"
HIFI_DIR="/path_to_reads"
###############################################################################
if [ ! -d $OUT_DIR/$NAME/q2_default ]
then
  mkdir -p $OUT_DIR/$NAME/q2_default
fi

#Merge hap_.fa and ALT 
##c2 is before purging so it is in another folder
cd $OUT_DIR/$NAME/q2_default
cat $ALT >> c2p2.fa; cat /hifiasm/purge_dups/c1.l2.k31/hap_default.fa >> c2p2.fa

# align HiFi reads to c2p2 and generate BAM files, then calculate read depth histogram and base-level read depth
cd $OUT_DIR/$NAME/q2_default
for i in $(ls $HIFI_DIR/*_noAdapt.fastq.gz)
do
        READ_NAME=$(basename $i .fastq.gz)
        minimap2 -xasm20 -I 4G -t 24 c2p2.fa $i | gzip -c - > $OUT_DIR/$NAME/q2_default/$READ_NAME\.paf.gz
done
pbcstat *.paf.gz
calcuts PB.stat > cutoffs_default 2>calcults_default.log
python3 /path/hist_plot.py -c cutoffs_default PB.stat PB_default.png

# Split an assembly and do a self-self alignment
cd $OUT_DIR/$NAME/q2_default
split_fa c2p2.fa > c2p2.split
minimap2 -t 18 -xasm5 -DP c2p2.split c2p2.split | pigz -p 6 -c - > c2p2.split.self.paf.gz

# Purge haplotigs and overlaps
cd $OUT_DIR/$NAME/q2_default
purge_dups -2 -T cutoffs_default -c PB.base.cov c2p2.split.self.paf.gz > dups_default.bed 2> purge_dups_default.log


# Get purged and haplotig sequences from draft assembly
cd $OUT_DIR/$NAME/q2_default
get_seqs -e dups_default.bed c2p2.fa
mv purged.fa purged_default.fa
mv hap.fa hap_default.fa
#########################################################################
module purge
module load lang/Anaconda3/2020.11
source activate seqtk-env

export PATH=$PATH:/work/dominguez/YelCar_GenomeAssembly_MPI/

asmstats $OUT_DIR/$NAME/q2_default/purged_default.fa > $OUT_DIR/$NAME/q2_default/purged_default.fa.stats
asmstats $OUT_DIR/$NAME/q2_default/hap_default.fa > $OUT_DIR/$NAME/q2_default/hap_default.fa.stats

## requested --cpus-per-task=24, --mem=64G, and it took 15 min to run

# 2. IPA (https://github.com/PacificBiosciences/pbipa)
#script for the assembly of the primmary (PRI) and alternative (ALT) assemblies:

module purge
module load lang/Anaconda3/2020.11
source activate snakemake-env2
source activate ipa-env

ipa local --nthreads 60 --njobs 1 -i /path/reads/filename_ccs_reads.fa

# It was ran requesting 50G of mem and lasted 2 days to finish.

# First round of PURGING (purge PRI p1 to obtain c1):
module purge
module load lang/Anaconda3/2020.11
source activate purge_dups_env

## Change this: ##########################################################
OUT_DIR="/path/IPA/purge_dups"
NAME="IPA_c1.purged"
ASM="/path/IPA/final.p_ctg.fasta"
ASM_NAME=$(basename $ASM .fa)
HIFI_DIR="/path_to_reads"
###############################################################################

if [ ! -d $OUT_DIR/$NAME ]
then
  mkdir -p $OUT_DIR/$NAME
fi

# Align HiFi reads and generate BAM files, then calculate read depth histogram and base-level read depth
cd $OUT_DIR/$NAME

for i in $(ls $HIFI_DIR/*_noAdapt.fastq.gz)
do
	READ_NAME=$(basename $i .fastq.gz)
	minimap2 -xasm20 -I 4G -t 24 $ASM $i | gzip -c - > $OUT_DIR/$NAME/$READ_NAME\.paf.gz
done
## -xmap-hifi is not in latest release. See https://github.com/lh3/minimap2/issues/739
## we will use the Reads mapped to contigs pairwise mapping format (PAF) file for calculating some statistics required in a later stage
## purge_dups initially produces a read-depth histogram from base-level coverages. This information is used for estimating the coverage cutoffs, taking into account that collapsed haplotype contigs will lead to reads from both alleles mapping to those contigs, whereas if the alleles have assembled as separate contigs, then the reads will be split over the two contigs, resulting in half the read-depth (Roach et al. 2018).

pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts PB.stat > cutoffs_default 2>calcults_default.log
python3 /path/hist_plot.py -c cutoffs_default PB.stat PB_default.png

# Split an assembly and do a self-self alignment:
cd $OUT_DIR/$NAME
split_fa $ASM > $ASM_NAME\.split
minimap2 -t 18 -xasm5 -DP $ASM_NAME\.split $ASM_NAME\.split | pigz -p 6 -c - > $ASM_NAME\.split.self.paf.gz

# Purge haplotigs and overlaps:
cd $OUT_DIR/$NAME
purge_dups -2 -T cutoffs_default -c PB.base.cov $ASM_NAME\.split.self.paf.gz > dups_default.bed 2> purge_dups_default.log

# Get purged primary and haplotig sequences from draft assembly:
cd $OUT_DIR/$NAME
get_seqs -e dups_default.bed $ASM
mv purged.fa purged_default.fa
mv hap.fa hap_default.fa

############################################################################
module purge
module load lang/Anaconda3/2020.11

source activate seqtk-env

export PATH=$PATH:/work/dominguez/YelCar_GenomeAssembly_MPI/

asmstats $OUT_DIR/$NAME/hap_default.fa > $OUT_DIR/$NAME/IPA_hap_default.fa.stats
asmstats $OUT_DIR/$NAME/purged_default.fa > $OUT_DIR/$NAME/IPA_purged_default.fa.stats

#Next step: Merge hap.fa and alt_asm and redo the above steps to get the haplotig set.

# It was ran requesting --cpus-per-task=24, --mem=64G. It took 10 min.

# Second round of PURGING (purge ALT c2p2 to obtain q2):
module purge
module load lang/Anaconda3/2020.11
source activate purge_dups_env

## Change this: ##########################################################
OUT_DIR="/path/reads_trimmed/IPA/purge_dups/"
NAME="IPA_c1.purged"
ASM="/path_to_PRIasm/final.p_ctg.fasta"
ASM_NAME=$(basename $ASM .fa)
ALT="/path_to_ALTasm/final.a_ctg.fasta"
HIFI_DIR="/path_to_reads"
###############################################################################

if [ ! -d $OUT_DIR/$NAME/q2 ]
then
  mkdir -p $OUT_DIR/$NAME/q2
fi

# Merge hap_.fa and ALT 
##c2 is before purging so it is in another folder
cd $OUT_DIR/$NAME/q2
cat $ALT >> c2p2.fa; cat /path/IPA/purge_dups/IPA_c1.purged/hap_default.fa >> c2p2.fa

# Align HiFi reads to c2p2 and generate BAM files, then calculate read depth histogram and base-level read depth:
cd $OUT_DIR/$NAME/q2
for i in $(ls $HIFI_DIR/*_noAdapt.fastq.gz)
do
        READ_NAME=$(basename $i .fastq.gz)
        minimap2 -xasm20 -I 4G -t 24 c2p2.fa $i | gzip -c - > $OUT_DIR/$NAME/q2/$READ_NAME\.paf.gz
done

pbcstat *.paf.gz

calcuts PB.stat > cutoffs_default 2>calcults_default.log
python3 /path/hist_plot.py -c cutoffs_default PB.stat PB_default.png

# Split an assembly and do a self-self alignment:
cd $OUT_DIR/$NAME/q2
split_fa c2p2.fa > c2p2.split
minimap2 -t 18 -xasm5 -DP c2p2.split c2p2.split | pigz -p 6 -c - > c2p2.split.self.paf.gz

# Purge haplotigs and overlaps:
cd $OUT_DIR/$NAME/q2
purge_dups -2 -T cutoffs_default -c PB.base.cov c2p2.split.self.paf.gz > dups_default.bed 2> purge_dups_default.log

# Get purged and haplotig sequences from draft assembly:
cd $OUT_DIR/$NAME/q2
get_seqs -e dups_default.bed c2p2.fa
mv purged.fa purged_default.fa
mv hap.fa hap_default.fa
#########################################################################

module purge
module load lang/Anaconda3/2020.11

source activate seqtk-env

export PATH=$PATH:/work/dominguez/YelCar_GenomeAssembly_MPI/

asmstats $OUT_DIR/$NAME/q2/purged_default.fa > $OUT_DIR/$NAME/q2/purged_default.fa.stats
asmstats $OUT_DIR/$NAME/q2/hap_default.fa > $OUT_DIR/$NAME/q2/hap_default.fa.stats

# It was ran with --cpus-per-task=24, #SBATCH --mem=24G and took 10 min.

### GENOME ASSEMBLY EVALUATION ###
##    GEP    ##

## purge dups


## GEP (James Sullivan - BeGenDiv- pipeline) to compare asms before and after purging

# following https://git.imp.fu-berlin.de/cmazzoni/GEP/-/tree/master
# 1. I git cloned the workflow 2. created a GEP conda environment (with snakemake) 3. created a sample sheet (.tsv) with the info of where are the reads, etc : /home/dominguez/GEP/YelCar_hifi_sampleSheet.tsv 4. changed the config.yaml to describe to Results file (can exist or it will be created) and the busco lineage*.Content of the config.yaml
Results:               "/work/dominguez/YelCar_GenomeAssembly_MPI/asm/reads_trimmed/GEP_Results/20032022"

samplesTSV:             "/home/dominguez/GEP/YelCar_hifi_sampleSheet.tsv"

busco5Lineage:          "vertebrata"

buscoMode:              "genome"

#* the name of the BUSCO lineage is choosen from https://busco.ezlab.org/list_of_lineages.html I chose vertebrata as that is what VGP does.

#modify /home/dominguez/GEP/SUBMIT_CONFIG/slurm/cluster.yaml to change the name of the partition to "all", there I can also change if I want the --mem and --time of each subjob, I only changed the partition name and left all the rest as default.

#for STEP 1: building the meryl kmer database
module load lang/Anaconda3
source activate GEP

nohup snakemake --profile SUBMIT_CONFIG/slurm --cores 20 &

#4 of 7 steps (57%) done... but those are the important steps anyway. fastqc did not worked!

#for STEP 2: run the genome asm evaluation
#first, make sure the old execution is finished or cancel it with scancel jobID. Also check in htop -u dominguez to see if any dormant processes show up that have the word snakemake or fastqc in them

#create a genome evaluation.tsv and change the config.yaml to put the path and name of this file instead of the file with the reads

#the file is C:\Users\Marisol\Documents\YeyCarGenomeAssembly\MPI-Dresden\Scripts_newHPC_Potsdam\runEval_YelCar_hifiasm_IPA.tsv

#the modified config.yaml:

#Results:               "/work/dominguez/YelCar_GenomeAssembly_MPI/asm/reads_trimmed/GEP_Results/200$

#samplesTSV:             "/home/dominguez/GEP/runEval_YelCar_hifiasm_IPA.tsv"

#busco5Lineage:          "vertebrata"

#buscoMode:              "genome"

module load lang/Anaconda3
source activate GEP
nohup snakemake --profile SUBMIT_CONFIG/slurm --cores 62 &


# it run it all, except from BUSCO
#17.03.2022 (see slack with James)

#Replace file GEP/envs/busco_and_assembly.yaml with the ine James sent by slack, and modify GEP/rules/run_05_11_21.smk , at line 114 replace the busco command with: (busco -m {params.mode} --offline --in {input.assembly} -o {params.assemblyName} -l {input.lineage} -c {threads} -f --limit 5) &> {log}
#(the only different is that the --offline and {params.mode} need to be switched places)

#then run in dry mode (snakemake -n) to see which steps it will repeat (if the pipeline recognizes the files that were already generated it will skip this steps and only do busco and reports).

#dryrun:
snakemake -n

#to see the table at the beginning of the output which jobs will be done, there you can see if it will run busco (plus the periphery steps to aggregate the full results)

#it will need to re-create a new conda environment since this was modified (it does this during the dryrun

#run pipeline again with more cores because it each busco job will take 16 threads
nohup snakemake --profile SUBMIT_CONFIG/slurm --cores 80 --jobs 5 &

#I put 80 cores (16x5) since each busco job will use 16 cores.

#18.03.2022

# finished running. It ran busco now, but did not produce the FinalReport because it could not create the PDF.
# the Final report should be in: 1_evaluation/finalResults

# Also, I checked /work/dominguez/YelCar_GenomeAssembly_MPI/asm/reads_trimmed/GEP_Results/20032022/1_evaluation/hifiasm_l1_k31/logs/hifiasm_l1_k31_busco5.log and looked good, but in the nohup.out I saw that 2 buscos failed, two markdown files missing,

#checked that there are no processes running 
ps -u dominguez
#there is a snakemake process still running
kill -9 42303
#check if it was killed
ps -u dominguez

# Erase and re-installed environmet with the new file James sent me: installGEP.yaml
conda env list
rm -r /home/dominguez/.conda/envs/GEP

module load lang/Anaconda3
source activate GEP
conda env create -f installGEP.yaml

#Run again (check before with dryrun which steps it still needs to do)
source activate GEP
snakemake -n
snakemake --unlock

nohup snakemake --profile SUBMIT_CONFIG/slurm --cores 5 --jobs 2 &
#or if I want to see what is happening:

snakemake --profile SUBMIT_CONFIG/slurm --cores 5 --jobs 2

#for the problem of pdf generation:
#modify FullMarkdown.md to FullMarkdown_modif.md (erase everything after \pagebreak (at the end))
tectonic
pandoc -o test.pdf FullMarkdown_modif.md --pdf-engine=tectonic


## all this will be fixed, so next time if I installed pipeline again , the PDF will be automatically run, and busco will be run with all the other programs.

#apart from the PDF with graphs for each assembly, there a table with key results, in .md and .tsv(this I can import in Excel):
/work/dominguez/YelCar_GenomeAssembly_MPI/asm/reads_trimmed/GEP_Results/20032022/1_evaluation/finalResults/FullTableMarkdown.md

#I opened it with notepad, copy and paste in a new HackMD note to see the table (because it is in markdown format)
https://hackmd.io/qHE-xrchRoWe6Ro0fPFj2Q?view

#achtung! the kmer multiplicity PRI (Stacked) is actually the flattened pri+alt (it is better to have ther .st but it is just a matter of taste), but those are in the individual merqury results folders
#ex: /work/dominguez/YelCar_GenomeAssembly_MPI/asm/reads_trimmed/GEP_Results/20032022/1_evaluation/hifiasm_l0_nok/04_merquryQVandKAT/hifiasm_l0_nok_merqOutput.hifiasm_l0_nok.PRI.spectra-cn.st.png


########################################

# Other assemblers I tried: Hicanu
# https://github.com/marbl/canu/releases

module purge
module load lang/Anaconda3/2020.11
source activate java11-env
export PATH=/work/dominguez/YelCar_GenomeAssembly_MPI/canu-2.1.1/bin:$PATH

canu -assemble -d /path/hicanu_FatNode/ -p YelCar_hicanu_assemble genomeSize=1.06g -pacbio-hifi /path_to_reads/reads_noAdapt.fastq useGrid=false

#change genomeSize for what was predicted by genomescope2
#Requested --mem=150G, --cpus-per-task=8, and took xx to finish.

###
# Other assemblers I tried: LJA
#https://github.com/AntonBankevich/LJA

module purge
module load lib/zlib/1.2.11-GCCcore-10.2.0
module load lang/Anaconda3
source activate cmake-env


cd /work/dominguez/YelCar_GenomeAssembly_MPI/asm/reads_trimmed/LJA/
/work/dominguez/LJA/bin/lja -o /path/LJA_fatnode_Cintel/ --reads /path_to_reads/reads_noAdapt.fastq --threads 8 --diploid

 #the job to assemble the hifi reads seems not to use much memory over most of its runtime. It seems that it is only at the end that it performs a huge allocation. And it is killed because of that.  So, I then ran it in our "fat" node (where I can request --mem=1500G). In this case, the job failed after 14 days of running due to a time limitation. I, then, tested it with a shorter version of the infile, requesting 24 CPUs (--mem-per-cpu=16G). It took 1h and then also failed in the step of uncompressing homopolymers in contigs. It produced the expected output files except from the final fasta assembly and final gfa. 
  
###
# Other assemblers: Rust-mdbg
#https://github.com/ekimb/rust-mdbg
#I did not try it given that it does not resolve both haplotypes.
###

#NETX STEPS:
#1. assembly mitochondrial genome
#2. Repeats annotation
#3. Genes Annotation
#4. CAFE analysis
