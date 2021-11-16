
## bam to fq (using bedtools v.2.29.1) ##
module purge
module load Anaconda/3.7
source activate bedtools-env2

bedtools bamtofastq -i /path/filename.ccs.bam -fq  /path/reads/filename_ccs.fq

### READS STATISTICS ###

## Plot reads length ##
module purge
module load Anaconda/3.7
source activate biopython-env

python -m plot_fasta_length_v2.py /path/reads/filename_ccs.fq /path/reads_statistics/filename_ccs.Plotlength.png

## Reads statistics with asmstats ##
module purge
module load Anaconda/3.7
source activate seqtk-env

asmstats /path/filename.ccs.bam  > /path/reads_statistics/filename_ccs.fq.stats

#Alternative: reads stats with gaas
module purge
module load Anaconda/3.7
source activate gaas-env

sed -n '1~4s/^@/>/p;2~4p' /path/reads/filename_ccs.fq > /path/reads/filename_ccs_reads.fa

gaas_fasta_statistics.pl -f /path/reads/filename_ccs_reads.fa

## Count reads on fasta ##
grep '^>'/path/reads/filename_ccs_reads.fa | wc -l

## Contamination check ##
#https://kat.readthedocs.io/en/latest/
module purge
module load Anaconda/3.7
source activate kat-env2

kat gcp -o /path/reads_statistics/filename_ccs_fasta.kat -t 40 -p png /path/reads/filename_ccs_reads.fa 

## Kmer analysis ##
#k=21
module purge
module load Anaconda/3.7
source activate jellyfish-env

jellyfish count -C -m 21 -s 1000 -t 24 -o /path/reads_statistics/filename_ccs.k21.jf /path/reads/filename_ccs_reads.fa

#k31
#change -3 to 31

#open the .histo file created in genomescope2 website 

### GENOME ASSEMBLY ###

# 1. hicanu
# https://github.com/marbl/canu/releases
/path/canu-2.1.1/bin/canu -d /path/asm/ -p filename_hicanu genomeSize=1040m -pacbio-hifi /path1/reads/filename_ccs_reads.fa useGrid=false

# change genomeSize for what was predicted by genomescope2

#2. hifiasm
# https://github.com/chhylp123/hifiasm
module purge
module load Anaconda/3.7
source activate hifiasm-env

hifiasm -t60 -o /path/asm/filename.hifiasm /path/reads/filename_ccs_reads.fa

# .gfa file can be viewed with Bandage

#3. IPA
#https://github.com/PacificBiosciences/pbipa

module purge
module load Anaconda/3.7
source activate snakemake-env2
source activate ipa-env

ipa local --nthreads 60 --njobs 1 -i /path/reads/filename_ccs_reads.fa

### GENOME ASSEMBLY EVALUATION ###

#in construction

## Purge dups + minimap
## merqury
## blobtools
## busco
## mummer to compare with other assemblies
