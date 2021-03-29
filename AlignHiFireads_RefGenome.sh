#Align HiFi reads to Reference Genome of close species

#create an index of the reference genome (generates a .mmi file)
pbmm2 index /hpc-cloud/users/mathnat/ibb/dominguez/ReferenceGenomes/GCF_000277835.1_GeoFor_1.0_genomic.fasta /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.mmi

#align the HiFi reads to the index (generates a .bam file)
pbmm2 align /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.mmi /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/m64064_210131_082109.hifi.subreads.fasta /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.movie.bam

#align the HiFi reads to the reference genome using bam file
pbmm2 align /hpc-cloud/users/mathnat/ibb/dominguez/ReferenceGenomes/GCF_000277835.1_GeoFor_1.0_genomic.fasta /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/m64064_210131_082109.hifi.subreads.fasta /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.movie.bam --preset CCS --sort -j 0 -J 0


#To extract accuracy metrics  
samtools view /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.movie.bam | awk '{ mc=""; for(i=12;i<=NF;i++) { split($i,TAG,":"); if(TAG[1]=="mc") {mc=TAG[3]; break; } } if(mc != "") { print $1 "\t" mc; } }' > /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/MappedConcordance.HiFiReadSet.GenomeGfortis.out

#To extract read length metrics
samtools view /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.movie.bam | head -n 4763 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/MappedRL.HiFiReadSet.GenomeGfortis.out

#coverage metrics
samtools depth -a /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/ref_Gfortis.movie.bam > /mnt/scratch/mathnat/ibb/dominguez/YelCar_GenomeAssembly/Analysis/Assembly_QC/MapReadsRefGen/HiFiReadSet.REFGfortis.sorted.Depth.out
