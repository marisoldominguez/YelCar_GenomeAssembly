#Align HiFi reads to Reference Genome of close species

#create an index of the reference genome (generates a .mmi file)
pbmm2 index /path/refgenome.fasta /path/ref_Gfortis.mmi
#align the HiFi reads to the index (generates a .bam file)
pbmm2 align /path/ref_Gfortis.mmi /path/hifi.subreads.fasta /path/ref_Gfortis.movie.bam

#align the HiFi reads to the reference genome using bam file
pbmm2 align /path/refgenome.fasta /path/hifi.subreads.fasta /path/ref_Gfortis.movie.bam --preset CCS --sort -j 0 -J 0 --log-level INFO


#To extract accuracy metrics  
samtools view /path/ref_Gfortis.movie.bam | awk '{ mc=""; for(i=12;i<=NF;i++) { split($i,TAG,":"); if(TAG[1]=="mc") {mc=TAG[3]; break; } } if(mc != "") { print $1 "\t" mc; } }' > /path/MappedConcordance.HiFiReadSet.GenomeGfortis.out

#To extract read length metrics
samtools view /path/ref_Gfortis.movie.bam | head -n 4763 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > /path/MappedRL.HiFiReadSet.GenomeGfortis.out

#coverage metrics
samtools depth -a /path/ref_Gfortis.movie.bam > /path/HiFiReadSet.REFGfortis.sorted.Depth.out

