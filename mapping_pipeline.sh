#!/bin/sh

# Mapping pipeline for phage deep sequencing. Map to reference then look at diversity among treatments

# Using samtools mpileup for SNP calling

# Get reads from European Nucleotide Archive (EBI) accession: PRJEB25016

# Stick all files in same directory and make filenames more simple:
mkdir Trimmed/all_reads/

for f in Trimmed/Sample_*/*.gz; do mv "$f" Trimmed/all_reads/; done

# Rename to lose all extra barcodes etc.
rename -v 's/[A-Z]+-[A-Z]+_L001_//' Trimmed/all_reads/*.fastq.gz

# Merge sequences w/ FLASH:
mkdir Trimmed/all_reads/flashed/

for i in `ls Trimmed/all_reads/*_R1_001.fastq.gz`
do
  echo $i
  base=$(basename $i "R1_001.fastq.gz")
  echo $base
  flash 20170801_Meaden1/reads/${base}_1_trimmed.fastq.gz 20170801_Meaden1/reads/${base}_2_trimmed.fastq.gz > ${dir}/${base}.sam
  ~/programs/FLASH-1.2.11/flash --output-prefix $base -d Trimmed/all_reads/flashed/ --max-overlap 150 Trimmed/all_reads/${base}R1_001.fastq.gz Trimmed/all_reads/${base}R2_001.fastq.gz
done

# Map to reference w/ BWA-MEM. 

#ref=DMS3_NC_008717.1.fna
# Edited DMS3_NC_008717.1 genome to match Cady et al. 2012. doi:  10.1128/JB.01184-12
ref=DMS3_mVir_Cady.fna
bwa index $ref

# Map to ref:
mkdir alignments/
for i in `ls Trimmed/all_reads/flashed/*.extendedFrags.fastq`
do
  echo $i
  base=$(basename $i ".extendedFrags.fastq")
  echo $base
bwa mem -t 16 $ref ${i} > alignments/${base}.sam
done

# Do conversions and tidy up files.
# Edit. Don't dedup reads- pooled data!

for sample in alignments/*.sam
do
  echo $sample
# Convert SAM to BAM
  describer=$(echo ${sample} | sed 's/.sam//')
  echo $describer
    samtools view -bT $ref $sample > ${describer}.uns.bam
# Sort BAM file
    samtools sort ${describer}.uns.bam ${describer}
# Remove intermediates:
    rm ${describer}.uns.bam
# Index bam file
    samtools index ${describer}.bam
# !Don't! remove duplicate reads.
done


# Call mpileup for each sample individually
for sample in alignments/*.bam
  do
    echo $sample
    samtools mpileup -Q35 -A -f $ref -d 150000  $sample > ${sample}.raw.bcf
  done


# Parse with Awk to get altref read counts per site so can assess SNP frequency
for sample in alignments/*.raw.bcf
  do
    echo $sample
    base=$(basename $sample ".raw.bcf")
    awk 'BEGIN {FS="\t"; OFS="\t"} $1 ~ /NC_mVir_Cady/ {
       ref1=gsub(/[.,]/,"",$5);
       alt1=gsub(/[AGCTagct]/,"",$5);
       n1=gsub(/[Nn]/,"",$5);
       indel1=gsub(/[+-]/,"",$5);
       print $1,$2,$3,$4,ref1,alt1,n1,indel1}' $sample > alignments/${base}_pileup_out.txt
done

# Map spacer sequences to phage:. Use BIM5 spacers as BIM2 and WT are subsets of this
bwa mem $ref bim5.fna > target_seqs1.sam

# Get coordinates: Field 4 = left most aligned base. Length = 32.
# Zero based as SAM files start at 1.
awk 'BEGIN {OFS="\t"} $1 ~ /BIM5_PA14_spacer/ {pos1=$4-1;pos2=$4+33;print $1,pos1,pos2}'  target_seqs1.sam > target_seqs1.bed

# Finish analysis in R
