#!/bin/sh
## Ancient data processing and mapping using aln
## $1 - Threads
## $2 - Raw reads folder
## $3 - Code name for sample
## $4 - Results folder
## $5 - Reference
## $6 - Minimum read length
## $7 - Mismatch parameter

### R1 adapter trimmer and remove reads shorter than 30 and PCR duplicates
~/Software/fastp-0.23.2/bin/fastp --dedup -g -j $4/$3 -h $4/$3 -i $2/$3.damaged_s1.fq.gz -I $2/$3.damaged_s2.fq.gz -V -w $1 -q 25 -l 30 --merged_out $4/$3-merged.fastq -m

### Index the reference (prepare for mapping)
#bwa index $5

### Map reads R1 with BWA, parse output through SAMtools (index ref file first)
bwa aln -l 999 -o 2 -n $7 -t $1 $5 $4/$3-merged.fastq | bwa samse $5 - $4/$3-merged.fastq | samtools view -F 4 -q 20 -@ $1 -uS - | samtools sort -@ $1 -o $4/$3.sort.bam

## Remove duplicates sort and index bam
samtools rmdup -S $4/$3.sort.bam $4/$3.rmdup.bam
samtools sort -@ $1 $4/$3.rmdup.bam -o $4/$3.rmdup.sort.bam
samtools index $4/$3.rmdup.sort.bam

# Add read groups
picard AddOrReplaceReadGroups I=$4/$3.rmdup.sort.bam O=$4/$3.rmdup.sort_RG.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$3
samtools index $4/$3.rmdup.sort_RG.bam

## Check number of reads mapping before and after duplicate removal and depth

rm $4/$3_records.txt
## how many reads originally
zcat $2/$3.damaged_s1.fq.gz | echo $((`wc -l`/4)) >> $4/$3_records.txt
## how many reads post adapter and short read removal
cat $4/$3-merged.fastq | echo $((`wc -l`/4)) >> $4/$3_records.txt
## how many reads mapped pre duplicate removal
samtools view -c $4/$3.sort.bam >> $4/$3_records.txt
## how many reads mapped post duplicate removal
samtools view -c $4/$3.rmdup.sort.bam >> $4/$3_records.txt
## Coverage and bp mapped
samtools depth -a $4/$3.rmdup.sort.bam | awk '{sum+=$3;cnt++}END{print sum/cnt "\n" sum}' >> $4/$3_records.txt

## mapdamage
mapDamage -i $4/$3.rmdup.sort.bam --merge-reference-sequences --no-stats -r $5 -d $4/$3_mapdamage

rm $4/$3.sort.bam* $4/$3.rmdup.bam* $4/$3-trimmed.fastq $4/$3.rmdup.sort.bam*

