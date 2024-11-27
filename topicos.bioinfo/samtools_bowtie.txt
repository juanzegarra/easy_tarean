#!/bin/bash
for i in *.fasta
do 
bowtie -f -S -v 0 -a -p 20 -t ../index_cdhitest_control_viral_contigs $i > ${i}_a_v0.sam 2> ${i}_a_v0_bowtie.log
samtools view -S -h -F 4 ${i}_a_v0.sam > ${i}_a_v0_onlymapped.sam
samtools sort -O SAM -o ${i}_a_v0_onlymapped_sorted.sam ${i}_a_v0_onlymapped.sam
samtools view -Sb ${i}_a_v0.sam > ${i}_a_v0.bam
done
