#!/usr/bin/bash
Ref=~/Bowtie2_Index/GRCh38
# 路径：/mnt/data/WuJun/Collaboration/LaiYunPing/WXX/ATAC_SEQ/IL17D_ATACSeq_KCs/Cleandata
for fq1 in `ls *.Cleandata.R1.fq.gz`
do
	name=${fq1%.Cleandata.R1.fq.gz}
	echo $name
	# 1. run mapping
	fq1=${name}/${name}_clean.R1.fq.gz
	fq2=${name}.Cleandata.R2.fq.gz
	bowtie2 -x $Ref -1 $fq1 -2 $fq2 -S ${name}.sam -p 50 2>&1 | tee -a ${name}.bowtie2.log
	samtools flagstat ${name}.sam > ${name}_stat.txt
	samtools view -h -@ 10 -f 2 -Sbu ${name}.sam | samtools sort -@ 10 -o ${name}.bam
	samtools index ${name}.bam
	samtools rmdup ${name}.bam ${name}.rmdup.bam
	samtools index ${name}.rmdup.bam
	rm ${name}.sam
	rm ${name}.bam
	rm ${name}.bam.bai
done
