#!/usr/bin/bash
macs2=~/Softwares/MACS2-2.2.9.1/bin/macs2
chromSize=~/RefDatabase/Bowtie2_Index/GRCh38.chrom.sizes
for bam in `ls *.rmdup.bam`
do
	name=${bam%.rmdup.bam}
	echo $name
	$macs2 callpeak -f BAMPE -B --SPMR -t $bam -n $name -q 0.05 --shift -100 --extsize 200 --nomodel -g hs --outdir ./ 2> ${name}.macs2.log
	bedGraphToBigWig ${name}}_treat_pileup.bdg ${chromSize} ${name}}.bw
done
