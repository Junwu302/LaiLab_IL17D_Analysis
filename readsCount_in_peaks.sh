#!usr/bin/bash
chip_path=~ChIP_SEQ/IL17D_ChIP_SEQ_KCs/Cleandata
Con_la_bam1=${chip_path}/con-la-1/con-la-1.rmdup.bam
Con_la_bam2=${chip_path}/con-la-3/con-la-1.rmdup.bam
IL17D_la_bam1=${chip_path}/IL17D-la-1/IL17D-la-1.rmdup.bam
IL17D_la_bam2=${chip_path}/IL17D-la-3/IL17D-la-1.rmdup.bam

Con_me_bam1=${chip_path}/con-me-1/con-me-1.rmdup.bam
Con_me_bam2=${chip_path}/con-me-2/con-me-2.rmdup.bam
IL17D_me_bam1=${chip_path}/IL17D-me-1/IL17D-me-1.rmdup.bam
IL17D_me_bam2=${chip_path}/IL17D-me-2/IL17D-me-2.rmdup.bam

atac_path=~ATAC_SEQ/IL17D_ATACSeq_KCs/Cleandata
Con_atac_bam1=${atac_path}/Con1.rmdup.bam
Con_atac_bam2=${atac_path}/Con2.rmdup.bam
IL17D_atac_bam1=${atac_path}/D1.rmdup.bam
IL17D_atac_bam2=${atac_path}/D2.rmdup.bam

multiBamSummary BED-file -p 30 --ignoreDuplicates --smartLabels --BED ATAC_Unified_Peaks.bed --bamfiles ${Con_atac_bam1} ${Con_atac_bam2} ${IL17D_atac_bam1} ${IL17D_atac_bam2} --extendReads --centerReads --outRawCounts ATACSeq_Count.txt
multiBamSummary BED-file -p 30 --ignoreDuplicates --smartLabels --BED La_Unified_Peaks.bed --bamfiles ${Con_la_bam1} ${Con_la_bam2} ${IL17D_la_bam1} ${IL17D_la_bam2} --extendReads --centerReads --outRawCounts H3K18la_Count.txt
multiBamSummary BED-file -p 30 --ignoreDuplicates --smartLabels --BED Me_Unified_Peaks.bed --bamfiles ${Con_me_bam1} ${Con_me_bam2} ${IL17D_me_bam1} ${IL17D_me_bam2} --extendReads --centerReads --outRawCounts H3K4me1_Count.txt
