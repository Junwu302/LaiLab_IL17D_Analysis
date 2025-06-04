computeMatrix reference-point --referencePoint center -S con-la-1.rmdup.bw con-la-3.rmdup.bw IL17D-la-1.rmdup.bw IL17D-la-3.rmdup.bw -R H3K18la_UpBed.bed --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o H3K18la_Up_matrix.mat.gz
computeMatrix reference-point --referencePoint center -S con-me-1.rmdup.bw con-me-2.rmdup.bw IL17D-me-1.rmdup.bw IL17D-me-2.rmdup.bw -R H3K4me1_DownBed.bed --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o H3K4me1_Down_matrix.mat.gz
computeMatrix reference-point --referencePoint center -S Con1.rmdup.bw Con2.rmdup.bw D1.rmdup.bw D2.rmdup.bw -R ATAC_DEBed.bed --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o ATAC_DE_matrix.mat.gz

plotHeatmap -m H3K18la_Up_matrix.mat.gz -out H3K18la_Up_peaks.png --colorList '#177cb0,#ffffff,#c93756'
plotHeatmap -m H3K4me1_Down_matrix.mat.gz -out H3K4me1_Down_peaks.png --colorList '#177cb0,#ffffff,#c93756'
plotHeatmap -m ATAC_DE_matrix.mat.gz -out ATAC_DEpeak.png --colorList '#177cb0,#ffffff,#c93756'
