de_peaks = function(Counts, Group){
  require(edgeR)
  require(ChIPseeker)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(org.Hs.eg.db)
  metaData = data.frame(Group, row.names = colnames(Counts)[-c(1:3)])
  DPeaks = DGEList(counts = Counts[,-c(1:3)], group = metaData$Group)
  CPM = cpm(DPeaks)
  keep = rowSums(CPM > 1) >= 2
  CPM = cbind(Counts[,1:3], CPM)
  DPeaks = DPeaks[keep, ,keep.lib.sizes = FALSE]
  DPeaks = calcNormFactors(DPeaks, method = 'TMM')
  design = model.matrix(~metaData$Group)
  DPeaks = estimateDisp(DPeaks, design, robust = TRUE)
  DPeaks = glmFit(DPeaks, design, robust = TRUE)
  DPeaks = glmLRT(DPeaks)$table
  DPeaks = cbind(Counts[keep,1:3], DPeaks)
  
  # annotation
  gr = GRanges(seqnames = Rle(DPeaks$chr), ranges = IRanges(DPeaks$start, DPeaks$end))
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  peakAnno = annotatePeak(gr, tssRegion=c(-2000, 2000), TxDb=txdb, annoDb="org.Hs.eg.db")
  peakAnno = data.frame(peakAnno@anno)
  return(list(DPeaks = DPeaks, CPM = CPM, peakAnno = peakAnno))
}

H3K4me1_DEPeaks = de_peaks(Counts = H3K4me1_Counts, Group = factor(rep(c('Con','IL17D'),each = 2), levels = c('Con','IL17D')))
H3K18la_DEPeaks = de_peaks(Counts = H3K18la_Counts, Group = factor(rep(c('Con','IL17D'),each = 2), levels = c('Con','IL17D')))
ATAC_DEPeaks = de_peaks(Counts = ATAC_Counts, Group = factor(rep(c('Con','IL17D'),each = 2), levels = c('Con','IL17D')))
