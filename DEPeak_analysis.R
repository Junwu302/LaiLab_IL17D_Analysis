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

enrich_func = function(genelist){
  require(clusterProfiler, quietly = T)
  require(org.Hs.eg.db, quietly = T)
  kegg_res = setReadable(enrichKEGG(gene = genelist, organism = 'hsa', pAdjustMethod = 'none', pvalueCutoff = 1),
                         OrgDb = 'org.Hs.eg.db',keyType ='ENTREZID')@result
  kegg_res$qvalue = p.adjust(kegg_res$pvalue, method = 'BH')
  go_res = enrichGO(gene = genelist, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID', ont = 'ALL', pAdjustMethod = 'none',
                    pvalueCutoff = 1,qvalueCutoff=1, readable = TRUE)@result
  go_res$qvalue = p.adjust(go_res$pvalue, method = 'BH')
  return(list(GO=go_res, KEGG = kegg_res))
}

H3K4me1_DEPeaks = de_peaks(Counts = H3K4me1_Counts, Group = factor(rep(c('Con','IL17D'),each = 2), levels = c('Con','IL17D')))
H3K18la_DEPeaks = de_peaks(Counts = H3K18la_Counts, Group = factor(rep(c('Con','IL17D'),each = 2), levels = c('Con','IL17D')))
ATAC_DEPeaks = de_peaks(Counts = ATAC_Counts, Group = factor(rep(c('Con','IL17D'),each = 2), levels = c('Con','IL17D')))

# H3K18la
H3K18la_anno = merge(H3K18la_DEPeaks$DPeaks[,c('chr','start','end','logFC','PValue')], 
                     H3K18la_DEPeaks$peakAnno[!is.na(H3K18la_DEPeaks$peakAnno$SYMBOL),c('seqnames','start','end','annotation','SYMBOL','geneId')],
                     by.x = c('chr','start','end'), by.y = c('seqnames','start','end'))
genelist = unique(H3K18la_anno$geneId[grepl('Promoter', H3K18la_anno$annotation) & H3K18la_anno$PValue < 0.05])
H3K18la_enrichment = enrich_func(genelist = genelist)

# H3K4me1
H3K4me1_anno = merge(H3K4me1_DEPeaks$DPeaks[,c('chr','start','end','logFC','PValue')], 
                     H3K4me1_DEPeaks$peakAnno[!is.na(H3K4me1_DEPeaks$peakAnno$SYMBOL),c('seqnames','start','end','annotation','SYMBOL','geneId')],
                     by.x = c('chr','start','end'), by.y = c('seqnames','start','end'))
genelist = unique(H3K4me1_anno$geneId[grepl('Promoter', H3K4me1_anno$annotation) & H3K4me1_anno$PValue < 0.05])
H3K4me1_enrichment = enrich_func(genelist = genelist)

# ATAC-seq
ATAC_anno = merge(ATAC_DEPeaks$DPeaks[,c('chr','start','end','logFC','PValue')], 
                  ATAC_DEPeaks$peakAnno[!is.na(ATAC_DEPeaks$peakAnno$SYMBOL),c('seqnames','start','end','annotation','SYMBOL','geneId')],
                     by.x = c('chr','start','end'), by.y = c('seqnames','start','end'))
genelist = unique(ATAC_anno$geneId[grepl('Promoter', ATAC_anno$annotation) & ATAC_anno$PValue < 0.05])
ATAC_enrichment = enrich_func(genelist = genelist)
