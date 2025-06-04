hg38_blacklist = rtracklayer::import.bed('~/RefDatabase/Bowtie2_Index/GRCh38_unified_blacklist.bed')
## 1. ATAC
atac_path = '~/ATAC_SEQ/IL17D_ATACSeq_KCs/Cleandata'
peak_fls = list.files(atac_path, pattern = 'ATAC_.*_peaks.narrowPeak$', full.names = T)
names(peak_fls) = gsub('_peaks.narrowPeak','',basename(peak_fls))
print(names(peak_fls))
ATAC_Con_peaks = load_peaks(peak_fls = peak_fls[1], blacklist = hg38_blacklist)
ATAC_IL17D_peaks = load_peaks(peak_fls = peak_fls[2], blacklist = hg38_blacklist)

ATAC_Unified_Peaks = c(ATAC_Con_peaks, ATAC_IL17D_peaks)
ATAC_Unified_Peaks = reduce(ATAC_Unified_Peaks[order(ATAC_Unified_Peaks)])
export(ATAC_Unified_Peaks, "ATAC_Unified_Peaks.bed")

## 2. H3K18la
chip_path = '~/chip_seq_kbd'
peak_fls = list.files(chip_path, pattern = '*_La_peaks.narrowPeak$', full.names = T)
names(peak_fls) = gsub('_peaks.narrowPeak','',basename(peak_fls))
print(names(peak_fls))
Con_La_Peak = load_peaks(peak_fls[1], blacklist = hg38_blacklist)
IL17D_La_Peak = load_peaks(peak_fls[2], blacklist = hg38_blacklist)
La_Unified_Peaks = c(Con_La_Peak, IL17D_La_Peak)
La_Unified_Peaks = reduce(La_Unified_Peaks[order(La_Unified_Peaks)])
export(La_Unified_Peaks, "La_Unified_Peaks.bed")

# H3K4me1
peak_fls = list.files(chip_path, pattern = '*_Me_peaks.narrowPeak$', full.names = T)
names(peak_fls) = gsub('_peaks.narrowPeak','',basename(peak_fls))
print(names(peak_fls))
Con_Me_Peak = load_peaks(peak_fls[1], blacklist = hg38_blacklist)
IL17D_Me_Peak = load_peaks(peak_fls[2], blacklist = hg38_blacklist)
Me_Unified_Peaks = c(Con_Me_Peak, IL17D_Me_Peak)
Me_Unified_Peaks = reduce(Me_Unified_Peaks[order(Me_Unified_Peaks)])
export(Me_Unified_Peaks, "Me_Unified_Peaks.bed")
