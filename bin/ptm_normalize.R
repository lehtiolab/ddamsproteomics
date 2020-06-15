#!/usr/bin/env Rscript

library(DEqMS)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)
psms = args[1]
#filtered_features = args[2]
setname = args[2]
# FIXME parse args, denomcols and samplenames are both optional
denomcols = args[3]
#samplenames = args?

psms = read.table(psms, header=T, sep="\t", comment.char="", quote="")
colnames(psms)[2] = "Feature"
#features = read.table(filtered_features, header=T, sep="\t", comment.char="", quote="")

# this is to ensure we get values for each l
psms[psms == 0] <- NA
psms = na.omit(psms)
lastcol = dim(psms)[2]
psms[,3:lastcol] = log2(psms[,3:lastcol])
psm.counts = as.data.frame(table(psms$Feature))
countout = data.frame(feat=psm.counts$Var1, matrix(psm.counts$Freq, nrow(psm.counts), 1))
write.table(countout, 'psmcounts', sep='\t', quote=F, row.names=F, col.names=F)
rownames(psm.counts) = psm.counts$Var1

# Normalize and filter features against the passed filtered_features table
medianfile = sprintf("%s_channelmedians", setname)
if (is.na(denomcols)) {
  print('Median sweeping')
  # Sweep from DEqMS, but without equalMedianNormalization
  group_col = 2
  dat.ratio = psms
  dat.ratio[,3:ncol(psms)] = dat.ratio[,3:ncol(psms)] - 
    matrixStats::rowMedians(as.matrix(dat.ratio[,3:ncol(psms)]),na.rm = TRUE)
  dat.summary = plyr::ddply(dat.ratio, colnames(psms)[group_col],
    function(x)matrixStats::colMedians(as.matrix(x[,3:ncol(psms)]),na.rm = TRUE))
  colnames(dat.summary)[2:ncol(dat.summary)] = colnames(psms)[3:ncol(psms)]
  peptides.nm = dat.summary[,-1]
  rownames(peptides.nm) = dat.summary[,1]
  #peptides.nm = peptides.nm[rownames(peptides.nm) %in% features[,1], ]
} else {
  print('Median summarizing with denominators')
  denomcols = as.numeric(strsplit(denomcols, ",")[[1]])
  peptides.nm = medianSummary(psms, group_col=2, ref_col=denomcols)
}
# FIXME colnames proteins.nm should get sample name prefixed, so CTRL_tmt10plex_126, TREAT_tmt10plex_127N etc
#colnames(proteins.nm)
featcol = data.frame(feats=rownames(peptides.nm))
write.table(cbind(featcol, peptides.nm), 'normalized_feats', sep='\t', quote=F, row.names=F)
