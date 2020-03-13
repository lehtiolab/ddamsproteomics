#!/usr/bin/env Rscript

library(DEqMS)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)
psms = args[1]
filtered_features = args[2]
setname = args[3]
# FIXME parse args, denomcols and samplenames are both optional
denomcols = args[4]
#samplenames = args?

psms = read.table(psms, header=T, sep="\t", comment.char="", quote="")
colnames(psms)[2] = "Feature"
features = read.table(filtered_features, header=T, sep="\t", comment.char="", quote="")

# this is to ensure we get values for each l
psms[psms == 0] <- NA
psms = na.omit(psms)
lastcol = dim(psms)[2]
psms[,3:lastcol] = log2(psms[,3:lastcol])
psm.counts = as.data.frame(table(psms$Feature))
print(head(psm.counts)) # REMOCVE
countout = data.frame(feat=psm.counts$Var1, matrix(psm.counts$Freq, nrow(psm.counts), 1))
write.table(countout, 'psmcounts', sep='\t', quote=F, row.names=F, col.names=F)
rownames(psm.counts) = psm.counts$Var1

# Normalize and filter features against the passed filtered_features table
medianfile = sprintf("%s_channelmedians", setname)
if (is.na(denomcols)) {
  print('Median sweeping')
  proteins.nm = medianSweeping(psms, group_col=2)
  proteins.nm = proteins.nm[rownames(proteins.nm) %in% features[,1], ]
} else {
  print('Median summarizing with denominators')
  denomcols = as.numeric(strsplit(denomcols, ",")[[1]])
  proteins = medianSummary(psms, group_col=2, ref_col=denomcols)

  ##### This could be built into yafeng/deqms#1, deprecate when released, call sweep/equalmediannorm with medianfile instead
  sizefactor = matrixStats::colMedians(as.matrix(proteins),na.rm = TRUE)
  write.table(data.frame(channels=colnames(proteins), medians=sizefactor), medianfile, sep='\t', quote=F, row.names=F)
  #######

  ### Add setnames to medians table. This should be done also for sweep when deqms has release:
  chmedians = read.table(medianfile, sep='\t', header=T)
  write.table(cbind(data.frame(setname=setname), chmedians), medianfile, sep='\t', quote=F, row.names=F, col.names=F)
  ###

  proteins.nm = equalMedianNormalization(proteins)
}
# FIXME colnames proteins.nm should get sample name prefixed, so CTRL_tmt10plex_126, TREAT_tmt10plex_127N etc
colnames(proteins.nm)
featcol = data.frame(feats=rownames(proteins.nm))
write.table(cbind(featcol, proteins.nm), 'normalized_feats', sep='\t', quote=F, row.names=F)
