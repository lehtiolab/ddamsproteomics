#!/usr/bin/env Rscript

library(DEqMS)
library(reshape2)
library(matrixStats)

sampletable = read.table('sampletable', header=F, sep='\t', comment.char='', quote='', colClasses=c('character'))
# Add an X to sample groups that start with a number because R does that to header fields of input feats
colnames(sampletable) = c('ch', 'set', 'sample', 'group')
sampletable$group = sub('^([0-9])', 'X\\1', sampletable$group)
lookup = sampletable$group
names(lookup) = apply(cbind(sampletable[c('group', 'sample', 'set', 'ch')]), 1, paste, collapse='_')
names(lookup) = gsub('[^a-zA-Z0-9_]', '_', names(lookup))

feats = read.table('feats', header=T, sep="\t", comment.char="", quote="")
colnames(feats) = sapply(colnames(feats), function(x) sub('q.value', 'q-value', x))
featcol = colnames(feats)[1]
rownames(feats) = feats[,1]

# Remove possible internal standard
feats.nostd = feats[, !grepl('^X__POOL', colnames(feats))]
sampletable = sampletable[sampletable$group != 'X__POOL',]

# Get all features with more than 1 measurement in ALL sample groups, discard the rest
# First regex for all channels, exclude the Amount PSMs columns with lookbehind regex (need PCRE)
feats.quant = feats.nostd[, grepl('plex_.*(?<!PSMs)$', colnames(feats.nostd), perl=T)]
feats.quant$feat = rownames(feats.quant)
feat.quantcount = melt(feats.quant, id.vars='feat')
feat.quantcount$group = sapply(as.character(feat.quantcount$variable), function(x) lookup[[ sub('_[a-z]+[0-9]+plex', '', x) ]])
feat.quantcount = dcast(aggregate(value~feat+group, na.omit(feat.quantcount), length), feat~group)
tmpfeat = feat.quantcount$feat
feat.quantcount$feat = NULL
feat.quantcount[feat.quantcount < 2 ] = NA
feat.quantcount$feat = tmpfeat
filtered_feats_quantcount = na.omit(feat.quantcount)

# With those features, filter the quant feats and median(PSM counts) > 0 and not NA
names(feats)[1] = 'feat'
names(feats.nostd)[1] = 'feat'
quantpsmcols = grepl('_Amount.fully.quanted.PSMs', colnames(feats.nostd))
feats.nostd$median_psmcount = round(rowMedians(as.matrix(feats.nostd[quantpsmcols]), na.rm=T))
feats.filt = merge(filtered_feats_quantcount['feat', drop=F], feats.nostd, by='feat')
feats.filt = feats.filt[feats.filt$median_psmcount > 0,]
rownames(feats.filt) = feats.filt$feat
feats.psms = feats.filt$median_psmcount
feats.filt = feats.filt[, grepl('plex_.*(?<!PSMs)$', colnames(feats.filt), perl=T)]

# Take median PSMs, do lmFit
rownames(sampletable) = 1:nrow(sampletable)
design = model.matrix(~0+group, sampletable)
colnames(design) = gsub('group', '', colnames(design))
fit1 = lmFit(as.matrix(feats.filt), design=design)

# Get all contrasts, do eBayes
combinations = combn(as.character(unique(sampletable$group)), 2)
contrasts = c()
for (ix in 1:ncol(combinations)) {
  comb = combinations[,ix]
  contrasts = append(contrasts, paste(comb, collapse='-'))
}
cont <- makeContrasts(contrasts=contrasts, levels = design)
fit2 = eBayes(contrasts.fit(fit1, contrasts = cont))
fit2$count = feats.psms
fit3 = spectraCounteBayes(fit2)

# Report
outcols = c('logFC', 'count', 'sca.P.Value', 'sca.adj.pval')
outfeats = feats
for (col in 1:length(contrasts)) {
  cond_report = outputResult(fit3, coef_col=col)[outcols]
  names(cond_report) = sapply(names(cond_report), function(x) paste(contrasts[col], x, sep='_'))
  cond_report$feat = rownames(cond_report)
  outfeats = merge(outfeats, cond_report, by='feat', all.x=T)
}
names(outfeats)[1] = featcol
write.table(outfeats, 'deqms_output', sep='\t', row.names=F, quote=F)
