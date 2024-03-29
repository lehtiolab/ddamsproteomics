#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
args = commandArgs(trailingOnly=TRUE)
nrsets = as.numeric(args[1])
has_fractions = args[2] == TRUE
oldmzmlfn = args[3]
oldmzmldef = ifelse(oldmzmlfn == 'nofile', FALSE, TRUE)
plateids = args[4:length(args)]

# Fraction needs to be a factor and not numeric, since it can contain strings (e.g. A2)
# That would be fine but when mixing oldmzmls from a rerun in, the join on Fraction op may fail
# when there is both the factor 02, and the numeric 2 in different tables
feats = read.table("psms", colClasses=c('Fraction'='factor'), header=T, sep="\t", comment.char = "", quote = "")
ycol = 'value'
if (has_fractions) {
  width = 4
  xcol ='plateID'
  feats$plateID = paste(feats$Biological.set, feats$Strip, sep='_')
  amount_ms2 = read.table("platescans", sep='\t', header=F)
  feats$Fraction = as.factor(feats$Fraction)
} else {
  width = 14
  xcol ='SpectraFile'
  amount_ms2 = read.table("filescans", sep="|", header=F)[,2:3]
}


##### PSM-scans
set_amount_psms = aggregate(SpecID~Biological.set, feats, length)
names(set_amount_psms) = c('Set', 'psmcount')
write.table(set_amount_psms, 'psmtable_summary.txt', row.names=F, quote=F, sep='\t')

amount_psms = aggregate(SpecID~get(xcol), feats, length)
mscol = 'MS2 scans'
psmcol = 'PSMs IDed'
names(amount_ms2) = c(xcol, mscol)
names(amount_psms) = c(xcol, psmcol)
amount_id=merge.data.frame(amount_ms2[c(xcol, mscol)], amount_psms[c(xcol, psmcol)], by.x=xcol, by.y=xcol, all=T)
amount_id = melt(amount_id, measure.vars=c(mscol, psmcol))

procents = dcast(amount_id, get(xcol)~variable, value.var=ycol)
colnames(procents)[1] = xcol
procents$p= procents$`PSMs IDed` / procents$`MS2 scans`
amount_id = merge(amount_id, procents[, c(xcol, 'p')])
amount_id$labeltext = ifelse(amount_id$variable == psmcol, paste(round(100 * amount_id$p, 2), '%'), '')

# FIXME move percent to be on tip of bar
svg('psm-scans', width=width, height=nrsets + 2) 
print(ggplot(amount_id) +
  geom_bar(aes(x=!!ensym(xcol), y=value, fill=variable, group=variable), stat='identity', position='dodge') + coord_flip() +
    ylab('# spectra') + theme_bw() + theme(axis.title.x=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
  geom_text(position=position_dodge(width=0.9), aes(y=value, x=!!ensym(xcol), group=variable, label=labeltext), hjust=0, colour="black", size=6))
dev.off()


# Missing isobaric values
if (length(grep('plex', names(feats)))) {
  channels = names(feats)[grepl('plex', names(feats))]
  psm_empty = melt(feats[c(xcol, channels)], id.vars=xcol)
  psm_empty = psm_empty[psm_empty$value==0,]
  if (nrow(psm_empty)) {
    psm_empty$value = 1
    psm_empty = aggregate(value~get(xcol)+ variable, psm_empty, sum)
    print('psm empty ready')
    names(psm_empty) = c(xcol, 'channels', 'nr_missing_values')
    psm_empty$channels = sub('.*plex_', '', psm_empty$channels)
    svg('missing-tmt', width=width, height=(nrsets + 2))
    print(ggplot(psm_empty) + 
      geom_bar(aes_string(x=xcol, y='nr_missing_values', fill='channels'), stat='identity', position="dodge") + ylab('# PSMs without quant') + coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()))
    dev.off()
  }
}

mcl = aggregate(as.formula(paste('SpecID~', xcol, '+ missed_cleavage')), feats, length)
mcl$missed_cleavage = as.factor(mcl$missed_cleavage)

mcl_am = subset(merge(mcl, amount_psms, by=xcol), missed_cleavage %in% c(0,1,2,3))
mcl_am$percent = mcl_am$SpecID / mcl_am$"PSMs IDed" * 100
mcl_am$textycoord = ifelse(mcl_am$missed_cleavage!=0, mcl_am$percent, 10)

svg('miscleav', width=width, height=(nrsets + 2))
mcplot = ggplot(mcl_am) +
    geom_bar(aes(x=!!ensym(xcol), y=percent, fill=missed_cleavage, group=missed_cleavage), position='dodge', stat='identity') +
    # 0.9 is the default dodge (90% of 1, 1 used bc all same value) but when not spec -> no dodge at all?
    geom_text(position=position_dodge(width=0.9), aes(x=!!ensym(xcol), y=textycoord+3, group=missed_cleavage, label=paste(SpecID, 'PSMs,', round(percent, 2), '%')), hjust=0, colour="black", size=4, inherit.aes=T) +
    ylim(c(0, 100)) + ylab('% of PSMs') +
    theme_bw() +
    theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
    coord_flip() 
print(mcplot)
dev.off()

# Now the per-fraction or per-file stats
if (has_fractions) {
  xcol = 'Fraction'
  fryield_form = 'SpecID ~ Fraction + plateID'
  # plateID not necessary because we take subfeats?
} else { 
  xcol = 'SpectraFile'
  fryield_form = 'SpecID ~ SpectraFile'
}

allfilefractions = read.table("mzmlfrs", colClasses=c(NA, NA, NA, 'factor'), header=F, sep='\t')
colnames(allfilefractions) = c('SpectraFile', 'setname', 'strip', 'Fraction')
if (oldmzmldef) {
  oldmzmls = read.table(oldmzmlfn, colClasses=c(NA, NA, NA, NA, 'factor'), header=F, sep='\t')[,c(1,3,4,5)]
  # TODO columns maybe added/disappea in future old/mzmldef file
  colnames(oldmzmls) = c('SpectraFile', 'setname', 'strip', 'Fraction')
  oldmzmls$setname = trimws(oldmzmls$setname)
  allfilefractions = rbind(allfilefractions, oldmzmls)
}
allfilefractions$plateID = paste(allfilefractions$setname, allfilefractions$strip, sep='_')

ptypes = list(
  retentiontime=c('Retention.time.min.', 'time(min)'),
  precerror=c('PrecursorError.ppm.', 'Precursor error (ppm)'),
  fryield=c('SpecID', '# PSMs'),
  msgfscore=c('MSGFScore', 'MSGF Score'),
  fwhm=c('FWHM', 'FWHM'),
  pif=c('Precursor.ion.fraction', 'Precursor/all in window'))

for (plateid in plateids) {
  if (has_fractions) {
    subfeats = subset(feats, plateID==plateid) 
    subfiles = subset(allfilefractions, plateID==plateid)
    h = 4
    w = 14
  } else { 
    subfeats = feats
    subfiles = allfilefractions
    h = nrow(unique(feats[xcol])) + 1
    w = 14
  }
  for (ptype in names(ptypes)) {
    fn = paste('PLATE', plateid, ptype, sep="___")
    if (ptypes[[ptype]][1] %in% colnames(subfeats)) {
      if (ptype == 'fryield' && nrow(subfeats)) {
        yieldcounts = aggregate(as.formula(fryield_form), subfeats, length)
        plotdata = merge(yieldcounts, subfiles[xcol], all.y=T)
        p = ggplot(plotdata) + geom_bar(aes_string(x=xcol, y=ptypes[[ptype]][1]), stat='identity')
      } else {
        plotdata = subfeats
        p = ggplot(plotdata, aes_string(x=xcol, y=ptypes[[ptype]][1])) + geom_violin(trim=F) 
      }
      if (ptype == 'precerror') {
        p = p + geom_hline(yintercept=0, size=2)
      }
      svg(fn, height=h, width=w)
      p = p + ylab(ptypes[[ptype]][2]) + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10))
      if(!has_fractions) p = p + xlab('Sample') + coord_flip()
      print(p)
      dev.off()
    }
  }
}
