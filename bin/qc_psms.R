#!/usr/bin/env Rscript
library(plotly)
library(htmltools)
library(tidyr)
library(glue)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
has_fractions = args[1] == TRUE
newmzmlfn = ifelse(args[2] == 'FALSE', '', args[2])
oldmzmlfn = ifelse(args[3] == 'FALSE', '', args[3])

# Fraction needs to be a factor and not numeric, since it can contain strings (e.g. A2)
# That would be fine but when mixing oldmzmls from a rerun in, the join on Fraction op may fail
# when there is both the factor 02, and the numeric 2 in different tables
feats = read.table("psms_clean", colClasses=c('Fraction'='factor'), header=T, sep="\t", comment.char = "", quote = "")

scancol = 'SpecID' # MSGF (sage has scannr)
miscleavcol = 'missed_cleavage'
rtcol = 'Retention.time.min.' # 'rt'
ioncol = 'Ion.injection.time.ms.' # ?
precerrcol = 'PrecursorError.ppm.' # precursor_ppm
scorecol = 'MSGFScore' # 'sage_discriminant_score'
filenamecol = 'SpectraFile' # 'filename'


if (has_fractions) {
  xcol ='plateID'
  feats$plateID = paste(feats$Biological.set, feats$Strip, sep='_')
  amount_ms2 = read.table("platescans", sep='\t', header=F)
  feats$Fraction = as.factor(feats$Fraction)
} else {
  xcol = filenamecol
  amount_ms2 = read.table("filescans", sep="\t", header=F)[,2:3]
}

nr_verts = length(unique(amount_ms2[[1]]))
vert_height = 200 * nr_verts + 200

##### PSM-scans
# Count PSMs per set
set_amount_psms = aggregate(get(scancol)~Biological.set, feats, length)
names(set_amount_psms) = c('Set', 'psmcount')
set_amount_psms$Set = gsub('[^A-Za-z0-9_]', '.', set_amount_psms$Set)
write.table(set_amount_psms, 'psmtable__summary.txt', row.names=F, quote=F, sep='\t')

# Count PSMs per plate
amount_psms = aggregate(get(scancol)~get(xcol), feats, length)
mscol = 'MS2 scans'
psmcol = 'PSMs IDed'
names(amount_ms2) = c(xcol, mscol)
names(amount_psms) = c(xcol, psmcol)
amount_id=merge.data.frame(amount_ms2[c(xcol, mscol)], amount_psms[c(xcol, psmcol)], by.x=xcol, by.y=xcol, all=T)

procents = data.frame(x=amount_id[xcol],  p=amount_id[psmcol] / amount_id[mscol])
colnames(procents) = c(xcol, 'p')

amount_id = pivot_longer(amount_id, !any_of(xcol), values_to='count')
amount_id = merge(amount_id, procents[, c(xcol, 'p')], all.x=TRUE)
amount_id$labeltext = ifelse(amount_id$name == psmcol, paste(round(100 * amount_id$p, 2), '%'), '')

colnames(procents)[1] = xcol
ggp = ggplot(amount_id) +
  geom_bar(aes(x=.data[[xcol]], y=count, fill=name, group=name), stat='identity', position='dodge') + 
  geom_text(position=position_dodge(width=0.9), aes(y=max(amount_ms2[mscol]) * 2/6, x=.data[[ xcol ]], group=name, label=labeltext), colour="black", size=6) +
  coord_flip() + 
  ylab('# spectra') + 
  theme_bw() + 
  theme(axis.title.x=element_text(size=15), axis.text=element_text(size=10, angle=90), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank())
p = ggplotly(ggp, width=400, height=vert_height) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
# Work around since plotly does not honor above legend.title=element_blank call
p$x$layout$legend$title$text = ''

htmlwidgets::saveWidget(p, 'amount_psms.html', selfcontained=F)

# Missing isobaric values
if (length(grep('plex', names(feats)))) {
  channels = names(feats)[grepl('plex', names(feats))]
  psm_empty = pivot_longer(feats[c(xcol, channels)], !any_of(xcol), values_to='count')
  psm_empty = psm_empty[psm_empty$count == 0,]
  if (nrow(psm_empty)) {
    psm_empty$value = 1
    psm_empty = aggregate(value~get(xcol)+ name, psm_empty, sum)
    names(psm_empty) = c(xcol, 'channels', 'nr_missing_values')
    psm_empty$channels = sub('.*plex_', '', psm_empty$channels)
    ggp = ggplot(psm_empty, aes(x=.data[[ xcol ]], y=nr_missing_values, fill=channels)) + 
      geom_bar(stat='identity', position="dodge") + ylab('# PSMs without quant') + coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text=element_text(size=10), axis.text.y=element_text(angle=90), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank())
    p = ggplotly(ggp, width=400, height=vert_height) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
  p$x$layout$legend$title$text = ''
  htmlwidgets::saveWidget(p, 'iso_missing_vals.html', selfcontained=F)
  }
}

# Missed cleavages
mcl = aggregate(get(scancol)~get(xcol)+get(miscleavcol), feats, length)
colnames(mcl) = c(xcol, 'missed_cleavage', 'nrscan')
mcl_am = subset(merge(mcl, amount_psms, by=xcol), missed_cleavage %in% c(0,1,2))
mcl_am$percent = mcl_am$nrscan / mcl_am$"PSMs IDed" * 100
mc_text_y = max(mcl_am$percent) * 2/6
mcl_am$missed_cleavage = as.factor(mcl_am$missed_cleavage)

mcplot = ggplot(mcl_am) +
    geom_bar(aes(x=.data[[ xcol ]], y=percent, fill=missed_cleavage, group=missed_cleavage), position='dodge', stat='identity') +
    # 0.9 is the default dodge (90% of 1, 1 used bc all same value) but when not spec -> no dodge at all?
    geom_text(position=position_dodge(width=0.9), aes(x=.data[[xcol]], y=mc_text_y, group=missed_cleavage, label=glue('{nrscan} PSMs')), colour="black", size=4, inherit.aes=T) +
    ylim(c(0, 100)) + ylab('% of PSMs') +
    theme_bw() +
    theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
    coord_flip() 
p = ggplotly(mcplot, width=400, height=vert_height) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
p$x$layout$legend$title$text = ''
htmlwidgets::saveWidget(p, 'missed_cleavages.html', selfcontained=F)


# Now the per-fraction or per-file stats
if (has_fractions) {
  xcol = 'Fraction'
} else { 
  xcol = filenamecol
}
fryield_form = glue('{scancol} ~ {xcol}')

if (str_length(newmzmlfn)) {
  newmzmls = read.table(newmzmlfn, colClasses='character', header=T, sep='\t', fill=T)
}
if (str_length(oldmzmlfn)) {
  allfilefractions = read.table(oldmzmlfn, colClasses='character', header=T, sep='\t')
  if (str_length(newmzmlfn)) {
    allfilefractions = rbind(allfilefractions, newmzmls)
  }
}  else {
  allfilefractions = newmzmls
}
allfilefractions$setname = trimws(allfilefractions$setname)
allfilefractions$Fraction = allfilefractions$fraction
allfilefractions$plateID = paste(allfilefractions$setname, allfilefractions$plate, sep='_')
allfilefractions$filename = basename(allfilefractions$mzmlfile)

# List of plot names and their axis label:
ptypes = list(
  retentiontime=c(rtcol, 'time(min)'),
  ioninjtime=c(ioncol, 'time(ms)'),
  precerror=c(precerrcol, 'Precursor error (ppm)'),
  fryield=c(scancol, '# PSMs'),
  score=c(scorecol, scorecol), # Diff score name for diff SE
  pif=c('Precursor.ion.fraction', 'Precursor/all in window')
)

for (plateid in amount_ms2$plateID) {
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
    if (nrow(subfeats) && ptypes[[ptype]][1] %in% colnames(subfeats)) {
      if (ptype == 'fryield' && nrow(subfeats)) {
        yieldcounts = aggregate(as.formula(fryield_form), subfeats, length)
        plotdata = merge(yieldcounts, subfiles[xcol], all.y=T)
        p = ggplot(plotdata, aes(x=.data[[ xcol ]], y=.data[[ ptypes[[ptype]][1] ]])) + geom_bar(stat='identity')
      } else {
        plotdata = subfeats
        p = ggplot(plotdata, aes(x=.data[[ xcol ]], y=.data[[ ptypes[[ptype]][1] ]])) + geom_violin(trim=F) 
      }
      if (ptype == 'precerror') {
        p = p + geom_hline(yintercept=0, size=2)
      }
      if(!has_fractions) p = p + xlab('Sample') + coord_flip()
      ggp = p + ylab(ptypes[[ptype]][2]) + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10))

      if(has_fractions) {
        p = ggplotly(ggp, width=1300)
      } else {
        p = ggplotly(ggp, width=1300, height=vert_height)
      }
      htmlwidgets::saveWidget(p, glue('PLATE___{plateid}___{ptype}.html'), selfcontained=F)
    }
  }
}
