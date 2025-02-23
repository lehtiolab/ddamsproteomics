#!/usr/bin/env Rscript
library(plotly)
library(htmltools)
library(tidyr)
library(glue)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
has_fractions = args[1] == TRUE
search_engine = args[2]

# Fraction needs to be a factor and not numeric, since it can contain strings (e.g. A2)
# That would be fine but when mixing oldmzmls from a rerun in, the join on Fraction op may fail
# when there is both the factor 02, and the numeric 2 in different tables
feats = read.table("psms_clean", colClasses=c('Fraction'='factor'), header=T, sep="\t", comment.char = "", quote = "")

if (search_engine == 'sage') {
  scancol = 'scannr'
  miscleavcol = 'missed_cleavages'
  rtcol = 'rt'
  ioncol = 'Ion.injection.time.ms.'
  precerrcol = 'precursor_ppm'
  scorecol = 'sage_discriminant_score'
  filenamecol = 'filename'

} else if (search_engine == 'msgf') {
  scancol = 'SpecID' # MSGF (sage has scannr)
  miscleavcol = 'missed_cleavage'
  rtcol = 'Retention.time.min.' # 'rt'
  ioncol = 'Ion.injection.time.ms.' # ?
  precerrcol = 'PrecursorError.ppm.' # precursor_ppm
  scorecol = 'MSGFScore' # 'sage_discriminant_score'
  filenamecol = 'SpectraFile' # 'filename'
}

tmtcolors = c(
  "#b15928",
  "#e31a1c",
  "#fb9a99",
  "#ff7f00",
  "#fdbf6f",
  "gold4",
  "gold1",
  "#33a02c",
  "#b2df8a",
  "#01665e",
  "#80cdc1",
  "#1f78b4",
  "#a6cee3",
  "#6a3d9a",
  "#cab2d6",
  "#c51b7d",
  "#f1b6da",
  "magenta4"
  )

boxplot_stats = function(data, col) {
  summary_stats = data %>%
  summarise(lower = quantile({{col}}, 0.25, na.rm=T),
            middle = median({{col}}, na.rm=T),
            upper = quantile({{col}}, 0.75, na.rm=T),
            minval = min({{col}}, na.rm=T),
            maxval = max({{col}}, na.rm=T),
    )
  summary_stats$iqr = summary_stats$upper - summary_stats$lower
  summary_stats$whisk_min = pmax(summary_stats$minval, summary_stats$middle - summary_stats$iqr * 1.5)
  summary_stats$whisk_max = pmin(summary_stats$maxval, summary_stats$middle + summary_stats$iqr * 1.5)
  return(summary_stats)
}


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
amount_id$plate_fn = amount_id[[xcol]]
write.table(amount_id, 'psms_ids.txt', row.names=F, quote=F, sep='\t')

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
    allchannels = unique(psm_empty$channels)
    ggp = ggplot(psm_empty, aes(x=.data[[ xcol ]], y=nr_missing_values, fill=channels)) + 
      geom_bar(stat='identity', position="dodge") + 
      scale_fill_manual(values=tmtcolors[1:length(allchannels)]) +
      ylab('# PSMs without quant') + coord_flip() + theme_bw() +
      theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text=element_text(size=10), axis.text.y=element_text(angle=90), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()
    )
    p = ggplotly(ggp, width=400, height=vert_height) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
  p$x$layout$legend$title$text = ''
  htmlwidgets::saveWidget(p, 'iso_missing_vals.html', selfcontained=F)
  }
}

# Missed cleavages
mcl = aggregate(get(scancol)~get(xcol)+get(miscleavcol), feats, length)
colnames(mcl) = c(xcol, 'missed_cleavage', 'nrpsms')
mcl_am = subset(merge(mcl, amount_psms, by=xcol), missed_cleavage %in% c(0,1,2))
mcl_am$percent = mcl_am$nrpsms / mcl_am$"PSMs IDed" * 100
mc_text_y = max(mcl_am$percent) * 2/6
mcl_am$missed_cleavage = as.factor(mcl_am$missed_cleavage)

mcplot = ggplot(mcl_am) +
    geom_bar(aes(x=.data[[ xcol ]], y=percent, fill=missed_cleavage, group=missed_cleavage), position='dodge', stat='identity') +
    # 0.9 is the default dodge (90% of 1, 1 used bc all same value) but when not spec -> no dodge at all?
    geom_text(position=position_dodge(width=0.9), aes(x=.data[[xcol]], y=mc_text_y, group=missed_cleavage, label=glue('{nrpsms} PSMs')), colour="black", size=4, inherit.aes=T) +
    ylim(c(0, 100)) + ylab('% of PSMs') +
    theme_bw() +
    theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text=element_text(size=10), axis.text.y=element_text(angle=90), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
    coord_flip() 
p = ggplotly(mcplot, width=400, height=vert_height) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
p$x$layout$legend$title$text = ''
htmlwidgets::saveWidget(p, 'missed_cleavages.html', selfcontained=F)
mcl_am$plate_fn = mcl_am[[xcol]]
write.table(mcl_am, 'miscleav.txt', row.names=F, quote=F, sep='\t')


# Now the per-fraction or per-file stats
allfilefractions = read.table('mzmls', colClasses='character', header=F, sep='\t', fill=T)
colnames(allfilefractions) = c('mzmlfile', 'setname', 'plate', 'fraction')
allfilefractions$setname = trimws(allfilefractions$setname)
allfilefractions$Fraction = allfilefractions$fraction
allfilefractions$plateID = paste(allfilefractions$setname, allfilefractions$plate, sep='_')
allfilefractions[[ filenamecol ]]  = basename(allfilefractions$mzmlfile)

# List of plot names and their axis label:
ptypes = list(
  retentiontime=c(rtcol, 'time(min)'),
  ioninjtime=c(ioncol, 'time(ms)'),
  precerror=c(precerrcol, 'Precursor error (ppm)'),
  fryield=c(scancol, '# PSMs'),
  score=c(scorecol, scorecol), # Diff score name for diff SE
  pif=c('Precursor.ion.fraction', 'Precursor/all in window')
)

if (has_fractions) {
  frcol = 'Fraction'
  plate_ids = amount_ms2[[xcol]]
} else { 
  frcol = filenamecol
  plate_ids = c('noplates')
}
for (plateid in plate_ids) {
  if (has_fractions) {
    subfeats = subset(feats, get(xcol)==plateid) 
    subfiles = subset(allfilefractions, get(xcol)==plateid)
  } else { 
    subfeats = feats
    subfiles = allfilefractions
  }
  for (ptype in names(ptypes)) {
    if (nrow(subfeats) && ptypes[[ptype]][1] %in% colnames(subfeats)) {
      if (ptype == 'fryield' && nrow(subfeats)) {
        fryield_form = glue('{scancol} ~ {frcol}')
        yieldcounts = aggregate(as.formula(fryield_form), subfeats, length)
        plotdata = merge(yieldcounts, subfiles[frcol], all.y=T)
        p = ggplot(plotdata, aes(x=.data[[ frcol ]], y=.data[[ ptypes[[ptype]][1] ]])) + geom_bar(stat='identity')
      } else {
        summary_stats <- subfeats[c(frcol, ptypes[[ptype]][1])] %>%
	    group_by(.data[[frcol]]) %>%
	    boxplot_stats(.data[[ptypes[[ptype]][1]]])
        # Plot using geom_crossbar, geom_errorbar
        p = ggplot(summary_stats, aes(x=.data[[ frcol ]], y=middle)) +
          geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), position=position_dodge(width=1) ) +
          geom_crossbar(aes(ymin = lower, ymax=upper), fill='white', linewidth=0.15, position=position_dodge(width=1))
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
