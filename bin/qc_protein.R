#!/usr/bin/env Rscript

library(plotly)
library(htmltools)
library(tidyr)
library(glue)
library(stringr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--sets', dest='sets', type='character', nargs='+')
parser$add_argument('--feattype', type='character')
parser$add_argument('--peptable', type='character')
parser$add_argument('--sampletable', type='character', default=FALSE)
parser$add_argument('--normtable', type='character', default=FALSE)
parser$add_argument('--conflvl', type='double')
opt = parser$parse_args()

nrsets = length(opt$sets)
vert_height = 100 * nrsets + 100
vert_height_iso = 300 * nrsets + 150
setnames = opt$sets
setnames = gsub('[^a-zA-Z0-9_]', '.', setnames)
setnames = sub('^([0-9])', 'X\\1', setnames)
feattype = opt$feattype
peptable = opt$peptable
conflvl = opt$conflvl
sampletable = opt$sampletable
feats = read.table("feats", header=T, sep="\t", comment.char = "", quote = "")

featcol = list(peptides='Peptide.sequence', proteins='Protein.ID', genes='Gene.Name', ensg='Gene.ID')[[feattype]]

# Summary tables of overlap
qcols = colnames(feats)[grep('_q.value', colnames(feats))]
qvals_long = pivot_longer(feats, any_of(qcols), names_to='Set', values_to='qval')
qvals_long$Set = sub('_q.value', '', qvals_long$Set)
feats_nrsets = aggregate(qval~get(featcol), qvals_long, length)
colnames(feats_nrsets) = c('feat', 'nrsets')
feat_set_overlap = aggregate(feat~nrsets, feats_nrsets, length)
write.table(feat_set_overlap, glue('{feattype}__overlap'), row.names=F, quote=F, sep='\t')


is_isobaric = length(grep('plex', names(feats)))
if (is_isobaric) {
  tmtcols = colnames(feats)[setdiff(grep('plex', colnames(feats)), grep('[qQ]uanted', colnames(feats)))]
  nrpsmscols = colnames(feats)[grep('_Fully.quanted.PSM.count', colnames(feats))]
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

# nrpsms used for isobaric quant
if (is_isobaric) {
  nrpsms = pivot_longer(feats, any_of(nrpsmscols), names_to='psmset', values_to='psmcount') 
  nrpsms$Set = sub('_Fully.quanted.PSM.count', '', nrpsms$psmset)
  summary_psms = aggregate(psmcount~Set, nrpsms, median)
  colnames(summary_psms) = c('Set', paste('no_psm_', feattype, sep=''))
  nrpsms = aggregate(psmcount~get(featcol)+Set, nrpsms, max)
  nrpsms = transform(nrpsms, setrank=ave(psmcount, Set, FUN = function(x) rank(x, ties.method = "random")))
  nrpsms$Set = sub('^X([^a-zA-Z])', '\\1', nrpsms$Set)
  nrpsms_plot = nrpsms[c('Set', 'psmcount', 'setrank')] %>%
    group_by(Set, psmcount) %>%
    summarise(maxrank = max(setrank))
  ggp = ggplot(nrpsms_plot, aes(y=psmcount, x=maxrank)) +
    geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
    theme_bw() + 
    theme(axis.title=element_text(size=15), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
    scale_x_reverse()
  p = ggplotly(ggp, width=400) %>%
          layout(legend = list(title='', orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
  htmlwidgets::saveWidget(p, 'iso_nrpsms.html', selfcontained=F)
}

# Amount of features
overlap = na.exclude(feats[qcols])
overlap = length(overlap[apply(overlap, 1, function(x) any(x<conflvl)),])
if (feattype == 'peptides') {
  qvals_long$nrprots = lengths(regmatches(qvals_long$Protein.s., gregexpr(';', qvals_long$Protein.s.))) + 1
} else {
  pepcols = colnames(feats)[grep('Unique.peptide.count', colnames(feats))]
  pepprots = pivot_longer(feats, any_of(pepcols), names_to='Set', values_to='uniquecount')
  pepprots$Set = sub('_Unique.peptide.count', '', pepprots$Set)
  pepmed = aggregate(uniquecount~Set, pepprots, median)
  colnames(pepmed) = c('Set', glue('no_pep_{feattype}'))
}
am_prots = qvals_long[!is.na(qvals_long$qval),]
am_prots = am_prots[am_prots$qval < conflvl,]

if (feattype == 'peptides') {
  totalunique = length(unique(subset(am_prots, nrprots == 1)$Peptide.sequence))
  unipepprotnr = aggregate(nrprots ~ Set, subset(am_prots, nrprots == 1), length)
  am_prots = aggregate(Peptide.sequence ~ Set, am_prots, length)
  am_prots = merge(am_prots, unipepprotnr)
  names(am_prots) = c('Set', 'All', 'Non-shared (unique)')
  missing = setdiff(setnames, am_prots$Set)
  missingvals = vector(mode='integer', length=length(missing))
  missing_df = data.frame(Set=missing, All=missingvals, nonshared=missingvals)
  names(missing_df) = c('Set', 'All', 'Non-shared (unique)')
  am_prots = rbind(am_prots, missing_df)
  am_prots$Set = sub('^X([^a-zA-Z])', '\\1', am_prots$Set)
  write.table(am_prots[,c(1,3)], glue('{feattype}__summary.txt'), row.names=F, quote=F, sep='\t')
  am_prots = pivot_longer(am_prots, !Set, names_to='peptype', values_to='pepcount')

  ggp = ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15)) +
    geom_bar(aes(Set, y=pepcount, fill=peptype), stat='identity', position='dodge') +
    geom_text(data=subset(am_prots, peptype=='All'), aes(Set, pepcount/2, label=pepcount), colour="white", size=7, nudge_x=-0.25) + 
    geom_text(data=subset(am_prots, peptype=='Non-shared (unique)'), aes(Set, pepcount/2, label=pepcount), colour="white", size=7, nudge_x=+0.25)
    writeLines(c(glue('Overlap for all sets: {overlap}'), glue('Total uniques: {totalunique}')), 'nrfeats__text.html')

} else {
  am_prots = aggregate(get(featcol) ~ Set, am_prots, length)
  colnames(am_prots)[2] = 'accession'
  missing = setdiff(setnames, am_prots$Set)
  missingvals = vector(mode='integer', length=length(missing))
  missing_df = data.frame(Set=missing, accession=missingvals)
  am_prots = rbind(am_prots, missing_df)
  if (is_isobaric) {
    # if isobaric, then show summary table of feats 1%FDR AND quant
    not_fullna = feats[rowSums(is.na(feats[,tmtcols])) != length(tmtcols),]
    sum_prots = pivot_longer(not_fullna, any_of(qcols), names_to='Set', values_to='qval')
    sum_prots = sum_prots[!is.na(sum_prots$qval),]
    sum_prots = sum_prots[sum_prots$qval < conflvl,]
    sum_prots$Set = sub('_q.value', '', sum_prots$Set)
    sum_prots = aggregate(get(featcol) ~ Set, sum_prots, length)
    summary = merge(pepmed, sum_prots, by='Set', all.y=T)
    colnames(summary)[ncol(summary)] = paste('nr_', feattype, '_q', sep='')
    summary = merge(summary, summary_psms, by='Set', all.y=T)
  } else {
    # just nr of proteins 1% FDR, no quant
    summary = merge(pepmed, am_prots, by='Set', all.y=T)
    colnames(summary)[ncol(summary)] = paste('nr_', feattype, sep='')
  }
  summary[is.na(summary)] = 0
  summary$Set = sub('^X([^a-zA-Z])', '\\1', summary$Set)
  am_prots$Set = sub('^X([^a-zA-Z])', '\\1', am_prots$Set)
  write.table(summary, glue('{feattype}__summary.txt'), row.names=F, quote=F, sep='\t')
  ggp = ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), plot.title=element_text(size=15)) +
    geom_bar(aes(Set, y=accession), stat='identity') +
    geom_text(aes(Set, accession/2, label=accession), colour="white", size=10) 
  writeLines(c(glue('Overlap for all sets: {overlap}'), glue('Total identified: {nrow(feats)}')), 'nrfeats__text.html')
}

p = ggplotly(ggp, width=400, height=vert_height) %>%
          layout(title='', legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
htmlwidgets::saveWidget(p, 'nrfeats.html', selfcontained=F)


# coverage
if (feattype == 'proteins') {
  covmed = median(feats$Coverage)
  max_x = max(feats$Coverage)
  ggp = ggplot(feats) + geom_histogram(aes(Coverage), bins=50) + theme_bw()
  # Two steps so ggplot_build can return something
  ggp = ggp + geom_text(x=max(feats$Coverage) * 0.75, y=max(ggplot_build(ggp)$data[[1]]$count)*0.75, label=sprintf('Median: %.3f', covmed))
  p = ggplotly(ggp, width=400)
  htmlwidgets::saveWidget(p, 'coverage.html', selfcontained=F)
}

#isobaric
# first get a fullsamplename to set lookup, if we have a sampletable
use_sampletable = FALSE
if (is.character(sampletable)) {
  use_sampletable = TRUE
  sampletable = read.table(sampletable, header=F, sep='\t', comment.char='', quote='', colClasses=c('character'))
  colnames(sampletable) = c('ch', 'set', 'sample', 'group')
  setlookup = sampletable$set
  rownames(sampletable) = apply(sampletable[c('group', 'sample', 'set', 'ch')], 1, paste, collapse='_')
  rownames(sampletable) = gsub('[^a-zA-Z0-9_]', '_', rownames(sampletable))
  rownames(sampletable) = sub('^([0-9])', 'X\\1', rownames(sampletable))
  names(setlookup) = rownames(sampletable)
}


boxplot_stats = function(data, col) {
  summary_stats = data %>%
  summarise(lower = quantile({{col}}, 0.25),
            middle = median({{col}}),
            upper = quantile({{col}}, 0.75),
            minval = min({{col}}),
            maxval = max({{col}}),
    )
  summary_stats$iqr = summary_stats$upper - summary_stats$lower
  summary_stats$whisk_min = pmax(summary_stats$minval, summary_stats$middle - summary_stats$iqr * 1.5)
  summary_stats$whisk_max = pmin(summary_stats$maxval, summary_stats$middle + summary_stats$iqr * 1.5)
  return(summary_stats)
}


if (is_isobaric) {
  # First produce the boxplots
  overlap = na.exclude(feats[c(tmtcols, qcols)])
  overlap = dim(overlap[apply(overlap[qcols], 1, function(x) any(x<conflvl)),])[1]
  tmt = pivot_longer(feats, any_of(tmtcols), names_to='channelset', values_to='value')
  if (use_sampletable) {
    tmt$Set = apply(tmt, 1, function(x) { key = sub('_[a-z0-9]*plex', '', x[['channelset']]); return (setlookup[[key]]) })
  } else { 
    tmt$Set = sub('_[a-z0-9]*plex.*', '', tmt$channelset)
  }
  tmt$channel = sub('.*_[a-z].*[0-9]*plex_', '', tmt$channelset)
  # Do not use geom_boxplot, we dont want ALL of the data frame in the HTML, this only passes
  # the relevant values and thus keeps it smaller
  summary_stats <- na.exclude(tmt[c('Set', 'value', 'channel')]) %>%
	    group_by(Set, channel) %>%
	    boxplot_stats(value)
  # Plot using geom_crossbar, geom_errorbar
  allchannels = unique(summary_stats$channel)
  ggp = ggplot(summary_stats, aes(x=Set, y=middle, fill=channel, width=0.35)) +
    scale_fill_manual(values=tmtcolors[1:length(allchannels)]) +
    scale_colour_manual(values=tmtcolors[1:length(allchannels)]) +
    coord_flip() + ylab('Fold change') +
    xlab('Channels') + theme_bw() + 
    theme(axis.title=element_text(size=10), axis.text=element_text(size=10)) + 
    theme(legend.text=element_text(2), legend.title=element_blank()) +
    geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), width=0.07, position=position_dodge(width=1) ) +
    geom_crossbar(aes(ymin = lower, ymax=upper, color=channel), width=0.1, position=position_dodge(width=1)) +
    guides(linetype='none', size='none', fill='none', shape='none')
  
  if (min(na.exclude(tmt$value)) >= 0) { ggp = ggp + scale_y_log10() }
  p = ggplotly(ggp, width=400, height=vert_height_iso) %>%
    layout(boxmode='group',
      legend = list(title='', orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom')
    )
  writeLines(c(glue('Overlap with values in all {length(tmtcols)} channels: {overlap}')), 'isobaric__text.html')
  htmlwidgets::saveWidget(p, 'isobaric.html', selfcontained=F)

  # Also produce normalization factor plot
  if (opt$normtable != FALSE) {
    norms = read.table(opt$normtable, header=F, sep='\t', comment.char='', quote='')
    colnames(norms) = c('Set', 'variable', 'value')
    norms$variable = sub('[a-z].*[0-9]*plex_', '', norms$variable)
    ggp = ggplot(norms, aes(Set, value, group=variable)) + 
      geom_col(aes(fill=variable), position=position_dodge(width=1)) +
      scale_fill_manual(values=tmtcolors[1:length(allchannels)]) +
      geom_text(aes(y=min(value), label=round(value, 4)), position=position_dodge(width=1), colour="gray20", size=3, hjust=0) +
      coord_flip() + ylab('Normalizing factor') + xlab('Channels') + theme_bw() + 
      theme(axis.title=element_text(size=15), axis.text=element_text(size=10)) + 
      theme(legend.text=element_text(size=10), legend.position="top", legend.title=element_blank())
    p = ggplotly(ggp, width=400, height=vert_height_iso) %>%
            layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    htmlwidgets::saveWidget(p, 'normfactors.html', selfcontained=F)
  }

  #nrpsms in features that are overlapping, i.e. complete in every set
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<conflvl)),]
  if (nrow(overlap) > 0) {
    nrpsms = pivot_longer(overlap, any_of(nrpsmscols), names_to='Set', values_to='psmcount')
    nrpsms$Set = sub('_Fully.quanted.PSM.count', '', nrpsms$Set)
    if (feattype == 'peptides') {
      agg_nrpsms = aggregate(psmcount~Peptide.sequence+Set, nrpsms, max)
    } else {
      agg_nrpsms = aggregate(psmcount~get(featcol)+Set, nrpsms, max)
    }
    agg_nrpsms = transform(agg_nrpsms, setrank=ave(psmcount, Set, FUN = function(x) rank(x, ties.method = "random")))
    agg_nrpsms$Set = sub('^X([^a-zA-Z])', '\\1', agg_nrpsms$Set)
    agg_nrpsms_plot = agg_nrpsms[c('Set', 'psmcount', 'setrank')] %>%
      group_by(Set, psmcount) %>%
      summarise(maxrank = max(setrank))
    ggp = ggplot(agg_nrpsms_plot, aes(y=psmcount, x=maxrank)) +
      geom_step(aes(color=Set), linewidth=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
      theme_bw() + 
      theme(axis.title=element_text(size=15), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
      scale_x_reverse()
    p = ggplotly(ggp, width=400) %>%
            layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    htmlwidgets::saveWidget(p, 'nrpsmsoverlapping.html', selfcontained=F)
  }

  # percentage_onepsm
  if (nrow(overlap) > 0) {
    feats_in_set = aggregate(psmcount~Set, data=nrpsms, length) 
    feats_in_set$percent_single = aggregate(psmcount~Set, data=nrpsms, function(x) length(grep('[^01]', x)))$psmcount / feats_in_set$psmcount * 100
    feats_in_set$Set = sub('^X([^a-zA-Z])', '\\1', feats_in_set$Set)
    ggp = ggplot(feats_in_set, aes(Set, percent_single)) +
      geom_col(aes(fill=Set)) + theme_bw() + ylab('% of identifications') +
      theme(axis.title=element_text(size=15), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank())
    p = ggplotly(ggp, width=400) %>%
            layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    htmlwidgets::saveWidget(p, 'percentage_onepsm.html', selfcontained=F)
  }
}

# ranked step plot MS1 peptide per protein
if (feattype != 'peptides') {
  peps = read.table(peptable, header=T, sep='\t', comment.char='', quote='')
  ms1qcols = grep('MS1.area', colnames(peps))
  if (length(ms1qcols)) {
    nrpep_set = pivot_longer(peps, any_of(ms1qcols), names_to='Set', values_to='ms1')
  #, na.rm=T)
    if (nrow(na.exclude(nrpep_set['ms1']))) {
      nrpep_set$Set = sub('_MS1.area.*', '', nrpep_set$Set)
      nrpep_set$Set = sub('^X([^a-zA-Z])', '\\1', nrpep_set$Set)
      nrpep_set = aggregate(ms1~Protein.s.+Set, nrpep_set, length) 
      nrpep_set = transform(nrpep_set, setrank=ave(ms1, Set, FUN = function(x) rank(x, ties.method = "random")))
      nrpep_set_plot = nrpep_set[c('Set', 'ms1', 'setrank')] %>%
        group_by(Set, ms1) %>%
        summarise(psmcount = max(setrank))
      ggp = ggplot(nrpep_set_plot, aes(y=ms1, x=psmcount)) +
        geom_step(aes(color=Set), linewidth=2) + scale_y_log10() + xlab('Rank') + ylab('# peptides with MS1') +
        theme_bw() + 
        theme(axis.title=element_text(size=15), axis.text=element_text(size=10), legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
        scale_x_reverse()
      p = ggplotly(ggp, width=400) %>%
              layout(legend = list(title='', orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
      htmlwidgets::saveWidget(p, 'ms1nrpeps.html', selfcontained=F)
    }
  }
}

# precursorarea
precursorcols = colnames(feats)[grep('area', colnames(feats))]
if (length(precursorcols)) {
    if (nrow(feats[precursorcols] %>% filter(if_any(everything(), ~ !is.na(.))))) {
      parea = pivot_longer(feats, any_of(precursorcols), names_to='Set', values_to='ms1')
      parea$Set = sub('_MS1.*', '', parea$Set)
      parea$Set = sub('^X([^a-zA-Z])', '\\1', parea$Set)
      parea = parea[complete.cases(parea$ms1), ]

  summary_stats <- parea[c('Set', 'ms1')] %>%
	    group_by(Set) %>%
	    boxplot_stats(ms1)

  # Plot using geom_crossbar, geom_errorbar
  ggp = ggplot(summary_stats, aes(x=Set, y=middle)) +
    ylab('Intensity') +
    theme_bw() + 
    scale_y_log10() + coord_flip() +
    theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank()) +
    geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), position=position_dodge(width=1) ) +
    geom_crossbar(aes(ymin = lower, ymax=upper), fill='white', linewidth=0.15, position=position_dodge(width=1))
    } else {
      ggp = ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw() +
            geom_text(aes(5, 50, label=sprintf('No MS1 found in %s', feattype)))
    }
    p = ggplotly(ggp, width=400, height=vert_height) #%>%
    htmlwidgets::saveWidget(p, 'precursorarea.html', selfcontained=F)
}

# DEqMS volcano plots
deqpval_cols = grep('_sca.P.Value$', names(feats))
deqFC_cols = grep('_logFC$', names(feats))
names(feats)[1] = 'feat'
if (length(deqpval_cols)) {
  s_table = unique(sampletable[sampletable$group != 'NO__GROUP', 'group'])
  s_table = sub('[^a-zA-Z0-9_]', '_', s_table)
  s_table = sub('^([0-9])', 'X\\1', s_table)
  cartprod = expand.grid(s_table, s_table)
  cartprod = cartprod[cartprod$Var1 != cartprod$Var2,]
  for (comparison in paste(cartprod$Var1, cartprod$Var2, sep='.')) {
    re_logfcname = sprintf('^%s_logFC', comparison) 
    if (length(grep(re_logfcname, names(feats)))) {
      compnice = sub('[.]', ' vs. ', comparison)
      logpname = sprintf('%s_log.sca.pval', comparison)
      feats[, logpname] = -log10(feats[, sprintf('%s_sca.P.Value', comparison)])
      ggp = ggplot(feats, aes(x=get(sprintf('%s_logFC', comparison)), y=get(logpname), label=feat)) +
        geom_point(size=0.5 )+ theme_bw(base_size = 16) + # change theme
        theme(axis.title=element_text(size=15), axis.text=element_text(size=10)) +
        xlab(sprintf("log2 FC(%s)", compnice)) + # x-axis label
        ylab('-log10 P-value') + # y-axis label
        geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
        geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
        geom_vline(xintercept = 0, colour = "black") # Add 0 lines
      if (feattype != 'peptides' && nrow(na.omit(feats[logpname]))) {
	    topfeats = feats[order(feats[,logpname], decreasing=TRUE)[1:10], ]
        ggp = ggp + geom_text(data=topfeats)
      }
      png(glue('deqms_volcano_{comparison}.png'))
      print(ggp)
      dev.off()
    }
  }
}

run_pca = function(scoredf, colortype) {
    scoredf$type = sampletable[rownames(scoredf), colortype]
    #Scree plot
    contributions <- data.frame(contrib=round(summary(pca_ana)$importance[2,] * 100, 2)[1:20])
    contributions$pc = sub('PC', '', rownames(contributions))
    ggp = ggplot(data=contributions, aes(x=reorder(pc, -contrib), y=contrib)) +
      geom_bar(stat='identity') +
      theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10)) +
      ylab("Contribution (%)") + xlab('PC (ranked by contribution)')
    p = ggplotly(ggp, width=400) %>%
            layout(legend = list(title='', orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    htmlwidgets::saveWidget(p, glue('scree_{colortype}.html'), selfcontained=F)

    # PCA plot
    ggp = ggplot(data=scoredf, aes(x=PC1, y=PC2, label=rownames(scoredf), colour=type)) +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_point(size=4) +
      theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10),
  		       legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
      xlab(sprintf("PC1 (%s%%)", contributions$contrib[1])) + ylab(sprintf("PC2 (%s%%)", contributions$contrib[2]))
    p = ggplotly(ggp, width=400) %>%
            layout(legend = list(title='', orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    htmlwidgets::saveWidget(p, glue('pca_{colortype}.html'), selfcontained=F)
}


# PCA
if (is_isobaric && use_sampletable) {
  topca = na.omit(feats[,tmtcols])
  if (nrow(topca)) {
    pca_ana <- prcomp(t(topca), scale. = TRUE)
    score.df <- as.data.frame(pca_ana$x)
    rownames(score.df) = sub('_[a-z0-9]*plex', '', rownames(score.df))
    if (length(sampletable[sampletable$group != 'NO__GROUP', 'group'])) {
        colortypes = c('group', 'set')
    } else {
        colortypes = c('set')
    }
    for (colortype in colortypes) {
      run_pca(score.df, colortype)
    }
  }
}
