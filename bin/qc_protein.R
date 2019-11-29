#!/usr/bin/env Rscript

library(ggplot2)
library(forcats)
library(reshape2)
library(ggrepel)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--sets', dest='sets', type='character', nargs='+')
parser$add_argument('--feattype', type='character')
parser$add_argument('--peptable', type='character')
parser$add_argument('--sampletable', type='character', default=FALSE)
parser$add_argument('--normtable', type='character')
opt = parser$parse_args()

#args = commandArgs(trailingOnly=TRUE)
nrsets = length(opt$sets)
setnames = opt$sets
feattype = opt$feattype
peptable = opt$peptable
sampletable = opt$sampletable
feats = read.table("feats", header=T, sep="\t", comment.char = "", quote = "")

featcol = list(peptides='Peptide.sequence', proteins='Protein.ID', genes='Gene.ID', assoc='Gene.Name')[[feattype]]

# nrpsms
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[grep('quanted_psm_count', colnames(feats))]
  nrpsms = melt(feats, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_[a-z0-9]*plex_[0-9NC]*_quanted_psm_count', '', nrpsms$variable)
  nrpsms$Set = sub('_quanted_psm_count', '', nrpsms$Set)
  summary_psms = aggregate(value~Set, nrpsms, median)
  colnames(summary_psms) = c('Set', paste('no_psm_', feattype, sep=''))
  nrpsms = aggregate(value~get(featcol)+Set, nrpsms, max)
  nrpsms = transform(nrpsms, setrank=ave(value, Set, FUN = function(x) rank(x, ties.method = "random")))
  png('nrpsms')
  print(ggplot(nrpsms, aes(y=value, x=setrank)) +
    geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
    theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
    scale_x_reverse())
    dev.off()
}

qcols = colnames(feats)[grep('_q.value', colnames(feats))]
overlap = na.exclude(feats[qcols])
overlap = dim(overlap[apply(overlap, 1, function(x) any(x<0.01)),])[1]
if (feattype == 'peptides') {
  am_prots = melt(feats, id.vars=c(featcol, "Protein.s."), measure.vars=qcols)
  am_prots$nrprots = lengths(regmatches(am_prots$Protein.s., gregexpr(';', am_prots$Protein.s.))) + 1
  # aggregate feats and remove col with ; where?
} else {
  am_prots = melt(feats, id.vars=featcol, measure.vars=qcols)
  pepcols = colnames(feats)[grep('Amount.Unique', colnames(feats))]
  pepprots = melt(feats, id.vars=featcol, measure.vars=pepcols)
  pepprots$variable = sub('_Amount.Unique.peptides', '', pepprots$variable)
  pepmed = aggregate(value~variable, pepprots, median)
  colnames(pepmed) = c('Set', paste('no_pep_', feattype, sep=''))
}
am_prots = am_prots[!is.na(am_prots$value),]
am_prots = am_prots[am_prots$value < 0.01,]
am_prots$Set = sub('_q.value', '', am_prots$variable)
png('featyield', height=(nrsets + 2) * 72)
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
  write.table(am_prots[,c(1,3)], 'summary.txt', row.names=F, quote=F, sep='\t')
  am_prots = melt(am_prots)
  colnames(am_prots)[3] = 'accession'
  print(ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank(), legend.text=element_text(size=20), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=20)) +
    geom_bar(aes(fct_rev(Set), y=accession, fill=variable), stat='identity', position='dodge') +
    geom_text(data=subset(am_prots, variable=='All'), aes(fct_rev(Set), accession/2, label=accession), colour="white", size=7, nudge_x=-0.25) + 
    geom_text(data=subset(am_prots, variable=='Non-shared (unique)'), aes(fct_rev(Set), accession/2, label=accession), colour="white", size=7, nudge_x=+0.25) + 
	ggtitle(paste('Overlap for all sets: ', overlap, '\nTotal uniques: ', totalunique)))
} else {
  am_prots = aggregate(get(featcol) ~ Set, am_prots, length)
  colnames(am_prots)[2] = 'accession'
  missing = setdiff(setnames, am_prots$Set)
  missingvals = vector(mode='integer', length=length(missing))
  missing_df = data.frame(Set=missing, accession=missingvals)
  am_prots = rbind(am_prots, missing_df)
  if (length(grep('plex', names(feats)))) {
    # if isobaric, then show summary table of feats 1%FDR AND quant
    tmtcols = colnames(feats)[setdiff(grep('plex', colnames(feats)), grep('quanted', colnames(feats)))]
    not_fullna = feats[rowSums(is.na(feats[,tmtcols])) != length(tmtcols),]
    sum_prots = melt(not_fullna, id.vars=featcol, measure.vars=qcols)
    sum_prots = sum_prots[!is.na(sum_prots$value),]
    sum_prots = sum_prots[sum_prots$value < 0.01,]
    sum_prots$Set = sub('_q.value', '', sum_prots$variable)
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
  write.table(summary, 'summary.txt', row.names=F, quote=F, sep='\t')
  print(ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank(), plot.title=element_text(size=20)) +
    geom_bar(aes(fct_rev(Set), y=accession), stat='identity') +
    geom_text(aes(fct_rev(Set), accession/2, label=accession), colour="white", size=10) + ggtitle(paste('Overlap for all sets: ', overlap, '\nTotal identified: ', nrow(feats))))
}
dev.off()

# coverage
if (feattype == 'proteins') {
  covmed = median(feats$Coverage)
  png('coverage')
  plot = ggplot(feats) + geom_histogram(aes(Coverage), bins=50) + theme_bw()
  plot = plot + geom_label(x=max(ggplot_build(plot)$data[[1]]$x)*0.75, y=max(ggplot_build(plot)$data[[1]]$count)*0.75, label=sprintf('Median: %.3f', covmed))
  print(plot)
  dev.off()
}

#isobaric
# first get a fullsamplename to set lookup, if we have a sampletable
use_sampletable = FALSE
if (sampletable) {
  use_sampletable = TRUE
  sampletable = read.table('sampletable', header=F, sep='\t', comment.char='', quote='', colClasses=c('character'))
  colnames(sampletable) = c('ch', 'set', 'sample', 'group')
  lookup = sampletable$set
  rownames(sampletable) = apply(sampletable[c('group', 'sample', 'set', 'ch')], 1, paste, collapse='_')
  names(lookup) = apply(sampletable[c('group', 'sample', 'set', 'ch')], 1, paste, collapse='_')
  names(lookup) = gsub('[^a-zA-Z0-9_-]', '_', names(lookup))
}

if (length(grep('plex', names(feats)))) {
  tmtcols = colnames(feats)[setdiff(grep('plex', colnames(feats)), grep('quanted', colnames(feats)))]
  overlap = na.exclude(feats[c(tmtcols, qcols)])
  overlap = dim(overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),])[1]
  tmt = melt(feats, id.vars=featcol, measure.vars = tmtcols)
  if (use_sampletable) {
    tmt$Set = apply(tmt, 1, function(x) { key = sub('_[a-z0-9]*plex', '', x[["variable"]]); print(key); return (lookup[[key]]) })
  } else { 
    tmt$Set = sub('_[a-z0-9]*plex.*', '', tmt$variable)
  }
  tmt$variable = sub('.*_[a-z].*[0-9]*plex_', '', tmt$variable)
  outplot = ggplot(na.exclude(tmt)) + geom_boxplot(aes(fct_rev(Set), value, fill=fct_rev(variable)), position=position_dodge(width=1)) +
    coord_flip() + ylab('Fold change') + xlab('Channels') + theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), plot.title=element_text(size=30) ) + 
    theme(legend.text=element_text(size=20), legend.position="top", legend.title=element_blank()) +
    ggtitle(paste('Overlap with values in \nall ', length(tmtcols), 'channels: ', overlap))
  if (min(na.exclude(tmt$value)) >= 0) { outplot = outplot + scale_y_log10() }
  png('isobaric', height=(3 * nrsets + 1) * 72)
  print(outplot)
  dev.off()
}

#nrpsmsoverlapping
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[grep('quanted_psm_count', colnames(feats))]
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[grep('plex', colnames(feats))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),]
  if (nrow(overlap) > 0) {
    nrpsms = melt(overlap, id.vars=featcol, measure.vars = nrpsmscols)
    nrpsms$Set = sub('_[a-z0-9]*plex_[0-9NC]*_quanted_psm_count', '', nrpsms$variable)
    nrpsms$Set = sub('_quanted_psm_count', '', nrpsms$Set)
    if (feattype == 'peptides') {
      nrpsms = aggregate(value~Peptide.sequence+Set, nrpsms, max)
    } else {
      nrpsms = aggregate(value~get(featcol)+Set, nrpsms, max)
    }
    nrpsms = transform(nrpsms, setrank=ave(value, Set, FUN = function(x) rank(x, ties.method = "random")))
    png('nrpsmsoverlapping')
    print(ggplot(nrpsms, aes(y=value, x=setrank)) +
      geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
      theme_bw() + 
      theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
      scale_x_reverse())
      dev.off()
  }
}


# percentage_onepsm
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[grep('quanted_psm_count', colnames(feats))]
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[grep('plex', colnames(feats))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),]
  if (nrow(overlap) > 0) {
    nrpsms = melt(overlap, id.vars=featcol, measure.vars = nrpsmscols)
    nrpsms$Set = sub('_quanted_psm_count', '', nrpsms$variable) # this is enough for DEqMS pipeline
    # but need also following line for normal pipe bc nr-psms values are reported per channel
    nrpsms$Set = sub('_[a-z0-9]*plex.*', '', nrpsms$Set)
    feats_in_set = aggregate(value~Set, data=nrpsms, length) 
    feats_in_set$percent_single = aggregate(value~Set, data=nrpsms, function(x) length(grep('[^01]', x)))$value / feats_in_set$value * 100
    png('percentage_onepsm')
    print(ggplot(feats_in_set, aes(Set, percent_single)) +
      geom_col(aes(fill=Set)) + theme_bw() + ylab('% of identifications') +
      theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) )
    dev.off()
  }
}

# ranked step plot MS1 peptide per protein
if (feattype != 'peptides') {
  peps = read.table(peptable, header=T, sep='\t', comment.char='', quote='')
  ms1qcols = grep('MS1.area', colnames(peps))
  nrpep_set = melt(peps, id.vars=c("Protein.s."), measure.vars=ms1qcols, na.rm=T)
  if (dim(nrpep_set)[1] != 0) {
    nrpep_set$Set = sub('_MS1.area.*', '', nrpep_set$variable)
    nrpep_set = aggregate(variable~Protein.s.+Set, nrpep_set, length) 
    nrpep_set = transform(nrpep_set, setrank=ave(variable, Set, FUN = function(x) rank(x, ties.method = "random")))
    png('ms1nrpeps')
    print(ggplot(nrpep_set, aes(y=variable, x=setrank)) +
      geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# peptides with MS1') +
      theme_bw() + 
      theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
      scale_x_reverse())
    dev.off()
  }
}

# precursorarea
if (length(grep('area', names(feats)))) {
    parea = melt(feats, id.vars=featcol, measure.vars = colnames(feats)[grep('area', colnames(feats))])
    parea$Set = sub('_MS1.*', '', parea$variable)
    png('precursorarea', height=(nrsets + 1) * 72)
    print(ggplot(parea) + 
      geom_boxplot(aes(fct_rev(Set), value)) + scale_y_log10() + coord_flip() + ylab("Intensity") + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank()))
    dev.off()
}


# DEqMS volcano plots
deqpval_cols = grep('_sca.P.Value$', names(feats))
deqFC_cols = grep('_logFC$', names(feats))
names(feats)[1] = 'feat'
if (length(deqpval_cols)) {
  s_table = unique(sampletable[sampletable$group != 'X__POOL', 'group'])
  cartprod = expand.grid(s_table, s_table)
  cartprod = cartprod[cartprod$Var1 != cartprod$Var2,]
  for (comparison in paste(cartprod$Var1, cartprod$Var2, sep='.')) {
    logfcname = sprintf('%s_logFC', comparison) 
    if (length(grep(logfcname, names(feats)))) {
      compnice = sub('[.]', ' vs. ', comparison)
      logpname = sprintf('%s_log.sca.pval', comparison)
      feats[, logpname] = -log10(feats[, sprintf('%s_sca.P.Value', comparison)])
      png(sprintf('deqms_volcano_%s', comparison))
      plot = ggplot(feats, aes(x=get(sprintf('%s_logFC', comparison)), y=get(logpname), label=feat)) +
        geom_point(size=0.5 )+ theme_bw(base_size = 16) + # change theme
        theme(axis.title=element_text(size=25), axis.text=element_text(size=20)) +
        xlab(sprintf("log2 FC(%s)", compnice)) + # x-axis label
        ylab('-log10 P-value') + # y-axis label
        geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
        geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
        geom_vline(xintercept = 0, colour = "black") # Add 0 lines
      if (feattype != 'peptides') {
	topfeats = feats[order(feats[logpname], decreasing=TRUE)[1:10], ]
        plot = plot + geom_text_repel(data=topfeats)
      }
      print(plot)
      dev.off()
    }
  }
}


# PCA
if (length(deqpval_cols)) {
  pca_ana <- prcomp(t(na.omit(feats[,tmtcols])), scale. = TRUE)
  score.df <- as.data.frame(pca_ana$x)
  rownames(score.df) = sub('_[a-z0-9]*plex', '', rownames(score.df))
  score.df$type = sampletable[rownames(score.df), "group"]

  #Scree plot
  contributions <- data.frame(contrib=round(summary(pca_ana)$importance[2,] * 100, 2)[1:20])
  contributions$pc = rownames(contributions)
  png('scree')
  print(ggplot(data=contributions, aes(x=reorder(pc, -contrib), y=contrib)) +
    geom_bar(stat='identity') +
    theme_bw() + theme(axis.title=element_text(size=25), axis.text=element_text(size=15)) +
    ylab("Contribution (%)"))
  dev.off()
  png('pca')
  print(ggplot(data=score.df, aes(x =PC1, y =PC2, label=rownames(score.df), colour=type)) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_point(size=4) +
    theme_bw() + theme(axis.title=element_text(size=25), axis.text=element_text(size=20),
		       legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
    xlab(sprintf("PC1 (%s%%)", contributions$contrib[1])) + ylab(sprintf("PC2 (%s%%)", contributions$contrib[2]))
    )
  dev.off()
}
