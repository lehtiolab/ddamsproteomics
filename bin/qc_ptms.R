#!/usr/bin/env Rscript

library(plotly)
library(htmltools)
library(tidyr)
library(glue)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
psmfn = args[1]
pepfn = args[2]
search_engine = args[3]

psms = read.table(psmfn, header=T, sep="\t", comment.char = "", quote = "")
peptides = read.table(pepfn, header=T, sep="\t", comment.char="", quote="")

if (search_engine == 'sage') {
  scancol = 'scannr'
  pepcol = 'peptide'

} else if (search_engine == 'msgf') {
  scancol = 'SpecID' # MSGF (sage has scannr)
  pepcol = 'Peptide'
}


## Prepare a data frame with all the info for plots

# first split the sites into multi sites (e.g. GG and Phos, sep by _)
# FIXME dyn col SpecID / scannr
multiptms = strsplit(psms$Top.luciphor.PTM, '_', fixed=T)
repeater = sapply(multiptms, FUN=length)

sites = data.frame(sites=unlist(multiptms),
                   specid=rep(psms[[ scancol ]], repeater),
                   bioset=rep(psms$Biological.set, repeater),
                   peptide=rep(psms[[ pepcol ]], repeater)
)
if ('Master.protein.s.' %in% names(psms)) {
    sites$protein = rep(psms$Master.protein.s., repeater)
} else {
    sites$protein = 'NA'
}

# Split the PTM/sites into multiple sites (e.g. Phos T21,Y45 sep by ,)
sites = cbind(sites$bioset, sites$specid, sites$peptide, sites$protein, stringr::str_split_fixed(sites$sites, ':', 2))
multisites = strsplit(sites[,6], ',', fixed=T)
repeater = sapply(multisites, FUN=length)
sites = data.frame(site=unlist(multisites),
                   bioset=rep(sites[,1], repeater),
                   specid=rep(sites[,2], repeater),
                   peptide=sub('_.*', '', rep(sites[,3], repeater)),
                   protein_ptms=rep(sites[,4], repeater),
                   ptm=rep(sites[,5], repeater)
                   )

# PSM/peptide summaries, first create filtering and aggregator columns
sites$ptm_residue = paste(sites$ptm, sub('[0-9]+', '', sites$site))
num_psm_set = aggregate(specid~bioset + ptm_residue, sites, length)
num_pep_set = aggregate(peptide~bioset + ptm_residue, sites[!duplicated(sites[c('ptm', 'site', 'bioset', 'peptide')]),], length)
num_psm_any = aggregate(specid ~ ptm_residue, sites, length)
num_pep_any = aggregate(peptide ~ ptm_residue, sites[!duplicated(sites[c('ptm', 'site', 'peptide')]),], length)
summary = merge(num_psm_set, num_pep_set, all=T)
anysetsum = merge(num_psm_any, num_pep_any, all=T)

if ('Master.protein.s.' %in% names(psms)) {
    # For proteins remove the ones with multi hits (; for multiprot and / for multisite on a protein
    # to create unique site info, then filter duplicated protein/site rows (this is from PSM info)
    uniprot = subset(sites, !grepl(';', protein_ptms, fixed=T))
    unisite = subset(uniprot, !grepl('__.+/', protein_ptms))

    # Now same for proteins as above for peptide/PSM, split annotated protein into multiptms/multisites
    protsites = cbind(unisite$bioset, unisite$protein_ptms, stringr::str_split_fixed(unisite$protein_ptms, '__', 2))
    protmultiptms = strsplit(protsites[,4], '_', fixed=T)
    repeater = sapply(protmultiptms, FUN=length)
    protmultiptms = data.frame(sites=unlist(protmultiptms),
                               protein=rep(protsites[,3], repeater),
                           bioset=rep(unisite$bioset, repeater)
                           )
    site_ptm = stringr::str_split_fixed(protmultiptms$sites, ':', 2)
    protmultiptms[,c('ptm', 'site')] = c(site_ptm[,1], site_ptm[,2])
    protmultisites = strsplit(protmultiptms$site, ',', fixed=T)
    repeater = sapply(protmultisites, FUN=length)
    protein_ptms = data.frame(site=unlist(protmultisites),
                              protein=rep(protmultiptms$protein, repeater),
                              bioset=rep(protmultiptms$bioset, repeater),
                              ptm=rep(protmultiptms$ptm, repeater)
                              )
    protein_ptms$ptm_residue = paste(protein_ptms$ptm, sub('[0-9]+', '', protein_ptms$site))
    protein_ptms = protein_ptms[!duplicated(protein_ptms[c('protein', 'ptm', 'site', 'bioset')]),]

    # add proteins to summary
    num_uniprot_set = aggregate(protein~bioset + ptm_residue, protein_ptms, length)
    # do another duplicate but without bioset for the all-set protein sites
    num_uniprot_any = aggregate(protein ~ ptm_residue, protein_ptms[!duplicated(protein_ptms[c('protein', 'ptm', 'site')]),], length)
    summary = merge(summary, num_uniprot_set, all=T)
    anysetsum = merge(anysetsum, num_uniprot_any, all=T)
}


# Prepare simple counted PSM/pep/protein data
psm_count = aggregate(specid~bioset, sites[!duplicated(sites[c('specid', 'bioset')]),], length)
pep_count = aggregate(peptide~bioset, sites[!duplicated(sites[c('peptide', 'bioset')]),], length)
featcount_summ = merge(psm_count, pep_count, all=T)
colnames(featcount_summ) = c('bioset', 'ptmpsmcount', 'ptmpepcount')

nrsets = nrow(unique(featcount_summ))

# Overlap of sites table
if (nrsets > 1) {
    pep_setcount = aggregate(bioset~peptide+site+ptm_residue, sites[!duplicated(sites[c('ptm', 'site', 'peptide', 'bioset')]),], length)
    pep_setcount$pepsite = paste(pep_setcount$peptide, pep_setcount$site)
    site_overlap = aggregate(pepsite~bioset+ptm_residue, pep_setcount, length)
    prot_setcount = aggregate(bioset~protein+site+ptm_residue, protein_ptms, length)
    prot_setcount$protsite = paste(prot_setcount$protein, prot_setcount$site)
    site_overlap = merge(site_overlap, aggregate(protsite~bioset+ptm_residue, prot_setcount, length), all=T)
    colnames(site_overlap) = c('nr_sets', 'ptm_residue', 'peplvl', 'protlvl')
    write.table(site_overlap, 'ptm__overlap.txt', row.names=F, quote=F, sep='\t')
}


# nr PSMs with PTMs
ggp = ggplot(psm_count) +
    geom_bar(aes(bioset, y=specid), stat='identity') +
    coord_flip() + ylab('# PSMs with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
p = ggplotly(ggp, width=400) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
# Work around since plotly does not honor above legend.title=element_blank call
p$x$layout$legend$title$text = ''
htmlwidgets::saveWidget(p, 'psm__ptms.html', selfcontained=F)


# nr peptides with PTMs
# FIXME Helptext - identical sequence with different PTMs counts as different (multiple) peptides in this plot
ggp = ggplot(pep_count) +
    geom_bar(aes(bioset, y=peptide), stat='identity') +
    coord_flip() + ylab('# peptides with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
p = ggplotly(ggp, width=400) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
# Work around since plotly does not honor above legend.title=element_blank call
p$x$layout$legend$title$text = ''
htmlwidgets::saveWidget(p, 'peptide__ptms.html', selfcontained=F)



# nr. proteins with PTMs
if ('Master.protein.s.' %in% names(psms)) {
    prot_count = aggregate(protein~bioset, protein_ptms[!duplicated(protein_ptms[c('protein', 'bioset')]),], length)

  ggp = ggplot(prot_count) + 
        geom_bar(aes(bioset, y=protein), stat='identity') +
        coord_flip() + ylab('# proteins with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
  p = ggplotly(ggp, width=400) %>%
          layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
  # Work around since plotly does not honor above legend.title=element_blank call
  p$x$layout$legend$title$text = ''
  htmlwidgets::saveWidget(p, 'protein__ptms.html', selfcontained=F)

    colnames(prot_count) = c('bioset', 'ptmprotcount')
    featcount_summ = merge(featcount_summ, prot_count, all=T)
}


ggp = ggplot(summary) +
    geom_bar(aes(bioset, specid, fill=ptm_residue), stat='identity', position='dodge') +
    coord_flip() + ylab('# PTM residues found') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
p = ggplotly(ggp, width=400) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
# Work around since plotly does not honor above legend.title=element_blank call
p$x$layout$legend$title$text = ''
htmlwidgets::saveWidget(p, 'psm__residues.html', selfcontained=F)



ggp = ggplot(summary) +
    geom_bar(aes(bioset, peptide, fill=ptm_residue), stat='identity', position='dodge') +
    coord_flip() + ylab('# PTM residues found') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
p = ggplotly(ggp, width=400) %>%
        layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
# Work around since plotly does not honor above legend.title=element_blank call
p$x$layout$legend$title$text = ''
htmlwidgets::saveWidget(p, 'peptide__residues.html', selfcontained=F)


anysetsum$bioset = 'Combined sets'
summary = merge(summary, anysetsum, all=T)
write.table(summary, 'ptm__table.txt', row.names=F, quote=F, sep='\t')
write.table(featcount_summ, 'ptm__featcount_table.txt', row.names=F, quote=F, sep='\t')
