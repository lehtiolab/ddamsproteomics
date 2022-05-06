#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
nrsets = as.numeric(args[1])
usedptms = unlist(strsplit(args[2], ';'))
psmfn = args[3]
pepfn = args[4]

psms = read.table(psmfn, header=T, sep="\t", comment.char = "", quote = "")
peptides = read.table(pepfn, header=T, sep="\t", comment.char="", quote="")

width = 4

## Prepare a data frame with all the info for plots

# first split the sites into multi sites (e.g. GG and Phos, sep by _)
multiptms = strsplit(psms$Top.luciphor.PTM, '_', fixed=T)
repeater = sapply(multiptms, FUN=length)
sites = data.frame(sites=unlist(multiptms),
                   specid=rep(psms$SpecID, repeater),
                   bioset=rep(psms$Biological.set, repeater),
                   peptide=rep(psms$Peptide, repeater)
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

# For proteins remove the ones with multi hits (; for multiprot and / for multisite on a protein
# to create unique site info
if ('Master.protein.s.' %in% names(psms)) {
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

    # add proteins to summary
    num_uniprot_set = aggregate(protein~bioset + ptm_residue, protein_ptms[!duplicated(protein_ptms[c('protein', 'ptm', 'site', 'bioset')]),], length)
    num_uniprot_any = aggregate(protein ~ ptm_residue, protein_ptms[!duplicated(protein_ptms[c('protein', 'ptm', 'site')]),], length)
    summary = merge(summary, num_uniprot_set, all=T)
    anysetsum = merge(anysetsum, num_uniprot_any, all=T)
}


# Prepare simple counted PSM/pep/protein data
psm_count = aggregate(specid~bioset, sites[!duplicated(sites[c('specid', 'bioset')]),], length)
pep_count = aggregate(peptide~bioset, sites[!duplicated(sites[c('peptide', 'bioset')]),], length)
print(psm_count)
print(pep_count)
featcount_summ = merge(psm_count, pep_count, all=T)
colnames(featcount_summ) = c('bioset', 'ptmpsmcount', 'ptmpepcount')


# nr PSMs with PTMs
svg('ptmpsmfeats', width=width, height=nrsets+2)
  print(ggplot(psm_count) +
    geom_bar(aes(bioset, y=specid), stat='identity') +
    coord_flip() + ylab('# PSMs with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
  )
dev.off()

# nr peptides with PTMs
# FIXME Helptext - identical sequence with different PTMs counts as different (multiple) peptides in this plot
svg('ptmpepfeats', width=width, height=nrsets+2)
  print(ggplot(pep_count) +
    geom_bar(aes(bioset, y=peptide), stat='identity') +
    coord_flip() + ylab('# peptides with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
  )
dev.off()


# nr. proteins with PTMs
if ('Master.protein.s.' %in% names(psms)) {
    prot_count = aggregate(protein~bioset, protein_ptms[!duplicated(protein_ptms[c('protein', 'bioset')]),], length)

    svg('ptmprotfeats', width=width, height=nrsets+2)
      print(ggplot(prot_count) + 
        geom_bar(aes(bioset, y=protein), stat='identity') +
        coord_flip() + ylab('# proteins with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
      )
    dev.off()

    colnames(prot_count) = c('bioset', 'ptmprotcount')
    featcount_summ = merge(featcount_summ, prot_count, all=T)
}


svg('psmptmresidues', width=width, height=nrsets+2)
print(ggplot(summary) +
    geom_bar(aes(bioset, specid, fill=ptm_residue), stat='identity', position='dodge') +
    coord_flip() + ylab('# PTM residues found') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
)
dev.off()


svg('pepptmresidues', width=width, height=nrsets+2)
print(ggplot(summary) +
    geom_bar(aes(bioset, peptide, fill=ptm_residue), stat='identity', position='dodge') +
    coord_flip() + ylab('# PTM residues found') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
)
dev.off()


anysetsum$bioset = 'Combined sets'
summary = merge(summary, anysetsum, all=T)
write.table(summary, 'summary.txt', row.names=F, quote=F, sep='\t')
write.table(featcount_summ, 'featcount_summary.txt', row.names=F, quote=F, sep='\t')
