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

# nr PSMs with PTMs
set_amount_psms = aggregate(SpecID~Biological.set, psms, length)
svg('ptmpsmfeats', width=width, height=nrsets+2)
  print(ggplot(set_amount_psms) +
    geom_bar(aes(Biological.set, y=SpecID), stat='identity') +
    coord_flip() + ylab('# PSMs with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
  )
dev.off()

# nr peptides with PTMs
qcols = colnames(peptides)[grep('_q.value', colnames(peptides), fixed=T)]
set_amount_pep = melt(peptides, id.vars='Peptide.sequence', measure.vars=qcols)
set_amount_pep = set_amount_pep[!is.na(set_amount_pep$value),]
set_amount_pep$variable = sub('_q.value', '', set_amount_pep$variable)
pepnums = aggregate(value~variable, set_amount_pep, length)
svg('ptmpepfeats', width=width, height=nrsets+2)
  print(ggplot(pepnums) +
    geom_bar(aes(variable, y=value), stat='identity') +
    coord_flip() + ylab('# peptides with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
  )
dev.off()

# nr. proteins with PTMs
if ('Master.protein.s.' %in% names(psms)) {
    prots = psms[c('Master.protein.s.', 'Biological.set')]
    prots = subset(prots, !grepl(';', Master.protein.s., fixed=T))
    prots = prots[!duplicated(prots),]
    set_amount_prots = aggregate(Master.protein.s.~Biological.set, prots, length)
    svg('ptmprotfeats', width=width, height=nrsets+2)
      print(ggplot(set_amount_prots) +
        geom_bar(aes(Biological.set, y=Master.protein.s.), stat='identity') +
        coord_flip() + ylab('# proteins with PTM') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
      )
    dev.off()
}


sites = psms[c('Biological.set', 'Top.luciphor.PTM')]
multiptms = strsplit(sites[,2], ';', fixed=T)
sites = data.frame(sites=unlist(multiptms), bioset=rep(psms$Biological.set, sapply(multiptms, FUN=length)))
sites = cbind(sites$bioset, stringr::str_split_fixed(sites$sites, ':', 2))
multisites = strsplit(sites[,3], ',', fixed=T)
sites = data.frame(site=unlist(multisites),
                   bioset=rep(sites[,1], sapply(multisites, FUN=length)),
                   ptm=rep(sites[,2], sapply(multisites, FUN=length))
                   )
sites = subset(sites, ptm %in% usedptms)
sites$site = sub('[0-9]+', '', sites$site)
sites$value=1
sites$site = paste(sites$ptm, sites$site)
sitenums = aggregate(value~site+bioset, sites, sum)
svg('psmptmresidues', width=width, height=nrsets+2)
print(ggplot(sitenums) +
    geom_bar(aes(bioset, value, fill=site), stat='identity', position='dodge') +
    coord_flip() + ylab('# PTM residues found') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
)
dev.off()


set_amount_pep$sites = sub('.*_', '', set_amount_pep$Peptide.sequence)
multiptms = strsplit(set_amount_pep$sites, ';', fixed=T)
sites = data.frame(sites=unlist(multiptms), bioset=rep(set_amount_pep$variable, sapply(multiptms, FUN=length)))
sites = cbind(sites$bioset, stringr::str_split_fixed(sites$sites, ':', 2))
multisites = strsplit(sites[,3], ',', fixed=T)
sites = data.frame(site=unlist(multisites),
                   bioset=rep(sites[,1], sapply(multisites, FUN=length)),
                   ptm=rep(sites[,2], sapply(multisites, FUN=length))
                   )
sites = subset(sites, ptm %in% usedptms)
sites$site = sub('[0-9]+', '', sites$site)
sites$value=1
sites$site = paste(sites$ptm, sites$site)
pepsitenums = aggregate(value~site+bioset, sites, sum)
svg('pepptmresidues', width=width, height=nrsets+2)
print(ggplot(pepsitenums) +
    geom_bar(aes(bioset, value, fill=site), stat='identity', position='dodge') +
    coord_flip() + ylab('# PTM residues found') + theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=10), axis.title.y=element_blank(), legend.text=element_text(size=10), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=15))
)
dev.off()
