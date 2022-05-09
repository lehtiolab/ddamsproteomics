#!/usr/bin/env python3

import sys
import os
import re
import argparse
from Bio import SeqIO

from luciphor_prep import PSM, Mods


# PTM input/output header fields
TOPPTM = 'Top luciphor PTM'
TOPFLR = 'Top PTM FLR'
TOPSCORE = 'Top PTM score'
OTHERPTMS = 'High-scoring PTMs'
SE_PEPTIDE = 'SearchEnginePeptide'
PEPTIDE = 'Peptide'
MASTER_PROTEIN = 'Master protein(s)'
FLANKING_SEQS = 'PTM flanking seq(s)'
PTMFIELDS = [SE_PEPTIDE, TOPPTM, TOPSCORE, TOPFLR, OTHERPTMS, FLANKING_SEQS]

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--minscore', type=float)
    parser.add_argument('-o', dest='outfile')
    parser.add_argument('--modfile')
    parser.add_argument('--labileptms', nargs='+', default=[])
    parser.add_argument('--stabileptms', nargs='+', default=[])
    parser.add_argument('--mods', nargs='+', default=[])
    parser.add_argument('--fasta')
    args = parser.parse_args(sys.argv[1:])

    minscore_high = args.minscore
    labileptms = [x.lower() for x in args.labileptms]
    stabileptms = [x.lower() for x in args.stabileptms]
    mods = [x.lower() for x in args.mods]

    # First prepare a residue + PTM weight -> PTM name dict for naming mods
    msgfmods = Mods()
    msgfmods.parse_msgf_modfile(args.modfile, labileptms, [*stabileptms, *mods])
    luci_modmap = msgfmods.lucimass_mod_dict()
    msgf_mod_map = msgfmods.msgfmass_mod_dict()

    # load sequences
    tdb = SeqIO.index(args.fasta, 'fasta')
    # Now go through scores, luciphor and PSM table
    with open('all_scores.debug') as scorefp, open('luciphor.out') as fp, open('psms') as psms, open(args.outfile, 'w') as wfp:
        header = next(fp).strip('\n').split('\t')
        scoreheader = next(scorefp).strip('\n').split('\t')
        psmheader = next(psms).strip('\n').split('\t')
        scorepep = {'specId': False}
        outheader = psmheader + PTMFIELDS
        wfp.write('\t'.join(outheader))
        lucptms = {}
        for line in fp:
            line = {k: v for k,v in zip(header, line.strip('\n').split('\t'))}
            lucipsm = PSM()
            lucipsm.parse_luciphor_peptide(line, luci_modmap, labileptms, stabileptms)

            # Get other highscoring permutations
            extrapeps = []
            if scorepep['specId'] == lucipsm.lucispecid:
                lucipsm.parse_luciphor_scores(scorepep, minscore_high)
            for scorepep in scorefp:
                scorepep = scorepep.strip('\n').split('\t')
                scorepep = {k: v for k,v in zip(scoreheader, scorepep)}
                if int(scorepep['isDecoy']):
                    continue
                if scorepep['specId'] != lucipsm.lucispecid:
                    break
                lucipsm.parse_luciphor_scores(scorepep, minscore_high)
            lucptms[lucipsm.lucispecid] = lucipsm

        # With that parsed, take the MSGF PSMs
        # First see if we have a stabile PTM that is on a fixed mod residue, 
        # because that needs parsing from PSM table (luciphor throws them out)
        for psm in psms:
            psm = psm.strip('\n').split('\t')
            psm = {k: v for k,v in zip(psmheader, psm)}
            psmid = '{}.{}.{}.{}'.format(os.path.splitext(psm['SpectraFile'])[0], psm['ScanNum'], psm['ScanNum'], psm['Charge'])
            if psmid in lucptms:
                luciptm = lucptms[psmid]
                ptm = {
                    TOPFLR: luciptm.top_flr,
                    TOPSCORE: luciptm.top_score,
                    OTHERPTMS: luciptm.format_alt_ptm_locs(),
                    }
                if msgfmods.has_varmods_on_fixmod_residues:
                    # If there are PTMs which are on fixed mod residues, we need to parse
                    # those from the residues from sequences in the PSM table because luciphor
                    # does not report those residues because they have a fixed mod
                    msgfpsm = PSM()
                    msgfpsm.parse_msgf_peptide(psm[PEPTIDE], msgf_mod_map,
                            labileptms, stabileptms)
                    luciptm.add_ptms_from_psm(msgfpsm.mods)
                ptm[TOPPTM] = luciptm.topptm_output()

                # Get protein site location of mods
                if MASTER_PROTEIN in psm:
                    annotate_protein_and_flanks(psm, luciptm, tdb, [*labileptms, *stabileptms])
                outpsm = {k: v for k,v in psm.items()}
                outpsm.update(ptm)
                outpsm[SE_PEPTIDE] = outpsm.pop(PEPTIDE)
                outpsm[PEPTIDE] = '{}_{}'.format(luciptm.sequence, ptm[TOPPTM])
                wfp.write('\n{}'.format('\t'.join([outpsm[k] for k in outheader])))


def annotate_protein_and_flanks(psm, ptmpsm, tdb, ptmnames):
    flankseqs = set()
    proteins_peploc = {}
    for p in psm[MASTER_PROTEIN].split(';'):
        proteins_peploc[p] = [x.start() for x in re.finditer(ptmpsm.sequence, str(tdb[p].seq))]
    proteins_loc = {p: [] for p, peplocs in proteins_peploc.items() if len(peplocs)}
    for p, peplocs in proteins_peploc.items():
        # peplocs = [4, 120, ..] # unusual to have multiple mathces?
        for ptm in ptmpsm.mods:
            if ptm['name_lower'] not in ptmnames:
                continue
            protseq = tdb[p].seq
            protptms = []
            site_protlocs = [ptm['site'][1] + x for x in peplocs]
            protlocs = '/'.join([str(x) for x in site_protlocs])
            protptms.append(f'{ptm["site"][0]}{protlocs}')
            flankpos = [(max(x-7, 0) , min(x+8, len(protseq))) for x in site_protlocs]
            flankseqs.update([str(protseq[x[0]:x[1]]) for x in flankpos])
            proteins_loc[p].append('{}:{}'.format(ptm['name'], ','.join(protptms)))
    psm[MASTER_PROTEIN] = ';'.join(['{}__{}'.format(p, '_'.join(ptmloc)) for p, ptmloc in proteins_loc.items()])
    psm[FLANKING_SEQS] = ';'.join(flankseqs)


if __name__ == '__main__':
    main()
