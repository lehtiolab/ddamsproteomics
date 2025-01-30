#!/usr/bin/env python3

import sys
import os
import re
import argparse
from Bio import SeqIO

from luciphor_prep import PSM, Mods


# PTM input/output header fields
def get_headerfields(search_engine):
    fields = {'TOPPTM': 'Top luciphor PTM', 'TOPFLR': 'Top PTM FLR', 'TOPSCORE': 'Top PTM score',
            'OTHERPTMS': 'High-scoring PTMs', 'SE_PEPTIDE': 'SearchEnginePeptide',
            'MASTER_PROTEIN': 'Master protein(s)', 'FLANKING_SEQS': 'PTM flanking seq(s)',

            'PEPTIDE': {'sage': 'peptide', 'msgf': 'Peptide'}[search_engine],
            'CHARGE': {'sage': 'charge', 'msgf': 'Charge'}[search_engine],
            'FILENAME': {'sage': 'filename', 'msgf': '#SpecFile'}[search_engine],
            'SCANNR': {'sage': 'scannr', 'msgf': 'ScanNum'}[search_engine],
            }
    fields['PTMFIELDS'] = [fields['SE_PEPTIDE'], fields['TOPPTM'], fields['TOPSCORE'], fields['TOPFLR'], fields['OTHERPTMS'], fields['FLANKING_SEQS']]
    return fields


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--minscore', type=float)
    parser.add_argument('-o', dest='outfile')
    parser.add_argument('--luci_in', dest='lucioutput')
    parser.add_argument('--luci_scores', dest='luciscores')
    parser.add_argument('--psms', dest='psms')
    parser.add_argument('--modfile')
    parser.add_argument('--labileptms', nargs='+', default=[])
    parser.add_argument('--stabileptms', nargs='+', default=[])
    parser.add_argument('--mods', nargs='+', default=[])
    parser.add_argument('--fasta')
    parser.add_argument('--searchengine', default='msgf', type=str)
    args = parser.parse_args(sys.argv[1:])

    minscore_high = args.minscore
    labileptms = [x.lower() for x in args.labileptms]
    stabileptms = [x.lower() for x in args.stabileptms]
    mods = [x.lower() for x in args.mods]
    fields = get_headerfields(args.searchengine)

    # First prepare a residue + PTM weight -> PTM name dict for naming mods
    msgfmods = Mods()
    msgfmods.parse_msgf_modfile(args.modfile, [*args.labileptms, *args.stabileptms, *args.mods])
    luci_modmap = msgfmods.lucimass_mod_dict()
    msgf_mod_map = msgfmods.msgfmass_mod_dict(args.searchengine)

    # load sequences
    tdb = SeqIO.index(args.fasta, 'fasta')
    # Now go through scores, luciphor and PSM table
    # TODO integrate luciphor and detect crash (exit 0, stderr msg), should fail with 1
    lucptms = {}
    scorepep = {'specId': False}
    try:
        with open(args.luciscores) as scorefp, open(args.lucioutput) as fp:
            header = next(fp).strip('\n').split('\t')
            scoreheader = next(scorefp).strip('\n').split('\t')
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
    except FileNotFoundError:
        pass

    # With that parsed, take the MSGF PSMs
    # First see if we have a stabile PTM that is on a fixed mod residue, 
    # because that needs parsing from PSM table (luciphor throws them out)
    with open(args.psms) as psms, open(args.outfile, 'w') as wfp:
        psmheader = next(psms).strip('\n').split('\t')
        outheader = psmheader + fields['PTMFIELDS']
        wfp.write('\t'.join(outheader))
        for psm in psms:
            psm = psm.strip('\n').split('\t')
            psm = {k: v for k,v in zip(psmheader, psm)}
            msgfpsm = PSM()
            msgfpsm.parse_search_peptide(psm[fields['PEPTIDE']], msgf_mod_map,
                    labileptms, stabileptms, args.searchengine)
            if not msgfpsm.has_labileptms():
                continue
            print(psm)
            if args.searchengine == 'sage':
                scan = re.sub('.*scan=', '', psm[fields['SCANNR']])
            else:
                scan = psm[fields['SCANNR']]
            psmid = '{}.{}.{}.{}'.format(os.path.splitext(psm[fields['FILENAME']])[0], scan, scan, psm[fields['CHARGE']])
            if psmid in lucptms:
                luciptm = lucptms[psmid]
                ptm = {
                    fields['TOPFLR']: luciptm.top_flr,
                    fields['TOPSCORE']: luciptm.top_score,
                    fields['OTHERPTMS']: luciptm.format_alt_ptm_locs(),
                    }
            else:
                luciptm = False
                ptm = {
                    fields['TOPFLR']: 'NA',
                    fields['TOPSCORE']: 'NA',
                    fields['OTHERPTMS']: 'NA',
                    }
            # If there are PTMs which are on fixed mod residues, we need to parse
            # those from the residues from sequences in the PSM table because luciphor
            # does not report those residues because they have a fixed mod
            # Also need to parse in case no luci PTM exist
            if not luciptm:
                luciptm = msgfpsm
            elif msgfmods.has_varmods_on_fixmod_residues:
                luciptm.add_ptms_from_psm(msgfpsm.mods)
            ptm[fields['TOPPTM']] = luciptm.topptm_output()

            # Get protein site location of mods
            if fields['MASTER_PROTEIN'] in psm:
                annotate_protein_and_flanks(psm, luciptm, tdb, [*labileptms, *stabileptms], fields)
            outpsm = {k: v for k,v in psm.items()}
            outpsm.update(ptm)
            outpsm[fields['SE_PEPTIDE']] = outpsm.pop(fields['PEPTIDE'])
            outpsm[fields['PEPTIDE']] = '{}_{}'.format(luciptm.sequence, ptm[fields['TOPPTM']])
            wfp.write('\n{}'.format('\t'.join([outpsm[k] for k in outheader])))


def annotate_protein_and_flanks(psm, ptmpsm, tdb, ptmnames, fields):
    flankseqs = set()
    proteins_peploc = {}
    for p in psm[fields['MASTER_PROTEIN']].split(';'):
        proteins_peploc[p] = [x.start() for x in re.finditer(ptmpsm.sequence, str(tdb[p].seq))]
    proteins_loc = {p: [] for p, peplocs in proteins_peploc.items() if len(peplocs)}
    for p, peplocs in proteins_peploc.items():
        # peplocs = [4, 120, ..] # unusual to have multiple mathces?
        for ptm in ptmpsm.mods:
            if ptm['name_lower'] not in ptmnames:
                continue
            protseq = tdb[p].seq
            protptms = []
            site_protlocs = [ptm['site_report'] + x for x in peplocs]
            protlocs = '/'.join([str(x) for x in site_protlocs])
            protptms.append(f'{ptm["aa"]}{protlocs}')
            flankpos = [(max(x-7, 0) , min(x+8, len(protseq))) for x in site_protlocs]
            flankseqs.update([str(protseq[x[0]:x[1]]) for x in flankpos])
            proteins_loc[p].append('{}:{}'.format(ptm['name'], ','.join(protptms)))
    psm[fields['MASTER_PROTEIN']] = ';'.join(['{}__{}'.format(p, '_'.join(ptmloc)) for p, ptmloc in proteins_loc.items()])
    psm[fields['FLANKING_SEQS']] = ';'.join(flankseqs)


if __name__ == '__main__':
    main()
