#!/usr/bin/env python3

import sys
import os
import re
from Bio import SeqIO

from luciphor_prep import aa_weights_monoiso
from create_modfile import get_msgfmods


# PTM input/output header fields
TOPPTM = 'Top luciphor PTM'
TOPFLR = 'Top PTM FLR'
TOPSCORE = 'Top PTM score'
OTHERPTMS = 'High-scoring PTMs'
SE_PEPTIDE = 'SearchEnginePeptide'
PEPTIDE = 'Peptide'
MASTER_PROTEIN = 'Master protein(s)'
PTMFIELDS = [SE_PEPTIDE, TOPPTM, TOPSCORE, TOPFLR, OTHERPTMS]

def main():
    minscore_high = float(sys.argv[1])
    outfile = sys.argv[2]
    modfile = sys.argv[3]
    ptms = sys.argv[4].split(';')
    mods = sys.argv[5].split(';')
    fasta = sys.argv[6]

    ptmmods = ptms + mods
    # First prepare a residue + PTM weight -> PTM name dict for naming mods
    ptmmasses = {}
    msgfmods = get_msgfmods(modfile)
    for ptmname in ptmmods:
        for modline in msgfmods[ptmname.lower()]:
            modline = modline.split(',')
            if modline[2] == 'fix':
                continue
            if modline[1] == '*':
                mass = float(modline[0])
                ptmmasses[int(round(mass, 0))] = ptmname
            else:
                for res in modline[1]:
                    mass = aa_weights_monoiso[res] + float(modline[0])
                    ptmmasses[int(round(mass, 0))] = ptmname
    # load sequences
    tdb = SeqIO.index(fasta, 'fasta')
    # Now go through scores, luciphor and PSM table
    with open('all_scores.debug') as scorefp, open('luciphor.out') as fp, open('psms') as psms, open(outfile, 'w') as wfp:
        header = next(fp).strip('\n').split('\t')
        scoreheader = next(scorefp).strip('\n').split('\t')
        psmheader = next(psms).strip('\n').split('\t')
        scorepep = {'specId': False}
        outheader = psmheader + PTMFIELDS
        wfp.write('\t'.join(outheader))
        lucptms = {}
        for line in fp:
            out = {}
            line = line.strip('\n').split('\t')
            line = {k: v for k,v in zip(header, line)}
            specid = line['specId']
            ptm = {
                TOPFLR: line['globalFLR'],
                TOPSCORE: line['pep1score'],
                OTHERPTMS: 'NA',
                }
            # match the modified residues and group:
            modresidues = {x: [] for x in ptmmods}
            barepep, start = '', 0
            modpep = line['predictedPep1']
            for x in re.finditer('([A-Z]){0,1}\[([0-9]+)\]', modpep):
                if x.group(1) is not None: # check if residue (or protein N-term)
                    barepep += modpep[start:x.start()+1]
                else:
                    barepep += modpep[start:x.start()]
                start = x.end()
                modresidues[ptmmasses[int(x.group(2))]].append((x.group(1), len(barepep)))
            barepep += modpep[start:]
            ptm.update({'barepep': barepep, 'modres': modresidues})
            ptm[TOPPTM] = ';'.join(['{}:{}'.format(name, ','.join(['{}{}'.format(x[0], x[1]) for x in resmods])) for 
                    name, resmods in modresidues.items() if len(resmods)])
            # Get other highscoring permutations
            extrapeps = []
            if scorepep['specId'] == specid:
                line['likeScoredPep'] = re.sub(r'([A-Z])\[[0-9]+\]', lambda x: x.group(1).lower(), line['predictedPep1'])
                if scorepep['curPermutation'] != line['likeScoredPep'] and float(scorepep['score']) > minscore_high:
                    extrapeps.append('{}:{}'.format(','.join(['{}{}'.format(x.group().upper(), x.start() + 1)
                                for x in re.finditer('[a-z]', scorepep['curPermutation'])]),
                                scorepep['score']))
            for scorepep in scorefp:
                scorepep = scorepep.strip('\n').split('\t')
                scorepep = {k: v for k,v in zip(scoreheader, scorepep)}
                if int(scorepep['isDecoy']):
                    continue
                if scorepep['specId'] != specid:
                    break
                line['likeScoredPep'] = re.sub(r'([A-Z])\[[0-9]+\]', lambda x: x.group(1).lower(), line['predictedPep1'])
                if scorepep['curPermutation'] != line['likeScoredPep'] and float(scorepep['score']) > minscore_high:
                    extrapeps.append('{}:{}'.format(','.join(['{}{}'.format(x.group().upper(), x.start() + 1)
                            for x in re.finditer('[a-z]', scorepep['curPermutation'])]),
                            scorepep['score']))
                if len(extrapeps):
                    ptm[OTHERPTMS] = ';'.join(extrapeps)
            lucptms[specid] = ptm
        # With that parsed, take the PSMs
        for psm in psms:
            psm = psm.strip('\n').split('\t')
            psm = {k: v for k,v in zip(psmheader, psm)}
            psmid = '{}.{}.{}.{}'.format(os.path.splitext(psm['SpectraFile'])[0], psm['ScanNum'], psm['ScanNum'], psm['Charge'])
            if psmid in lucptms:
                ptm = lucptms[psmid]
                # Get protein site location of mod
                if MASTER_PROTEIN in psm:
                    proteins = psm[MASTER_PROTEIN].split(';')
                    proteins = {p: tdb[p].seq.find(ptm['barepep']) for p in proteins}
                    proteins_loc = {p: [] for p, peploc in proteins.items() if peploc > -1}
                    for p, peploc in proteins.items():
                        for ptmname, ptmlocs in ptm['modres'].items():
                            if ptmname not in ptms:
                                continue
                            protptms = []
                            for res_loc in ptmlocs:
                                protptms.append('{}{}'.format(res_loc[0], res_loc[1] + peploc))
                                proteins_loc[p].append('{}_{}'.format(ptmname, ','.join(protptms)))
                    psm[MASTER_PROTEIN] = ';'.join(['{}:{}'.format(p, ':'.join(ptmloc)) for p, ptmloc in proteins_loc.items()])
                outpsm = {k: v for k,v in psm.items()}
                outpsm.update(ptm)
                outpsm[SE_PEPTIDE] = outpsm.pop(PEPTIDE)
                outpsm[PEPTIDE] = '{}_{}'.format(ptm['barepep'], ptm[TOPPTM])
                wfp.write('\n{}'.format('\t'.join([outpsm[k] for k in outheader])))


if __name__ == '__main__':
    main()
