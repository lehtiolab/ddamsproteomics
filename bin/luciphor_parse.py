#!/usr/bin/env python3

import sys
import os
import re
import argparse
from Bio import SeqIO

from luciphor_prep import aa_weights_monoiso, parse_mods_msgf_pep, add_mods_translationtable
from create_modfile import get_msgfmods, parse_cmd_mod, categorize_mod, modpos, NON_BLOCKING_MODS


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
    labileptms = args.labileptms

    # First prepare a residue + PTM weight -> PTM name dict for naming mods
    mass_to_name, ptmmasses = {}, {}
    msgfmods = get_msgfmods(args.modfile)
    for ptmname in labileptms + args.stabileptms + args.mods:
        for modline in msgfmods[ptmname.lower()]:
            modline = modline.split(',')
            if modline[2] == 'fix':
                continue
            modmass = float(modline[0])
            if modline[1] == '*':
                ptmmasses[int(round(modmass, 0))] = ptmname
            else:
                for res in modline[1]:
                    mass = aa_weights_monoiso[res] + modmass
                    ptmmasses[int(round(mass, 0))] = ptmname
            mass_to_name[modmass] = ptmname

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
            modresidues = {x: [] for x in labileptms}
            barepep, start = '', 0
            modpep = line['predictedPep1']
            for x in re.finditer('([A-Z]){0,1}\[([0-9]+)\]', modpep):
                if x.group(1) is not None: # check if residue (or protein N-term)
                    barepep += modpep[start:x.start()+1]
                else:
                    barepep += modpep[start:x.start()]
                start = x.end()
                ptmname = ptmmasses[int(x.group(2))]
                if ptmname in labileptms:
                    modresidues[ptmname].append((x.group(1), len(barepep)))
            # Remove PTMs if not found on PSM
            modresidues = {pn: res for pn, res in modresidues.items() if len(res)}

            barepep += modpep[start:]
            ptm.update({'barepep': barepep, 'modres': modresidues})
            ptm[TOPPTM] = '_'.join(['{}:{}'.format(name, ','.join(['{}{}'.format(x[0], x[1]) for x in resmods])) for 
                    name, resmods in modresidues.items()])
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
        # First see if we have a stabile PTM that if on a fixed mod residue, 
        # because that needs parsing from PSM table (luciphor throws them out)
        fixedmods, varmods, masslookup, varfixmods_masses = {}, [], {}, []
        for modname in args.stabileptms + args.mods:
            modlines = parse_cmd_mod(modname, msgfmods)
            categorize_mod(modlines, fixedmods, varmods)
        for modlines in fixedmods.values():
            for mod in modlines:
                add_mods_translationtable(mod, masslookup)
        varmod_on_fixmod_res = False
        for mod in varmods:
            mmass, mres, mfm, mprotpos, mname = mod
            mp = modpos(mod)
            if mp in fixedmods:
                varmod_on_fixmod_res = True
                nonblocked_fixed = NON_BLOCKING_MODS.get(mname, [])
                blocked_fixmass = sum([float(x[0]) for x in fixedmods[mp] if x[-1] not in nonblocked_fixed])
                mod[0] = str(round(-(blocked_fixmass - float(mod[0])), 6))
                varfixmods_masses.append(add_mods_translationtable(mod, masslookup))
            else:
                add_mods_translationtable(mod, masslookup)
        for ptmname in labileptms:
            for ptmdef in parse_cmd_mod(ptmname, msgfmods):
                # FIXME OPTIONAL competition for multi-residue/line spec, how to know which 
                # modifications can compete?
                ptmdef = ptmdef.split(',')
                if modpos(ptmdef) in fixedmods:
                    fixmass = sum([float(x[0]) for x in fixedmods[modpos(ptmdef)]])
                    ptmdef[0] = str(round(-(fixmass - float(ptmdef[0])), 5))
                add_mods_translationtable(ptmdef, masslookup)

        for psm in psms:
            psm = psm.strip('\n').split('\t')
            psm = {k: v for k,v in zip(psmheader, psm)}
            psmid = '{}.{}.{}.{}'.format(os.path.splitext(psm['SpectraFile'])[0], psm['ScanNum'], psm['ScanNum'], psm['Charge'])
            if psmid in lucptms:
                ptm = lucptms[psmid]
                if varmod_on_fixmod_res:
                    parsedpep = parse_mods_msgf_pep(psm[PEPTIDE], masslookup, varfixmods_masses)
                    if parsedpep:
                        sites = {masslookup[x]: [] for x in varfixmods_masses}
                        for site in parsedpep['sites']:
                            for modmass in site[2]:
                                sites[modmass].append(f'{parsedpep["barepeptide"][site[0]]}{site[0]}')
                        sites = {mass: ','.join(s) for mass, s in sites.items()}
                        stabilefixptms = '_'.join([f'{mass_to_name[mass]}:{s}' for mass, s in sites.items()])
                        ptm[TOPPTM] = f'{ptm[TOPPTM]}_{stabilefixptms}'

                # Get protein site location of mods
                if MASTER_PROTEIN in psm:
                    flankseqs = set()
                    proteins_peploc = {}
                    for p in psm[MASTER_PROTEIN].split(';'):
                        proteins_peploc[p] = [x.start() for x in re.finditer(ptm['barepep'], str(tdb[p].seq))]
                    proteins_loc = {p: [] for p, peplocs in proteins_peploc.items() if len(peplocs)}
                    for p, peplocs in proteins_peploc.items():
                        # peplocs = [4, 120, ..] # unusual to have multiple mathces?
                        for ptmname, ptmlocs in ptm['modres'].items():
                            if ptmname not in labileptms:
                                continue
                            protseq = tdb[p].seq
                            protptms = []
                            for res_loc in ptmlocs:
                                site_protlocs = [res_loc[1] + x for x in peplocs]
                                protlocs = '/'.join([str(x) for x in site_protlocs])
                                protptms.append(f'{res_loc[0]}{protlocs}')
                                flankpos = [(max(x-8, 0) , min(x+7, len(protseq))) for x in site_protlocs]
                                flankseqs.update([str(protseq[x[0]:x[1]]) for x in flankpos])

                            proteins_loc[p].append('{}:{}'.format(ptmname, ','.join(protptms)))
                    psm[MASTER_PROTEIN] = ';'.join(['{}__{}'.format(p, '_'.join(ptmloc)) for p, ptmloc in proteins_loc.items()])
                    psm[FLANKING_SEQS] = ';'.join(flankseqs)
                outpsm = {k: v for k,v in psm.items()}
                outpsm.update(ptm)
                outpsm[SE_PEPTIDE] = outpsm.pop(PEPTIDE)
                outpsm[PEPTIDE] = '{}_{}'.format(ptm['barepep'], ptm[TOPPTM])
                wfp.write('\n{}'.format('\t'.join([outpsm[k] for k in outheader])))


if __name__ == '__main__':
    main()
