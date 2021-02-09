#!/usr/bin/env python3

import sys
import re
from Bio import SeqIO

from luciphor_prep import aa_weights_monoiso
import luciphor_parse as lucp
from create_modfile import get_msgfmods, modpos, categorize_mod, parse_cmd_mod


def main():
    inpsms = sys.argv[1]
    outfile = sys.argv[2]
    modfile = sys.argv[3]
    ptms = sys.argv[4].split(';') if sys.argv[4] else []
    locptms = sys.argv[5].split(';') if sys.argv[5] else []
    mods = sys.argv[6].split(';') if sys.argv[6] else []
    fasta = sys.argv[7]

    msgfmods = get_msgfmods(modfile)
    fixedmods = {}
    varmods = []
    for passedmod in ptms + locptms + mods:
        mod = parse_cmd_mod(passedmod, msgfmods)
        try:
            categorize_mod(mod, fixedmods, varmods)
        except Exception:
            sys.stderr.write('Could not identify modification "{}", use one of [{}]\n'.format(passedmod, ', '.join(msgfmods.keys())))
            sys.exit(1)
    varmods = {x[4]: x for x in varmods}
    ptmmasses = {}
    for ptmname in ptms + locptms:
        for res in varmods[ptmname][1]:
            modline = varmods[ptmname][:]
            modline = [float(modline[0]), res, modline[2], modline[3]]
            if modpos(modline) in fixedmods:
                fixmass = sum([float(x[0]) for x in fixedmods[modpos(modline)]])
                ptmmasses[str(round(-(fixmass - modline[0]), 3))] = ptmname
            else:
                ptmmasses[str(round(modline[0], 3))] = ptmname
    tdb = SeqIO.index(fasta, 'fasta')
    with open(inpsms) as fp, open(outfile, 'w') as wfp:
        header = next(fp).strip('\n').split('\t')
        outheader = header + lucp.PTMFIELDS
        wfp.write('\t'.join(outheader))

        for psm in fp:
            psm = psm.strip('\n').split('\t')
            psm = {k: v for k,v in zip(header, psm)}
            proteins = psm[lucp.MASTER_PROTEIN].split(';')
            psm.update({x: 'NA' for x in lucp.PTMFIELDS})
            modresidues = {x: [] for x in ptms + locptms}
            barepep, start = '', 0
            for x in re.finditer('([A-Z]){0,1}[+-]([0-9\.+\-]+)', psm['Peptide']):
                if x.group(1) is not None:
                    # mod is on a residue
                    barepep += psm['Peptide'][start:x.start()+1]
                else:
                    # mod is on protein N-term
                    barepep += psm['Peptide'][start:x.start()]
                start = x.end()
                for mass in re.sub('([+-])', ' \\1', x.group(2)).replace('+', '').split(' '):
                    if mass in ptmmasses:
                        modresidues[ptmmasses[mass]].append((x.group(1), len(barepep)))
            modsfound = set([k for k,v in modresidues.items() if len(v)])
            if not modsfound.intersection(ptms) or modsfound.intersection(locptms):
                # remove peptides without PTM and those with labile PTM (the latter will
                # be treated by luciphor)
                continue
            barepep += psm['Peptide'][start:]
            #ptm.update({'modres': modresidues})
            psm[lucp.TOPPTM] = ';'.join(['{}:{}'.format(name, ','.join(['{}{}'.format(x[0], x[1]) for x in resmods])) for 
                    name, resmods in modresidues.items() if len(resmods)])
            #if psmid in lucptms:
            #    ptm = lucptms[psmid]
            proteins = {p: tdb[p].seq.find(barepep) for p in proteins}
            proteins_loc = {p: [] for p, peploc in proteins.items() if peploc > -1}
            for p, peploc in proteins.items():
                for ptmname, ptmlocs in modresidues.items():
                    if ptmname not in ptms:
                        continue
                    protptms = []
                    for res_loc in ptmlocs:
                        protptms.append('{}{}'.format(res_loc[0], res_loc[1] + peploc))
                        proteins_loc[p].append('{}_{}'.format(ptmname, ','.join(protptms)))
            psm[lucp.PROTEIN] = ';'.join(['{}:{}'.format(p, ':'.join(ptmloc)) for p, ptmloc in proteins_loc.items()])
            #outpsm = {k: v for k,v in psm.items()}
            #outpsm.update(psm)
            assert re.sub('[0-9.\[\]+-]', '', psm['Peptide']) == barepep
            psm[lucp.SE_PEPTIDE] = psm.pop(lucp.PEPTIDE)
            psm[lucp.PEPTIDE] = '{}_{}'.format(barepep, psm[lucp.TOPPTM])
            wfp.write('\n{}'.format('\t'.join([psm[k] for k in outheader])))


if __name__ == '__main__':
    main()

