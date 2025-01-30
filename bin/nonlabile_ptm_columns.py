#!/usr/bin/env python3

import sys
import re
from Bio import SeqIO
import argparse

from luciphor_prep import Mods, PSM, aa_weights_monoiso
import luciphor_parse as lucp


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o', dest='outfile')
    parser.add_argument('--psms')
    parser.add_argument('--modfile')
    parser.add_argument('--labileptms', nargs='+', default=[])
    parser.add_argument('--stabileptms', nargs='+', default=[])
    parser.add_argument('--mods', nargs='+', default=[])
    parser.add_argument('--fasta')
    parser.add_argument('--searchengine', default='msgf', type=str)
    args = parser.parse_args(sys.argv[1:])

    ptms = [x.lower() for x in args.stabileptms]
    locptms = [x.lower() for x in args.labileptms]
    mods = [x.lower() for x in args.mods]

    msgfmods = Mods()
    msgfmods.parse_msgf_modfile(args.modfile, [*locptms, *ptms, *mods])
    msgf_mod_map = msgfmods.msgfmass_mod_dict(args.searchengine)
    tdb = SeqIO.index(args.fasta, 'fasta')
    fields = lucp.get_headerfields(args.searchengine)
    print(fields)
    with open(args.psms) as fp, open(args.outfile, 'w') as wfp:
        header = next(fp).strip('\n').split('\t')
        outheader = header + fields['PTMFIELDS']
        wfp.write('\t'.join(outheader))

        for psm in fp:
            psm = psm.strip('\n').split('\t')
            psm = {k: v for k,v in zip(header, psm)}
            psm.update({x: 'NA' for x in fields['PTMFIELDS']})
            ptmpsm = PSM()
            ptmpsm.parse_search_peptide(psm[fields['PEPTIDE']], msgf_mod_map, locptms, ptms,
                    args.searchengine)
            if ptmpsm.has_stableptms() and not ptmpsm.has_labileptms():
                psm[fields['TOPPTM']] = ptmpsm.topptm_output()
                # Add protein location annotation
                if fields['MASTER_PROTEIN'] in psm:
                    lucp.annotate_protein_and_flanks(psm, ptmpsm, tdb, ptms, fields)
                psm[fields['SE_PEPTIDE']] = psm.pop(fields['PEPTIDE'])
                psm[fields['PEPTIDE']] = '{}_{}'.format(ptmpsm.sequence, psm[fields['TOPPTM']])
                wfp.write('\n{}'.format('\t'.join([psm[k] for k in outheader])))


if __name__ == '__main__':
    main()

