#!/usr/bin/env python3

import sys
import re
from os import environ

from jinja2 import Template

from create_modfile import get_msgfmods, categorize_mod, parse_cmd_mod, modpos, NON_BLOCKING_MODS


aa_weights_monoiso = { # From ExPASY
        'A': 71.03711,
        'R': 156.10111,
        'N': 114.04293,
        'D': 115.02694,
        'C': 103.00919,
        'E': 129.04259,
        'Q': 128.05858,
        'G': 57.02146,
        'H': 137.05891,
        'I': 113.08406,
        'L': 113.08406,
        'K': 128.09496,
        'M': 131.04049,
        'F': 147.06841,
        'P': 97.05276,
        'S': 87.03203,
        'T': 101.04768,
        'W': 186.07931,
        'Y': 163.06333,
        'V': 99.06841,
        'U': 150.953636,
        'O': 237.147727,
        }

def get_msgf_seq_mass(mass):
    return '{}{}'.format('+' if mass > 0 else '', str(round(mass, 3)))


def add_mods_translationtable(mod, ttab):
    mass = float(mod[0])
    msgf_seq_mass = get_msgf_seq_mass(mass)
    ttab[msgf_seq_mass] = mass
    return msgf_seq_mass


def get_luci_mod(mod):
    if mod[3] == 'N-term':
        residue = '['
    elif mod[3] == 'C-term':
        residue = ']'
    else:
        residue = mod[1]
    return '{} {}'.format(residue, mod[0])

def main():
    psmfile = sys.argv[1]
    template = sys.argv[2]
    modlibfile = sys.argv[3]
    mods = sys.argv[4]
    cmdptms = sys.argv[5].split(';')
    outfile = sys.argv[6]
    ms2tol = environ.get('MS2TOLVALUE')
    ms2toltype = {'ppm': 1, 'Da': 0}[environ.get('MS2TOLTYPE')]


    massconversion_msgf = {}
    fixedmods, varmods = {}, []
    target_mods, decoy_mods, lucifixed, lucivar = [], set(), [], []

    msgfmods = get_msgfmods(modlibfile)
    # Get nontarget variable mods from cmd line, classify as fixed/var
    for cmdmod in mods.split(';'):
        modlines = parse_cmd_mod(cmdmod, msgfmods)
        categorize_mod(modlines, fixedmods, varmods)
    # Now loop variable mods again and create lookup dict for PSMtable mod -> luciphor input mod
    #for cmdmod in mods.split(';'):
    #    modlines = parse_cmd_mod(cmdmod, msgfmods)
    for modlines in fixedmods.values():
        for mod in modlines:
            add_mods_translationtable(mod, massconversion_msgf)
    for mod in varmods:
        realmodmass = mod[0]
        mmass, mres, mfm, mprotpos, mname = mod
        mp = modpos(mod)
        if mp in fixedmods:
            nonblocked_fixed = NON_BLOCKING_MODS[mname] if mname in NON_BLOCKING_MODS else []
            blocked_fixmass = sum([float(x[0]) for x in fixedmods[mp] if x[-1] not in nonblocked_fixed])
            mod[0] = str(round(-(blocked_fixmass - float(mod[0])), 5))
        add_mods_translationtable(mod, massconversion_msgf)

    # Prep fixed mods for luciphor template
    for mod in [x for fm in fixedmods.values() for x in fm]:
        for residue in mod[1]:
            lucimod = get_luci_mod([mod[0], residue, *mod[2:]])
            if lucimod not in lucifixed:
                lucifixed.append(lucimod)
    # Var mods too, and add to mass list to filter PSMs on later (all var mods must be annotated on sequence input)
    varmods_mass = []
    for mod in varmods:
        varmods_mass.append(get_msgf_seq_mass(float(mod[0])))
        for residue in mod[1]:
            lucimod = get_luci_mod([mod[0], residue, *mod[2:]])
            if lucimod not in lucivar:
                lucivar.append(lucimod)

    # Get PTMs from cmd line and prep for template
    ptms_mass = []
    for cmdptm in cmdptms:
        for ptm in parse_cmd_mod(cmdptm, msgfmods):
            # FIXME OPTIONAL competition for multi-residue/line spec, how to know which 
            # modifications can compete?
            ptm = ptm.split(',')
            realmodmass = ptm[0]
            if modpos(ptm) in fixedmods:
                fixmass = sum([float(x[0]) for x in fixedmods[modpos(ptm)]])
                ptm[0] = str(round(-(fixmass - float(ptm[0])), 5))
                print(realmodmass, ptm)
            ptms_mass.append(add_mods_translationtable(ptm, massconversion_msgf))
            for residue in ptm[1]:
                luciptm = get_luci_mod([ptm[0], residue, *ptm[2:]])
                ## TODO
                ## Now we add e.g. -187 if acetyl mod competes with fixed TMT
                ## That does not work, but the below line also doesnt work, which
                ## Puts the actual +42 mass in the luciphor config
                ## Total mass of residue is in the luci input K+42 = 170
                #luciptm = get_luci_mod([realmodmass, residue, *ptm[2:]])
                if luciptm not in target_mods:
                    target_mods.append(luciptm)
                    decoy_mods.add(ptm[0])
        # Neutral loss only for phospho (and possibly glycosylation), add more if we need
        nlosses, decoy_nloss = [], []
        if ptm[4] == 'Phospho':
            nlosses.append('sty -H3PO4 -97.97690')
            decoy_nloss.append('X -H3PO4 -97.07690')


    with open(template) as fp, open('luciphor_config.txt', 'w') as wfp:
        lucitemplate = Template(fp.read())
        wfp.write(lucitemplate.render(
            outfile=outfile,
            fixedmods=lucifixed,
            varmods=lucivar,
            ptms=target_mods,
            ms2tol=ms2tol,
            ms2toltype=ms2toltype,
            dmasses=decoy_mods,
            neutralloss=nlosses,
            decoy_nloss=decoy_nloss
            ))

    # acetyl etc? # FIXME replace double notation 229-187 in PSM table with the actual mass (42)
    # translation table needed...
    # But how to spec in luciphor, it also wants fixed/var/target mods? Does it apply fixed regardless?
    # Or does it know there is competition?

    with open(psmfile) as fp, open('lucipsms', 'w') as wfp:
        header = next(fp).strip('\n').split('\t')
        pepcol = header.index('Peptide')
        spfile = header.index('SpectraFile')
        charge = header.index('Charge')
        scan = header.index('ScanNum')
        evalue = header.index('PSM q-value')
        wfp.write('srcFile\tscanNum\tcharge\tPSMscore\tpeptide\tmodSites')
        for line in fp:
            ptm_in_seq = False
            sites = []
            pep = ''
            start = 0
            line = line.strip('\n').split('\t')
            seq = line[pepcol]
            #for mod in re.finditer('[\+\-0-9]+.[0-9]+', seq):
            for mod in re.finditer('[\+\-0-9]+.[\+\-.0-9]+', seq):
                modtxt = mod.group()
                if modtxt.count('+') + modtxt.count('-') > 1:
                    # multi mod on single residue
                    multimods = re.findall('[\+\-][0-9.]+', modtxt)
                else:
                    multimods = [modtxt]
                modmass = sum([massconversion_msgf[x] for x in multimods])
                if set(multimods).intersection(ptms_mass):
                    ptm_in_seq = True
                    # IF fixed + target -> use the target mass from the msgfmods file
                    # target + target -> 
                # FIXME overlapping mods
    #                if modtxt in ptms_mass:
    #                    ptm_= True
                pep += seq[start:mod.start()]
                if not set(multimods).intersection([*ptms_mass, *varmods_mass]):
                    # Do not annotate fixed mods
                    start = mod.end()
                    continue
                if len(pep) == 0:
                    aamass = modmass
                else:
                    aamass = round(aa_weights_monoiso[pep[-1]] + modmass, 5)
                sites.append([len(pep)- 1, aamass])
                start = mod.end()
            if not ptm_in_seq:
                # only output PSMs with PTM for luciphor
                continue
            if start != mod.endpos:
                pep += seq[start:]
            if sites[0][0] == -1: sites[0][0] = -100
            # TODO add C-terminal mods (rare)
            wfp.write('\n{}\t{}\t{}\t{}\t{}\t{}'.format(line[spfile], line[scan], line[charge], line[evalue], pep, ','.join(['{}={}'.format(x[0], x[1]) for x in sites])))


if __name__ == '__main__':
    main()
