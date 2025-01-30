#!/usr/bin/env python3

import sys
import re
from os import environ
import argparse

from mods import Mods, aa_weights_monoiso


class PSM: 
    '''A PSM class containing mods, scores, FLR, etc
    Mods are defined as dicts, and apart from other info like aa, type, mass, etc,
    they contain two keys, site_lucin, site_report which are zero resp. one-based
    residue indices for luciphor input (zero based) and reporting to PSM tables
    (one based)
    '''

    def __init__(self):
        self.mods = []
        self.top_flr = False
        self.top_score = False
        self.lucispecid = False
        self.alt_ptm_locs = []
        self.sequence = False
        self.seq_in_scorepep_fmt = False


    def get_modtype(self, mod, labileptmnames, stableptmnames):
        if not mod['var']:
            mtype = 'fixed'
        elif mod['name_lower'] in labileptmnames:
            mtype = 'labile'
        elif mod['name_lower'] in stableptmnames:
            mtype = 'stable'
        else:
            mtype = 'variable'
        return mtype

    def get_mod_dict(self, residue, sitenum, modptm, labileptmnames, stableptmnames):
        return {'aa': residue, 'site_lucin': sitenum, 'site_report': sitenum + 1,
                'type': self.get_modtype(modptm, labileptmnames, stableptmnames),
                'mass': modptm['mass'], 'name': modptm['name'],
                'name_lower': modptm['name_lower'], 'adjusted_mass': modptm['adjusted_mass']}

    def parse_search_peptide(self, searchseq, msgf_mods, labileptmnames, stableptmnames, search_engine):
        self.mods = []
        barepep = ''
        start = 0
        if search_engine == 'sage':
            regex = r'([A-Z]){0,1}([\[\]0-9\.+\-]+)'
        else:
            regex = r'([A-Z]){0,1}([0-9\.+\-]+)'
        for x in re.finditer(regex, searchseq):
            if x.group(1) is not None:
                # mod is on a residue
                barepep = f'{barepep}{searchseq[start:x.start()+1]}'
                residue = barepep[-1]
                sitenum = len(barepep) - 1
            else:
                # mod is on protein N-term
                residue = '['
                sitenum = -100
            # TODO cterm = 100, ']'
            start = x.end()
            for mass in re.findall(r'[\+\-][0-9.]+', x.group(2)):
                try:
                    mod = msgf_mods[float(mass)][0] # only take first, contains enough info
                except KeyError:
                    if search_engine == 'sage':
                        # Sage turns 57.021464 into 57.021465 (carbamidomethyl)
                        # Not sure about others or why, but try to counteract here
                        alt_sagemass = float(f'{float(mass):.4f}')
                        mod = msgf_mods[alt_sagemass][0]
                    else:
                        raise
                self.mods.append(self.get_mod_dict(residue, sitenum, mod, labileptmnames,
                    stableptmnames))
        self.sequence = f'{barepep}{searchseq[start:]}'

    def parse_luciphor_peptide(self, luciline, ptms_map, labileptmnames, stableptmnames):
        '''From a luciphor sequence, create a peptide with PTMs
        ptms_map = {f'{residue}int(79 + mass_S/T/Y)': {'name': Phospho, etc}
        '''
        self.top_flr = luciline['globalFLR']
        self.top_score = luciline['pep1score']
        self.lucispecid = luciline['specId']
        self.mods = []
        barepep, start = '', 0
        modpep = luciline['predictedPep1']
        for x in re.finditer(r'([A-Z]){0,1}\[([0-9]+)\]', modpep):
            if x.group(1) is not None: # check if residue (or protein N-term)
                barepep += modpep[start:x.start()+1]
            start = x.end()
            ptm = ptms_map[f'{x.group(1)}{int(x.group(2))}']
            if ptm['name_lower'] in labileptmnames:
                sitenum = len(barepep) - 1 if len(barepep) else -100
                residue = barepep[-1] if len(barepep) else '['
                self.mods.append(self.get_mod_dict(residue, sitenum, ptm, labileptmnames,
                    stableptmnames))
        self.sequence = f'{barepep}{modpep[start:]}'
        self.seq_in_scorepep_fmt = re.sub(r'([A-Z])\[[0-9]+\]', lambda x: x.group(1).lower(), modpep)

    def parse_luciphor_scores(self, scorepep, minscore):
        permut = scorepep['curPermutation']
        if permut != self.seq_in_scorepep_fmt and float(scorepep['score']) > minscore:
            self.alt_ptm_locs.append([f'{x.group()}{x.start() + 1}:{scorepep["score"]}'
                for x in re.finditer('[a-z]', permut)])

    def format_alt_ptm_locs(self):
        alt_locs = [','.join(x).upper() for x in self.alt_ptm_locs]
        return ';'.join(alt_locs) if len(alt_locs) else 'NA'

    def has_labileptms(self):
        return any(m['type'] == 'labile' for m in self.mods)
    
    def has_stableptms(self):
        return any(m['type'] == 'stable' for m in self.mods)

    def luciphor_input_sites(self):
        lucimods = []
        for m in self.mods:
            if m['type'] != 'fixed':
                lucimods.append((m['site_lucin'], str(m['mass'] + aa_weights_monoiso[m['aa']])))
        return ','.join([f'{x[0]}={x[1]}' for x in lucimods])

    def add_ptms_from_psm(self, psmmods):
        existing_mods = {m['name']: m for m in self.mods}
        for psmmod in psmmods:
            if psmmod['name'] not in existing_mods:
                self.mods.append(psmmod)

    def topptm_output(self):
        ptmsites = {}
        output_types = {'labile', 'stable'}
        for ptm in self.mods:
            if ptm['type'] not in output_types:
                continue
            site = f'{ptm["aa"]}{ptm["site_report"]}'
            try:
                ptmsites[ptm['name']].append(site)
            except KeyError:
                ptmsites[ptm['name']] = [site]
        return '_'.join([f'{p}:{",".join(s)}' for p, s in ptmsites.items()])


def create_msgf_mod_lookup():
    lookup = {}



def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--psmfile')
    parser.add_argument('--template')
    parser.add_argument('--lucipsms', dest='lucipsms')
    parser.add_argument('--modfile')
    parser.add_argument('--labileptms', nargs='+', default=[])
    parser.add_argument('--mods', nargs='+', default=[])
    parser.add_argument('--searchengine', default='msgf', type=str)

    args = parser.parse_args(sys.argv[1:])

    labileptms = [x.lower() for x in args.labileptms]
    othermods = [x.lower() for x in args.mods]

    msgfmods = Mods()
    msgfmods.parse_msgf_modfile(args.modfile, [*args.labileptms, *args.mods])
    # Prep fixed mods for luciphor template
    lucifixed = []
    for mod in msgfmods.fixedmods:
        lucifixed.append(msgfmods.get_luci_input_mod_line(mod))

    # Var mods too, and add to mass list to filter PSMs on later (all var mods must be annotated on sequence input)
    lucivar = []
    for mod in msgfmods.varmods:
        if mod['name_lower'] in labileptms:
            continue
        lucivar.append(msgfmods.get_luci_input_mod_line(mod))

    # Get PTMs from cmd line and prep for template
    # Luciphor does not work when specifying PTMs on same residue as fixed mod
    # e.g. TMT and something else, because it throws out the residues with fixed mods
    # https://github.com/dfermin/lucXor/issues/11
    target_mods, decoy_mods = [], set()
    nlosses, decoy_nloss = [], []
    for mod in msgfmods.varmods:
        if mod['name_lower'] in labileptms:
            target_mods.append(msgfmods.get_luci_input_mod_line(mod))
            decoy_mods.add(msgfmods.get_mass_or_adj(mod))
        if mod['name'] == 'Phospho':
            nlosses.append('sty -H3PO4 -97.97690')
            decoy_nloss.append('X -H3PO4 -97.07690')
            
    with open(args.template) as fp, open('luciphor_config.txt', 'w') as wfp:
        wfp.write(fp.read())
        for fixmod in lucifixed:
            wfp.write(f'\nFIXED_MOD = {fixmod}')
        for varmod in lucivar:
            wfp.write(f'\nVAR_MOD = {varmod}')
        for ptmmod in target_mods:
            wfp.write(f'\nTARGET_MOD = {ptmmod}')
        for nloss in nlosses:
            wfp.write(f'\nNL = {nloss}')
        for dmass in decoy_mods:
            wfp.write(f'\nDECOY_MASS = {dmass}')
        for dnl in decoy_nloss:
            wfp.write(f'\nDECOY_NL = {dnl}')

    # acetyl etc? # FIXME replace double notation 229-187 in PSM table with the actual mass (42)
    # translation table needed...
    # But how to spec in luciphor, it also wants fixed/var/target mods? Does it apply fixed regardless?

    msgf_mod_map = msgfmods.msgfmass_mod_dict(args.searchengine)
    colummap = {'sage': {
        'pep': 'peptide',
        'fn': 'filename',
        'ch': 'charge',
        'scan': 'scannr',
        'value': 'PSM q-value',
        }, 'msgf': {
        'pep': 'Peptide',
        'fn': '#SpecFile',
        'ch': 'Charge',
        'scan': 'ScanNum',
        'value': 'PSM q-value',
            }
        }
            
    with open(args.psmfile) as fp, open(args.lucipsms, 'w') as wfp:
        header = next(fp).strip('\n').split('\t')
        pepcol = header.index(colummap[args.searchengine]['pep'])
        spfile = header.index(colummap[args.searchengine]['fn'])
        charge = header.index(colummap[args.searchengine]['ch'])
        scancol = header.index(colummap[args.searchengine]['scan'])
        evalue = header.index(colummap[args.searchengine]['value'])
        wfp.write('srcFile\tscanNum\tcharge\tPSMscore\tpeptide\tmodSites')
        for line in fp:
            line = line.strip('\n').split('\t')
            psm = PSM()
            psm.parse_search_peptide(line[pepcol], msgf_mod_map, labileptms, othermods, args.searchengine)
            if args.searchengine == 'sage':
                scan = re.sub('.*scan=', '', line[scancol])
            else:
                scan = line[scancol]
            # TODO add C-terminal mods (rare)
            if psm.has_labileptms():
                wfp.write('\n{}\t{}\t{}\t{}\t{}\t{}'.format(line[spfile], scan, line[charge], line[evalue], psm.sequence, psm.luciphor_input_sites()))


if __name__ == '__main__':
    main()
