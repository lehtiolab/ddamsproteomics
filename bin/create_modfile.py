#!/usr/bin/env python3

import sys


msgfmods = {
  'tmt10plex': ['229.162932,*,fix,N-term,TMT6plex', '229.162932,K,fix,any,TMT6plex'],
  'itraq8plex': ['304.205360,*,fix,N-term,iTRAQ8plex', '304.205360,K,fix,any,iTRAQ8plex'],
  'itraq4plex': ['144.102063,*,fix,N-term,iTRAQ4plex', '144.102063,K,fix,any,iTRAQ4plex'],
  'acetyl': ['42.010565,K,opt,any,Acetyl'],
  'phospho': ['79.966331,STY,opt,any,Phospho'],
  'methyl': ['14.015650,K,opt,any,Methyl'],
  'dimethyl': ['28.031300,K,opt,any,Dimethyl'],
  'trimethyl': ['42.046950,K,opt,any,Trimethyl'],
  'carbamidomethyl': ['57.021464,C,fix,any,Carbamidomethyl'],
  'oxidation': ['15.994915,M,opt,any,Oxidation'],
  }

def modpos(mod):
    """Return aminoacid__position of modification
    Because some fixed mods can be on aa * (e.g. N-term)"""
    return '{}__{}'.format(mod[1], mod[3])


def categorize_mod(mod, fixedmods, varmods):
    moddata = [x.split(',') for x in mod]
    for modline in moddata:
        if modline[2] == 'fix':
            try:
                if not modline in fixedmods[modpos(modline)]:
                    fixedmods[modpos(modline)].append(modline)
            except KeyError:
                fixedmods[modpos(modline)] = [modline]
        elif modline[2] == 'opt':
            varmods.append(modline)


nummods = sys.argv[1]
predefined_mods = sys.argv[2]
custom_mods = [] if len(sys.argv) < 4 else sys.argv[3].split(';')

fixedmods = {}
varmods = []

# Parse modifications passed
for modname in predefined_mods.split(';'):
    try:
        mod = msgfmods[modname.lower()]
    except KeyError:
        sys.stderr.write('Could not identify modification "{}", use one of [{}]\n'.format(modname, ', '.join(msgfmods.keys())))
        sys.exit(1)
    categorize_mod(mod, fixedmods, varmods)
for mod in custom_mods:
    categorize_mod(mod, fixedmods, varmods)


# Adjust variable mod weight if competing with fix mods on AA
# TODO make this selectable, since there are mods which CAN co-occur (e.g. TMT on the mod's amine group)
for modline in varmods:
    if modpos(modline) in fixedmods:
        fixmass = sum([float(x[0]) for x in fixedmods[modpos(modline)]])
        modline[0] = str(round(-(fixmass - float(modline[0])), 5))

with open('mods.txt', 'w') as fp:
    fp.write('NumMods={}'.format(nummods))
    for modline in [*[x for f in fixedmods.values() for x in f], *varmods]:
        fp.write('\n{}'.format(','.join(modline)))

#Mass diff with TMT
#Acetyl: ['-187.152367,K,opt,any,Acetyl'],
#Methyl: ['-215.147282,K,opt,any,Methyl'],
#Dimethyl: ['-201.131632,K,opt,any,Dimethyl'],
#Trimethyl: ['-187.115982,K,opt,any,Trimethyl'],
