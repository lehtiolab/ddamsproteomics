#!/usr/bin/env python3

import sys

# variable mods which do not obstruct a list of fixed mods:
NON_BLOCKING_MODS = {
        'GG': ['TMTpro', 'TMT6plex', 'iTRAQ8plex', 'iTRAQ4plex'],
        }

def get_msgfmods(modfile):
    msgfmods = {}
    with open(modfile) as fp:
        for line in fp:
            line = line.strip('\n')
            if line == '' or line[0] == '#':
                continue
            moddata = line.split(',')
            try:
                msgfmods[moddata[-1].lower()].append(line)
            except KeyError:
                msgfmods[moddata[-1].lower()] = [line]
    if 'tmt6plex' in msgfmods:
        msgfmods['tmt10plex'] = msgfmods['tmt6plex']
    return msgfmods


def modpos(mod):
    """Return aminoacid__position of modification
    Because some fixed mods can be on aa * (e.g. N-term)"""
    return '{}__{}'.format(mod[1], mod[3])


def categorize_mod(mods, fixedmods, varmods):
    """Puts mods into varmods/fixedmods lists"""
    moddata = [x.split(',') for x in mods]
    for modline in moddata:
        try:
            checkmass = float(modline[0])
        except ValueError:
            raise RuntimeError('Incorrect mod specification, need weight')
        if modline[2] == 'fix':
            try:
                if not modline in fixedmods[modpos(modline)]:
                    fixedmods[modpos(modline)].append(modline)
            except KeyError:
                fixedmods[modpos(modline)] = [modline]
        elif modline[2] == 'opt':
            if not modline in varmods:
                varmods.append(modline)
        else:
            raise RuntimeError('Incorrect mod specification')


def parse_cmd_mod(cmdmod, msgfmods):
    try:
        return msgfmods[cmdmod.lower()]
    except KeyError:
        # Mod is a line of MSGF mod by user
        return [cmdmod]


def main():
    nummods = sys.argv[1]
    modlibfile = sys.argv[2]
    passed_mods = sys.argv[3]
    fixedmods = {}
    varmods = []
    # Parse modifications passed
    msgfmods = get_msgfmods(modlibfile)
    for passedmod in passed_mods.split(';'):
        mod = parse_cmd_mod(passedmod, msgfmods)
        try:
            categorize_mod(mod, fixedmods, varmods)
        except Exception:
            sys.stderr.write('Could not identify modification "{}", use one of [{}]\n'.format(passedmod, ', '.join(msgfmods.keys())))
            sys.exit(1)
    # Adjust variable mod weight if competing with fix mods on AA
    # TODO make this selectable, since there are mods which CAN co-occur (e.g. TMT on the mod's amine group)
    for modline in varmods:
        mp = modpos(modline)
        mmass, mres, mfm, mprotpos, mname = modline
        if mp in fixedmods:
            nonblocked_fixed = NON_BLOCKING_MODS[mname] if mname in NON_BLOCKING_MODS else []
            blocked_fixmass = sum([float(x[0]) for x in fixedmods[mp] if x[-1] not in nonblocked_fixed])
            modline[0] = str(round(-(blocked_fixmass - float(modline[0])), 5))
    with open('mods.txt', 'w') as fp:
        fp.write('NumMods={}'.format(nummods))
        for modline in [*[x for f in fixedmods.values() for x in f], *varmods]:
            fp.write('\n{}'.format(','.join(modline)))


if __name__ == '__main__':
    main()
