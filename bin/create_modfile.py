#!/usr/bin/env python3

import sys
import json
from mods import Mods


def main():
    nummods = sys.argv[1]
    modlibfile = sys.argv[2]
    search_engine = sys.argv[3]
    passed_mods = sys.argv[4].split(';')
    fixedmods = {}
    varmods = []
    # Parse modifications passed
    parsemods = Mods()
    parsemods.parse_msgf_modfile(modlibfile, passed_mods)
    modfn = 'mods.txt'
    if search_engine == 'msgf':
        with open(modfn, 'w') as fp:
            fp.write('NumMods={}'.format(nummods))
            for modline in parsemods.get_msgf_modlines():
                fp.write(f'\n{modline}')
    elif search_engine == 'sage':
        with open(modfn, 'w') as fp:
            json.dump(parsemods.get_sage_mods(), fp)


if __name__ == '__main__':
    main()
