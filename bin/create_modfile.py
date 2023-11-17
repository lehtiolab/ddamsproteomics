#!/usr/bin/env python3

import sys
from mods import Mods


def main():
    nummods = sys.argv[1]
    modlibfile = sys.argv[2]
    passed_mods = sys.argv[3].split(';')
    fixedmods = {}
    varmods = []
    # Parse modifications passed
    msgfmods = Mods()
    msgfmods.parse_msgf_modfile(modlibfile, passed_mods)
    with open('mods.txt', 'w') as fp:
        fp.write('NumMods={}'.format(nummods))
        for modline in msgfmods.get_msgf_modlines():
            fp.write(f'\n{modline}')


if __name__ == '__main__':
    main()
