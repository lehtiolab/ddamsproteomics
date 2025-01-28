import sys

NON_BLOCKING_MODS = {
        'GG': ['TMTpro', 'TMT6plex', 'iTRAQ8plex', 'iTRAQ4plex'],
        }

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


class Mods:
    '''
    name is UNIMOD so tmt10plex is called TMT6plex
    name_lower is how it can be referred to by e.g. other programs e.g. OpenMS, tmt10plex
    TODO needs more attention
    '''

    def __init__(self):
        self.mods = []
        self.fixedmods = []
        self.varmods = []
        self.bymass = {}
        self.has_varmods_on_fixmod_residues = False

    def parse_msgf_modfile(self, modfile, mods_passed):
        # FIXME make sure parsing is only Unimod/mass, then set
        # fixed/var, pos, res yourself in this method
        mods_to_find = [x.lower() for x in mods_passed]
        with open(modfile) as fp:
            for line in fp:
                line = line.strip('\n')
                if line == '' or line[0] == '#' or 'NumMods' in line:
                    continue
                # TODO validate line
                msplit = line.split(',')
                name = msplit[4]
                pos = msplit[3]
                varfix = msplit[2]
                residues = set(msplit[1])
                # tmt6plex can be hidden tmt10plex, same UNIMOD mass/name
                if name.lower() == 'tmt6plex' and 'tmt10plex' in mods_to_find:
                    lowername = 'tmt10plex'
                elif name.lower() not in mods_to_find:
                    continue
                else:
                    lowername = name.lower()
                for res in residues:
                    self.mods.append({
                            'name': name, 'mass': float(msplit[0]),
                            'adjusted_mass': False,
                            'residue': res, 'var': varfix == 'opt',
                            'pos': pos, 'name_lower': lowername
                            })
        # See if user has defined own mods also
        for mtofind in mods_passed:
            moddef = mtofind.split(',')
            if len(moddef) == 5:
                # Found a custom mod defined by user
                name = moddef[4]
                pos = moddef[3]
                varfix = moddef[2]
                residues = set(moddef[1])
                for res in residues:
                    self.mods.append({
                        'name': name, 'mass': float(moddef[0]),
                        'adjusted_mass': False,
                        'residue': res, 'var': varfix == 'opt',
                        'pos': pos, 'name_lower': name.lower()
                        })

        fixedpos = {}
        for mod in self.mods:
            if mod['var']:
                self.varmods.append(mod)
            else:
                self.fixedmods.append(mod)
                res = mod['residue']
                if res in fixedpos:
                    fixedpos[res].append(mod)
                else:
                    fixedpos[res] = [mod]

        # get blocking/nonblocking mods and adjust mass (fake mass)
        for mod in self.varmods:
            nonblocked_fixed = NON_BLOCKING_MODS.get(mod['name'], []) 
            if mod['residue'] in fixedpos:
                self.has_varmods_on_fixmod_residues = True
            adjustment = 0
            for fmod in fixedpos.get(mod['residue'], []):
                if fmod['name'] not in nonblocked_fixed:
                    adjustment += fmod['mass']
            mod['adjusted_mass'] = round(-(adjustment - mod['mass']), 5)

    def get_grouped_mods_by_mass(self):
        grouped = {}
        for mod in self.mods:
            mass = self.get_mass_or_adj(mod)
            if mass not in grouped:
                grouped[mass] = [mod]
            else:
                grouped[mass].append(mod)
                check_names = set(x['name'] for x in grouped[mass])
                if len(check_names) != 1:
                    print('Cannot have two modifications of the same mass but different names')
                    sys.exit(1)
        return grouped

    def get_msgf_modlines(self):
        grouped = self.get_grouped_mods_by_mass()
        for mass, mods in grouped.items():
            name = mods[0]['name']
            line_res = {}
            for mod in mods:
                var = int(mod['var']) # T/F -> 1/0
                mid = f'{var}__{mod["pos"]}'
                if mid not in line_res:
                    line_res[mid] = [mod['residue']]
                else:
                    line_res[mid].append(mod['residue'])
            for mid, residues in line_res.items():
                var, pos = mid.split('__')
                vf = 'opt' if int(var) else 'fix'
                yield f'{mass},{"".join(residues)},{vf},{pos},{name}'

    def get_sage_mods(self):
        grouped = self.get_grouped_mods_by_mass()
        outmods = {'static_mods': {}, 'variable_mods': {}}
        posmap = {'N-term': '^',
                'C-term': '$',
                }
        mtypemap = {False: 'static_mods', True: 'variable_mods'}
        for mass, mods in grouped.items():
            line_res = {}
            for mod in mods:
                mtype = mtypemap[mod['var']]
                pos = '' if mod['pos'] == 'any' else posmap[mod['pos']]
                res = '' if mod['residue'] == '*' else mod['residue']
                sage_residue = f'{pos}{res}'
                try:
                    outmods[mtype][sage_residue].append(mass)
                except KeyError:
                    outmods[mtype][sage_residue] = [mass]
        outmods['static_mods'] = {k: sum(v) for k,v in outmods['static_mods'].items()}
        return outmods

    def get_mass_or_adj(self, mod):
        return mod['adjusted_mass'] or mod['mass']
 
    def get_luci_input_mod_line(self, mod):
        '''Doing this for each residue since the adjusted masses can differ
        per residue (due to competition)'''
        if mod['pos'] == 'N-term':
            residue = '['
        elif mod['pos'] == 'C-term':
            residue = ']'
        else:
            residue = mod['residue']
        return f'{residue} {self.get_mass_or_adj(mod)}'

    def msgfmass_mod_dict(self):
        '''Create MSGF output mass (round(x,3) ) to mod lookup'''
        mod_map = {}
        for mod in self.mods:
            try:
                mod_map[round(self.get_mass_or_adj(mod), 3)].append(mod)
            except KeyError:
                mod_map[round(self.get_mass_or_adj(mod), 3)] = [mod]
        return mod_map

    def lucimass_mod_dict(self):
        '''Create luciphor output mass (int(x+ aa) ) to mod lookup'''
        modmap = {}
        for mod in self.varmods:
            # round (, None) generates an int for luciphor
            mass = round(aa_weights_monoiso[mod['residue']] + self.get_mass_or_adj(mod))
            modmap[f'{mod["residue"]}{mass}'] = mod
        return modmap


