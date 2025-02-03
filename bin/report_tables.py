#!/usr/bin/env python3

import re
import os
import sys
import base64
import shutil
import argparse
from glob import glob
from collections import defaultdict
from datetime import datetime

from lxml import html
from jinja2 import Environment, FileSystemLoader

# Parse plots HTML file

date = datetime.strftime(datetime.now(), '%Y%m%d, %H:%M')

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--version', dest='version', help='WF version')
parser.add_argument('--doi', dest='doi', help='WF DOI') 
parser.add_argument('--templatedir', dest='templatedir', help='dir with template report')
parser.add_argument('--plates', dest='plates', nargs='+', default=['noplates'], help='platenames')
args = parser.parse_args(sys.argv[1:])

# Get template
templatefn = os.path.join(args.templatedir, 'report.html')
with open(os.path.join(args.templatedir, 'bulma.js')) as fp:
    bulma = fp.read()
with open(templatefn) as fp:
    template = Environment(loader=FileSystemLoader(args.templatedir)).from_string(fp.read())

# Plots from plotly
def get_plotly_html(fn):
    with open(fn) as fp:
        plothtml = html.parse(fp)
    plotbox = plothtml.find("body/div[@id='htmlwidget_container']")
    plotcode = plothtml.find("body/script[@type='application/json']")
    plot = html.tostring(plotbox, encoding=str).strip() + html.tostring(plotcode, encoding=str).strip().replace('\\n', '\\\\n').replace('\\\\\\n', '\\\\\\\\n')
    return plot

psmplots = {
        'amount_psms': 'amount_psms.html',
        'miscleav': 'missed_cleavages.html',
        'isomissvals': 'iso_missing_vals.html',
        }
plateplotnames = ['fryield', 'score', 'pif', 'retentiontime', 'precerror', 'fwhm', 'ioninjtime']
plots = defaultdict(dict)

plateplots = defaultdict(defaultdict)
pdir = 'psmplots'
for plotname, fileplot in psmplots.items():
    pfile = os.path.join(pdir, fileplot)
    if os.path.exists(pfile):
        plots[plotname] = get_plotly_html(pfile)
    else:
        plots[plotname] = False
for plotname in plateplotnames:
    for plate in args.plates:
        pfile = os.path.join(pdir, f'PLATE___{plate}___{plotname}.html')
        if os.path.exists(pfile):
            plateplots[plotname][plate] = get_plotly_html(pfile)
        else:
            plateplots[plotname][plate] = False
    if all(x is False for x in plateplots[plotname].values()):
        plateplots[plotname] = False

featnames = [('peptides', 'Peptides'), 
       ('proteins', 'Proteins'), 
       ('genes', 'Gene names'), 
        ('ensg', 'ENSGs'),
]
featplotnames = [('nrfeats', 'Identifications'),
          ('isobaric', 'Isobaric intensities'),
          ('normfactors', 'Isobaric normalization factors'),
          ('precursorarea', 'Precursor area intensity'),
          ('nrpsms', '# PSMs used for isobaric quantitation per identification'),
          ('nrpsmsoverlapping', '# PSMs used for isobaric quantitation per identification for only complete overlapping set'),
          ('percentage_onepsm', 'Percentage of identifications with >1 quantifying PSM in the complete overlapping set'),
          ('ms1nrpeps', '# peptides with MS1 quant per protein (top 3 used)'),
          ('coverage', 'Overall protein coverage'),
          ]

featplotfns = {
        'coverage': ('coverage.html', False),
        'ms1nrpeps': ('ms1nrpeps.html', False),
        'precursorarea': ('precursorarea.html', False),
        'nrfeats': ('nrfeats.html', 'nrfeats__text.html'),
        'nrpsms': ('iso_nrpsms.html', False),
        'nrpsmsoverlapping': ('nrpsmsoverlapping.html', False),
        'percentage_onepsm': ('percentage_onepsm.html', False),
        'isobaric': ('isobaric.html', 'isobaric__text.html'),
        'normfactors': ('normfactors.html', False),
        }
featplots = defaultdict(defaultdict)
for plotname, (pfn, textfn) in featplotfns.items():
    for featname, feattitle in featnames:
        pdir = f'{featname}__plothtml'
        pfile = os.path.join(pdir, pfn)
        if os.path.exists(pfile):
            featplots[plotname][featname] = get_plotly_html(pfile)
            if textfn:
                with open(os.path.join(pdir, textfn)) as fp:
                    featplots[plotname][f'{featname}__text'] = fp.read().strip().split('\n')
        else:
            featplots[plotname][featname] = False
    if all(x is False for x in featplots[plotname].values()):
        featplots[plotname] = False

expplotnames = [
          ('pcagroup', 'Principal component analysis'),
          ('screegroup', 'PCA Scree plot'),
          ]
pcaplotnames = [
          ('pcaset', 'Principal component analysis'),
          ('screeset', 'PCA Scree plot'),
          ]
expplotfns = [('pcaset', 'pca_set.html'),
        ('pcagroup', 'pca_group.html'),
        ('screeset', 'scree_set.html'),
        ('screegroup', 'scree_group.html'),
        ]
expplots = defaultdict(defaultdict)
for plotname, pfn in expplotfns:
    for featname, feattitle in featnames:
        pdir = f'{featname}__plothtml'
        pfile = os.path.join(pdir, pfn)
        if os.path.exists(pfile):
            expplots[plotname][featname] = get_plotly_html(pfile)
        else:
            expplots[plotname][featname] = False
    if all(x is False for x in expplots[plotname].values()):
        expplots[plotname] = False

deqmsplots = defaultdict(dict)
deqmscomps = set()
for featname, feattitle in featnames:
    pdir = f'{featname}__plothtml'
    for fn in glob(f'{pdir}/deqms_volcano_*.png'):
        comp = re.sub(f'^{pdir}/*deqms_volcano_', '', fn)
        comp = re.sub('.png$', '', comp)
        deqmscomps.add(comp)
        with open(fn, 'rb') as fp:
            deqmsplots[featname][comp] = base64.b64encode(fp.read()).decode('utf-8')


#PTMs
ptmplots = defaultdict(dict)
ptmtables = {'summary': [], 'featcount': [], 'overlap': []}
ptmsumtitles = {
        'nr_sets': 'IDed in # overlapping sets', 
        'bioset': 'Experiment set',
        'ptm_residue': 'PTM site',
        'specid': '# sites PSMs',
        'peptide': '# sites peptides',
        'protein': '# sites master proteins',
        'ptmpsmcount': '# PSMs with PTM',
        'ptmpepcount': '# peptides with PTM',
        'ptmprotcount': '# proteins with PTM',
        'protlvl': 'Proteins',
        'peplvl': 'Peptides',
        }
ptm_headers = {
        'summary': [ptmsumtitles[x]
            for x in ['bioset', 'ptmpsmcount', 'ptmpepcount', 'ptmprotcount']],
        'featcount': [ptmsumtitles[x]
            for x in ['bioset', 'ptm_residue', 'specid', 'peptide', 'protein']],
        'overlap': [ptmsumtitles[x]
            for x in ['ptm_residue', 'nr_sets', 'peplvl', 'protlvl']],
        }

pdir = 'ptmplots'
if os.path.exists(pdir) and os.path.isdir(pdir):
    # PTMs in analysis
    for featname, title in [('psm', 'PSMs'), ('peptide', 'Peptides'), ('protein', 'Proteins')]:
        pfile = os.path.join(pdir, f'{featname}__ptms.html')
        if os.path.exists(pfile):
            ptmplots['feats'][title] = get_plotly_html(pfile)
        pfile = os.path.join(pdir, f'{featname}__residues.html')
        if os.path.exists(pfile):
            ptmplots['residues'][title] = get_plotly_html(pfile)
    # PTM tables
    with open('ptm__featcount_table.txt') as fp:
        header = next(fp).strip('\n').split('\t')
        for line in fp:
            line = line.strip('\n').split('\t')
            fields = {ptmsumtitles[f]: line[ix] for ix, f in enumerate(header)}
            ptmtables['summary'].append(fields)
    with open('ptm__table.txt') as fp:
        header = next(fp).strip('\n').split('\t')
        for line in fp:
            line = line.strip('\n').split('\t')
            fields = {ptmsumtitles[f]: line[ix] for ix, f in enumerate(header)}
            ptmtables['featcount'].append(fields)
    if os.path.exists('ptm__overlap.txt'):
        with open('ptm__overlap.txt') as fp:
            header = next(fp).strip('\n').split('\t')
            for line in fp:
                line = line.strip('\n').split('\t')
                fields = {ptmsumtitles[f]: line[ix] for ix, f in enumerate(header)}
                ptmtables['overlap'].append(fields)




# NB Only taking libs from all peptide plot, same libs in allfr anyway
for dirp, dirnames, fns in os.walk('.', followlinks=True):
    for fn in fns:
        if fn.endswith('.min.js'):
            src = os.path.join(dirp, fn)
            dst = os.path.join(dirp, fn.replace('.min.js', '.js'))
            shutil.copy(src, dst)
libs = []
for dirp, dirnames, fns in os.walk('.', followlinks=True):
    for fn in fns:
        srcfn = os.path.join(dirp, fn)
        if not os.path.exists(srcfn) or fn.endswith('.min.js') or fn.endswith('.scss'):
            continue
        with open(srcfn) as fp:
            if fn.endswith('.js'):
                libs.append([fn, f'<script type="text/javascript">{fp.read()}</script>'])
            elif fn.endswith('.css') and not fn.endswith('.scss'):
                libs.append([fn, f'<style type="text/css">{fp.read()}</style>'])
if len(libs):
    with open('libs.js', 'w') as fp:
        for fn in ['htmlwidgets.js', 'plotly.js', 'typedarray.js', 'jquery.js', 'crosstalk.min.css',
                'crosstalk.js', 'plotly-htmlwidgets.css', 'plotly-latest.js']:
            lib = [x[1] for x in libs if x[0] == fn][0]
            fp.write(f'{lib}\n')


# Summary table
tabletitles = {
        'Set': 'Experiment set', 
        'no_pep_proteins': 'Peptides/protein (unique, median)',
        'no_pep_genes': 'Peptides/protein (genecentric, unique, median)',
        'no_psm_proteins': 'PSMs/protein for quant (median)',
        'no_psm_genes': 'PSMs/protein for quant (genecentric, median)',
        'nr_proteins': 'Proteins ID (1%FDR)',
        'nr_genes': 'Proteins ID (gene centric, 1%FDR)',
        'nr_assoc': 'Proteins ID (symbol centric, 1%FDR)',
        'nr_proteins_q': 'Proteins ID and quant. (1%FDR)',
        'nr_genes_q': 'Proteins ID and quant. (gene centric, 1%FDR)',
        'nr_assoc_q': 'Proteins ID and quant. (symbol centric, 1%FDR)',
        'Non-shared (unique)': 'Peptides (unique, 1%FDR)',
        'psmcount': 'PSMs (total)',
        }
summary_field_order = ['nr_proteins_q', 'nr_proteins', 'nr_genes_q', 'nr_genes', 'nr_assoc', 'nr_assoc_q', 'Non-shared (unique)', 'no_pep_proteins', 'no_pep_genes', 'psmcount', 'no_psm_proteins', 'no_psm_genes',
        ]

summary_table = defaultdict(dict)
for fn in glob('*__summary.txt'):
    with open(fn) as fp:
        header = next(fp).strip('\n').split('\t')
        for line in fp:
            line = line.strip('\n').split('\t')
            fields = {f: line[ix] for ix, f in enumerate(header)}
            setname = fields.pop('Set')
            for f in fields:
                summary_table[setname][f] = fields[f]
for _s, fields in summary_table.items():
    summary_fields = [x for x in summary_field_order if x in fields]
    break

# PSM tables
psmtables = {'ids': [], 'miscleav': []}
with open('psmids') as fp:
    head = next(fp).strip().split()
    plates = defaultdict(defaultdict)
    for line in fp:
        lnmap = {head[ix]: x for ix, x in enumerate(line.strip().split('\t'))}
        if lnmap['name'] == 'MS2 scans':
            plates[lnmap['plateID']]['scans'] = lnmap['count']
        elif lnmap['name'] == 'PSMs IDed':
            plates[lnmap['plateID']].update({'psms': lnmap['count'], 'pc': lnmap['labeltext']})
psmtables['ids'] = [[p, nms['scans'], nms['psms'], nms['pc']] for p, nms in plates.items()]

with open('miscleav') as fp:
    head = next(fp).strip().split()
    plates = defaultdict()
    for line in fp:
        lnmap = {head[ix]: x for ix, x in enumerate(line.strip().split('\t'))}
        print(lnmap)
        psmtables['miscleav'].append([lnmap['plateID'], lnmap['missed_cleavage'], lnmap['nrpsms'], lnmap['IDed']])


# Overlap
overlap = defaultdict(dict)
for feattype, _ft in featnames:
    fn = f'{feattype}__overlap'
    if os.path.exists(fn):
        with open(fn) as fp:
            _h = next(fp)
            for line in fp:
                overlapnrsets, nr_feats = line.strip().split('\t')
                overlap[feattype][overlapnrsets] = nr_feats
    else:
        overlap[feattype] = False

# Isobaric normalization factors
normfactable = defaultdict(list)
for feattype, _ft in featnames:
    pdir = f'{feattype}__plothtml'
    fn = os.path.join(pdir, 'allnormfacs')
    if os.path.exists(fn):
        with open(fn) as fp:
            for line in fp:
                normfactable[feattype].append(line.strip().split('\t'))


# Warning box
warnings = []
for fn in glob('warnings*'):
    if os.path.exists(fn):
        with open(fn) as fp:
            for line in fp:
                if warn := line.strip():
                    warnings.append(warn)


# Write to template
with open('report_groovy_template.html', 'w') as fp:
    fp.write(template.render(reportdate=date,
        version=args.version,
        doi=args.doi,
        plots=plots,
        plateplots=plateplots,
        featplots=featplots,
        featnames=featnames,
        featplotnames=featplotnames,
        plates=args.plates,
        expplots=expplots,
        expplotnames=expplotnames,
        pcaplotnames=pcaplotnames,
        deqmsplots=deqmsplots,
        deqmscomps=deqmscomps,
        tabletitles=tabletitles,
        psmtables=psmtables,
        summary_fields=summary_fields,
        summary_table=summary_table,
        overlap=overlap,
        ptmplots=ptmplots,
        ptmtables=ptmtables,
        ptmtitles=ptm_headers,
        normfacs=normfactable,
        warnings=warnings,
        ))
