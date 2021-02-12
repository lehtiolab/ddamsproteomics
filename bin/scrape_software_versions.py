#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
import os

regexes = {
    'lehtiolab/ddamsproteomics': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'MSGF+': ['v_msgf.txt', r"v(20[0-9\.]+)"],
    'Dinosaur': ['v_dino.txt', r"([0-9\.]+[0-9])"],
    'Hardklor': ['v_hk.txt', r"([0-9\.]+)"],
    'Kronik': ['v_kr.txt', r"([0-9\.]+)"],
    'Luciphor2': ['v_luci.txt', r"Version\: (.+)"],
    'Percolator': ['v_perco.txt', r"([0-9\.]+)"],
    'msstitch': ['v_mss.txt', r"(\S+)"],
    'OpenMS': ['v_openms.txt', r"Version: ([0-9A-Z\-\.]+)"],
    'DEqMS': ['v_deqms.txt', r"\[1\]..([0-9.]+)"],
}

results = OrderedDict()
results['lehtiolab/ddamsproteomics'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    if os.path.exists(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'lehtiolab/ddamsproteomics-software-versions'
section_name: 'lehtiolab/ddamsproteomics Software Versions'
section_href: 'https://github.com/lehtiolab/ddamsproteomics'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
