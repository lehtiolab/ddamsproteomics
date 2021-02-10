#!/usr/bin/env python

import sys
from Bio import SeqIO


def check_fasta(db):
    """Checks if there are no duplicates in the FASTA passed,
    to prevent crashes late in the pipeline (e.g. MSGF+ is fine with
    duplicates, but PSM table isnt)"""
    try:
        SeqIO.index(db, 'fasta')
    except ValueError as e:
        sys.stderr.write('Cannot parse FASTA: {}\n'.format(e))
        sys.exit(130)


if __name__ == '__main__':
    check_fasta(sys.argv[1])
