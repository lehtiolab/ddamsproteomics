#!/usr/bin/env bash

export NXF_VER=22.10.5

export testdir=tests/
export testdata=static-resources/test-data/ddamsproteomics

rm -r test_output

bash ${testdir}/labelfree.sh
bash ${testdir}/labelfree_phos.sh
bash ${testdir}/tmt16.sh
bash ${testdir}/tmt18.sh
bash ${testdir}/tmt16_18.sh
bash ${testdir}/tims.sh

# FIXME look through all the QCs to check if runs are correct
