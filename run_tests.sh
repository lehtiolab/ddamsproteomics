#!/usr/bin/env bash

export NXF_VER=22.10.5

rundir=$(pwd)
export repodir=$(dirname "$(realpath -s "$0")")
export testdir="${repodir}/tests/"
export testdata="${rundir}/static-resources/test-data/ddamsproteomics"

if [ -e "${testdata}" ]
then
    cd "${testdata}" && git checkout ddamsproteomics_test_data && cd "${rundir}"
else
    git clone --single-branch --branch ddamsproteomics_test_data https://github.com/lehtiolab/static-resources
fi

[ -e test_output ] && rm -r test_output

bash ${testdir}/labelfree.sh
bash ${testdir}/labelfree_phos.sh
bash ${testdir}/tmt16.sh
bash ${testdir}/tmt18.sh
bash ${testdir}/tmt16_18.sh
bash ${testdir}/tims.sh

# FIXME look through all the QCs to check if runs are correct
