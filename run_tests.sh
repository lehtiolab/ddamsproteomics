#!/usr/bin/env bash

red=$(tput setaf 1)
green=$(tput setaf 2)
reset=$(tput sgr0)


export NXF_VER=24.04.4

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
mkdir test_output

declare -A results

test_names=(
    labelfree
    labelfree_phos
    tmt16
    tmt18
    tmt16_18
    tims
)

if [ -z "$1" ]
then
       	for testname in ${test_names[@]}
       	do
	       	bash "${testdir}/${testname}.sh"
	       	results["$testname"]="$?"
       	done
	for testname in ${test_names[@]}
	do
	    if [[ ${results["$testname"]} == 0 ]]
	    then
	        echo ${green}"$testname" - SUCCESS ${reset}
	    else
	        echo ${red}"$testname" - FAIL ${reset}
	    fi
	done
else
       	bash "${testdir}/$1.sh"
       	if [[ "$?" == 0 ]]
       	then
	       	echo ${green}"$1" - SUCCESS ${reset}
       	else
	       	echo ${red}"$1" - FAIL ${reset}
       	fi
fi
