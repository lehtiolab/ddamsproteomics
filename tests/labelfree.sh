#!/usr/bin/env bash

set -eu

echo Normal labelfree test
name=lf
cat ${testdir}/lf_mzmls.txt | envsubst > test_output/mzmldef
nextflow run -resume -profile test ${repodir}/main.nf --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --genes \
    --tdb ${testdata}/lf.fa \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --mods 'carbamidomethyl;oxidation'

# Best PSMs
#setA T: 0.007246 / D: 0.00129
#setB T: 0.002949853 / D: 0.005882
# Keep high pep conf lvl
# Cannot only do no-decoy warning since decoys have higher q-value than targets
# No target in setA / no decoy setB: 0.005

echo Labelfree run with warnings: no decoy in setB, no target in setA
name=lf_notarget
cat ${testdir}/lf_mzmls.txt | envsubst > test_output/mzmldef
nextflow run -resume -profile test ${repodir}/main.nf --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --genes \
    --tdb ${testdata}/lf.fa \
    --psmconflvl 0.005 --pepconflvl 0.2 \
    --mods 'carbamidomethyl;oxidation'

echo Labelfree run without fractions, with warnings: no decoy in setB, no target in setA
name=lf_nofrac_notarget
cat ${testdir}/lf_mzmls_nofrac.txt | envsubst > test_output/mzmldef
nextflow run -resume -profile test ${repodir}/main.nf --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --genes \
    --tdb ${testdata}/lf.fa \
    --psmconflvl 0.005 --pepconflvl 0.2 \
    --mods 'carbamidomethyl;oxidation'


echo Single file labelfree test
name=lf_singlefile
cat <(head -n1 ${testdir}/lf_mzmls.txt) <(grep setA ${testdir}/lf_mzmls.txt) | envsubst > test_output/mzmldef
nextflow run -resume -profile test ${repodir}/main.nf --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --genes \
    --tdb ${testdata}/lf.fa \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --mods 'carbamidomethyl;oxidation'
