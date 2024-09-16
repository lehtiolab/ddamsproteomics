#!/usr/bin/env bash

set -eu

echo TMT 16/18 mix, multi DB run
name=tmt16_18mix
baseresults=test_output/${name}
cat "${testdir}/tmt16_mzmls.txt" <(sed 's/0set-A/20set-A/' "${testdir}/tmt18_mzmls.txt" | tail -n+2) | envsubst > test_output/mzmldef
nextflow run -resume -profile test ${repodir}/main.nf --name ${name} --outdir ${baseresults} \
    --input test_output/mzmldef \
    --sampletable "${testdir}/tmt18_setAB_samples.txt" \
    --isobaric '0set-A:tmt16plex:126:131N 20set-A:tmt18plex:131' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --psmconflvl 0.4 --pepconflvl 0.4
