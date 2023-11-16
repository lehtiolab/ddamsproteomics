#!/usr/bin/env bash

set -eu

echo TMT 16/18 mix, multi DB run
name=tmt16_18mix
baseresults=test_output/${name}
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir ${baseresults} \
    --mzmldef <(cat "${testdir}/tmt16_mzmls.txt" <(sed 's/0set-A/0setB/' "${testdir}/tmt18_mzmls.txt") | envsubst) \
    --sampletable "${testdir}/tmt18_setAB_samples.txt" \
    --isobaric '0set-A:tmt16plex:126:131N 0setB:tmt18plex:131' \
    --tdb "${testdata}/tmt16*_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --psmconflvl 0.4 --pepconflvl 0.4
