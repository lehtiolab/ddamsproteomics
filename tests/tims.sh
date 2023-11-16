#!/usr/bin/env bash

set -eu

echo TIMS TMT16, noms1quant

name=tims_tmt16
baseresults=test_output/${name}
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir ${baseresults} \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --noms1quant \
    --mzmldef <(cat "${testdir}/tims_mzmls.txt" | envsubst) \
    --genes \
    --isobaric '0set-A:tmtpro:126:131C' \
    --tdb "${testdata}/tims_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --psmconflvl 0.2 --pepconflvl 0.2 --proteinconflvl 0.2 \
    --activation cid --keepnapsmsquant
