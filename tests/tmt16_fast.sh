#!/usr/bin/env bash

set -eu
  
echo TMT16 test with hiRIEF, isobaric and Phospho
# no hirief, so it can be added in the addsetB test even if not used here
# Test TMT16, DEqMS w denominator, keepnapsmsquant, implicit normalizing (deqms forces normalize)
# Warning: no decoys for any set (aka only setA)
name=tmt16denomdeq
resultsdir=test_output/${name}
mkdir -p $resultsdir
cat "${testdir}/tmt16_mzmls.txt" | envsubst > ${resultsdir}/mzmldef
nextflow run -resume -profile test ${repodir}/main.nf --name ${name} --outdir ${resultsdir} \
    --input ${resultsdir}/mzmldef \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --deqms --keepnapsmsquant --genes \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt
    # FIXME cannot run with carbamyl +43 -> -261 and Phospho, luciprep crash \
