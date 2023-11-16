#!/usr/bin/env bash

set -eu

echo TMT16 test denom deqms hardklor keepnapsmsquant
# no hirief, so it can be added in the addsetB test even if not used here
# Test TMT16, DEqMS w denominator, keepnapsmsquant, implicit normalizing (deqms forces normalize)
# Warning: not enough q-values/linear model q-values for gene FDR -> using svm
name=tmt16denomdeq
baseresults=test_output/${name}
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir ${baseresults} \
    --input <(cat "${testdir}/tmt16_mzmls.txt" | envsubst) \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --deqms --keepnapsmsquant --genes \


echo TMT16 phospho + acetyl and total protnorm
name=tmt16_acetyl_phos
baseresults=test_output/${name}
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir ${baseresults} \
    --mzmldef <(grep fr07 "${testdir}/tmt16_mzmls.txt" | envsubst) \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho --ptms Acetyl \
    --proteinconflvl 0.03 --psmconflvl 0.5 --pepconflvl 0.5 \
    --totalproteomepsms "${baseresults}/target_psmtable.txt"

echo TMT16 add a set, carbamyl, set-coloring PCA
# Tests set adding, not using DEqMS, using own defined MSGF mod (carbamyl)
# change name, this mzML is already in  the existing data and spectraIDs will collide
ln -fs "$(pwd)/test-data/ddamsproteomics/tmt16_fr07>_1000@.mzML" "$(pwd)/test-data/ddamsproteomics/linked_tmt16_fr07_1000.mzML"
name=tmt16_addsetB
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir test_output/${name} \
    --mzmldef <(cat "${testdir}/tmt16_setB_mzmls.txt"  | envsubst) \
    --sampletable "${testdir}/tmt16_setAB_samples.txt" \
    --isobaric '0set-A:tmtpro:126:131N set-B:tmt16plex:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation;43.005814,*,opt,N-term,Carbamyl' \
    --genes \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --targetpsms "${baseresults}/target_psmtable.txt" \
    --decoypsms "${baseresults}/decoy_psmtable.txt" \
    --targetpsmlookup "${baseresults}/target_psmlookup.sql" \
    --decoypsmlookup "${baseresults}/decoy_psmlookup.sql" \
    --oldmzmldef ${testdir}/tmt16_mzmls.txt \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt

echo TMT16 frac and non frac mix
name=tmt16_mixfrac_lg
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir test_output/${name} \
    --mzmldef <(cat <(grep fr07 "${testdir}/tmt16_mzmls.txt") <(grep fr08 "${testdir}/tmt16_mzmls.txt" | cut -f1-3) | envsubst) \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --deqms --keepnapsmsquant --genes \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt


# Test single file so we dont get intro trouble with the treat-single-list-as-string NF behaviour
echo TMT16 single-file for escaping in listified steps
name=tmt16_singlefile
nextflow run -resume -profile docker ${repodir}/main.nf --name ${name} --outdir test_output/${name} \
    --mzmldef <(grep fr08 "${testdir}/tmt16_mzmls.txt" | envsubst) \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --genes \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt
