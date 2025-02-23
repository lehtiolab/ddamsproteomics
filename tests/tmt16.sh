#!/usr/bin/env bash

set -eu
  
echo TMT16 test denom deqms hardklor keepnapsmsquant
# no hirief, so it can be added in the addsetB test even if not used here
# Test TMT16, DEqMS w denominator, keepnapsmsquant, implicit normalizing (deqms forces normalize)
# Warning: no decoys for any set (aka only setA)
name=tmt16denomdeq
resultsdir=test_output/${name}
mkdir -p $resultsdir
cat "${testdir}/tmt16_mzmls.txt" | envsubst > ${resultsdir}/mzmldef
$NXFCMD --name ${name} --outdir ${resultsdir} \
    --input ${resultsdir}/mzmldef \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation;43.005814,*,opt,N-term,Unknown' \
    --maxmiscleav 2 \
    --pepconflvl 0.05 \
    --deqms --keepnapsmsquant --genes


echo TMT16 phospho + acetyl and total protnorm
name=tmt16_acetyl_phos
baseresults=${resultsdir}
resultsdir=test_output/${name}
mkdir -p $resultsdir
cat <(head -n1 "${testdir}/tmt16_mzmls.txt") <(grep fr07 ${testdir}/tmt16_mzmls.txt) | envsubst > ${resultsdir}/mzmldef
$NXFCMD --name ${name} --outdir test_output/${name} \
    --input ${resultsdir}/mzmldef \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho --ptms Acetyl \
    --proteinconflvl 0.03 --psmconflvl 0.5 --pepconflvl 0.5 \
    --genes \
    --mediannormalize --totalproteomepsms "${baseresults}/target_psmtable.txt"

### FIXME add anotehr totalproteome without median norm

echo TMT16 add a set, carbamyl, set-coloring PCA
# Tests set adding, not using DEqMS, using own defined MSGF mod (carbamyl)
# change name, this mzML is already in  the existing data and spectraIDs will collide
ln -fs "${testdata}/tmt16_fr07>_1000@.mzML" "${testdata}/linked_tmt16_fr07_1000.mzML"
name=tmt16_addsetB
baseresults=${resultsdir}
resultsdir=test_output/${name}
mkdir -p $resultsdir
cat "${testdir}/tmt16_mzmls.txt" | envsubst > ${resultsdir}/oldmzmls
cat "${testdir}/tmt16_setB_mzmls.txt" | envsubst > ${resultsdir}/mzmldef
$NXFCMD --name ${name} --outdir test_output/${name} \
    --input ${resultsdir}/mzmldef \
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
    --oldmzmldef test_output/${name}/oldmzmls \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt

echo TMT16 frac and non frac mix
name=tmt16_mixfrac_lg
baseresults=${resultsdir}
resultsdir=test_output/${name}
mkdir -p $resultsdir
cat <(head -n1 "${testdir}/tmt16_mzmls.txt") <(grep fr07 ${testdir}/tmt16_mzmls.txt) <(grep fr08 ${testdir}/tmt16_mzmls.txt | cut -f1-3) | envsubst > ${resultsdir}/mzmldef
$NXFCMD --name ${name} --outdir test_output/${name} \
    --input ${resultsdir}/mzmldef \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --deqms --keepnapsmsquant --genes \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt


# Test single file so we dont get intro trouble with the treat-single-list-as-string NF behaviour
echo TMT16 single-file for escaping in listified steps
name=tmt16_singlefile
resultsdir=test_output/${name}
mkdir -p $resultsdir
cat <(head -n1 ${testdir}/tmt16_mzmls.txt) <(grep fr08 ${testdir}/tmt16_mzmls.txt) | envsubst > ${resultsdir}/mzmldef
$NXFCMD --name ${name} --outdir test_output/${name} \
    --input ${resultsdir}/mzmldef \
    --sampletable "${testdir}/tmt16_samples.txt" \
    --hardklor --isobaric '0set-A:tmtpro:126:131N' \
    --tdb "${testdata}/tmt16_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --genes \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --hirief https://github.com/nf-core/test-datasets/raw/6defbf8a92a46b0ac48bb05f9ad96b62716b4a5d/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt
