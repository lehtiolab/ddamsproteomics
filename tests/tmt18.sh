#!/usr/bin/env bash

set -eu

echo TMT18 phos, no MS1 found somehow
# Test TMT18, Phos, 
#DEqMS w denominator, implicit normalizing (deqms forces normalize)
# Warning: not enough q-values/linear model q-values for gene FDR -> using svm
# Also no MS1 values matching PSMs, which is not specified
# Also no PTMs calculated in luciphor - warning but still plots, just no FLR
name=tmt18phos
baseresults=test_output/${name}
cat "${testdir}/tmt18_mzmls.txt" | envsubst > test_output/mzmldef
$NXFCMD --name ${name} --outdir ${baseresults} \
    --input test_output/mzmldef \
    --sampletable "${testdir}/tmt18_samples.txt" \
    --hardklor --isobaric '0set-A:tmt18plex:126:131N' \
    --tdb "${testdata}/tmt18_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho --psmconflvl 0.2 --pepconflvl 0.2 \
    --deqms --genes


echo TMT18 phos add a set
# Test TMT18, Phos, 
#DEqMS w denominator, keepnapsmsquant, implicit normalizing (deqms forces normalize)
# Not enough PSMs for decoy in new set -> No proteins/genes in QC at all, also not for old set
name=tmt18phos_addset
mkdir -p test_output/${name}
ln -fs "${testdata}/tmt18_fr06_1000.mzML" "${testdata}/linked_tmt18_fr06_1000.mzML"
cat "${testdir}/tmt18_mzmls.txt" | envsubst > test_output/${name}/oldmzmls
sed 's/tmt18_fr/linked_tmt18_fr/;s/0set-A/20set-A/' "${testdir}/tmt18_mzmls.txt" | envsubst > test_output/mzmldef
$NXFCMD --name ${name} --outdir test_output/${name} \
    --input test_output/mzmldef \
    --sampletable "${testdir}/tmt18_setAB_samples.txt" \
    --hardklor --isobaric '0set-A:tmt18plex:sweep 20set-A:tmt18plex:sweep' \
    --tdb "${testdata}/tmt18_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho --psmconflvl 0.02 --pepconflvl 0.05 \
    --targetpsms "${baseresults}/target_psmtable.txt" \
    --decoypsms "${baseresults}/decoy_psmtable.txt" \
    --targetpsmlookup "${baseresults}/target_psmlookup.sql" \
    --decoypsmlookup "${baseresults}/decoy_psmlookup.sql" \
    --ptmpsms "${baseresults}/ptm_psmtable.txt" \
    --oldmzmldef test_output/${name}/oldmzmls \
    --deqms --genes


echo TMT18 rerun with different settings post PSMs, --noms1quant
# No need for PSM conf lvl bc it is used in percolator before PSM table
# But pep conf level is used also in QC so needs to be here
name=tmt18phos_rerun
mkdir -p test_output/${name}
cat "${testdir}/tmt18_mzmls.txt" | envsubst > test_output/${name}/oldmzmls
$NXFCMD --name ${name} --outdir test_output/${name} \
    --sampletable "${testdir}/tmt18_samples.txt" \
    --isobaric '0set-A:tmt18plex:131' \
    --tdb "${testdata}/tmt18_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho \
    --noms1quant \
    --targetpsms "${baseresults}/target_psmtable.txt" \
    --decoypsms "${baseresults}/decoy_psmtable.txt" \
    --ptmpsms "${baseresults}/ptm_psmtable.txt" \
    --targetpsmlookup "${baseresults}/target_psmlookup.sql" \
    --decoypsmlookup "${baseresults}/decoy_psmlookup.sql" \
    --pepconflvl 0.05 \
    --oldmzmldef test_output/${name}/oldmzmls
