#!/usr/bin/env bash

set -eu

echo TMT18 phos
# Test TMT18, Phos, 
#DEqMS w denominator, implicit normalizing (deqms forces normalize)
# Warning: not enough q-values/linear model q-values for gene FDR -> using svm
name=tmt18phos
baseresults=test_output/${name}
nextflow run -resume -profile standard,docker main.nf --name ${name} --outdir ${baseresults} \
    --mzmldef "${testdir}/tmt18_mzmls.txt" \
    --sampletable "${testdir}/tmt18_samples.txt" \
    --hardklor --isobaric '0set-A:tmt18plex:126:131N' \
    --tdb "${testdata}/tmt18_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho --psmconflvl 0.2 --pepconflvl 0.2 \
    --deqms --genes


echo TMT18 phos add a set
# Test TMT18, Phos, 
#DEqMS w denominator, keepnapsmsquant, implicit normalizing (deqms forces normalize)
# Warning: not enough q-values/linear model q-values for gene FDR -> using svm
name=tmt18phos_addset
ln -fs "$(pwd)/test-data/ddamsproteomics/tmt18_fr06_1000.mzML" "$(pwd)/test-data/ddamsproteomics/linked_tmt18_fr06_1000.mzML"
nextflow run -resume -profile standard,docker main.nf --name ${name} --outdir test_output/${name} \
    --mzmldef <(sed 's/tmt18_fr/linked_tmt18_fr/;s/set-A/setB/' "${testdir}/tmt18_mzmls.txt") \
    --sampletable "${testdir}/tmt18_setAB_samples.txt" \
    --hardklor --isobaric '0set-A:tmt18plex:sweep 0setB:tmt18plex:sweep' \
    --tdb "${testdata}/tmt18_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho --psmconflvl 0.02 --pepconflvl 0.05 \
    --targetpsms "${baseresults}/target_psmtable.txt" \
    --decoypsms "${baseresults}/decoy_psmtable.txt" \
    --targetpsmlookup "${baseresults}/target_psmlookup.sql" \
    --decoypsmlookup "${baseresults}/decoy_psmlookup.sql" \
    --ptmpsms "${baseresults}/ptm_psmtable.txt" \
    --oldmzmldef ${testdir}/tmt18_mzmls.txt \
    --deqms --genes


echo TMT18 rerun with different settings post PSMs
# No need for PSM conf lvl bc it is used in percolator before PSM table
# But pep conf level is used also in QC so needs to be here
name=tmt18phos_rerun
nextflow run -resume -profile standard,docker main.nf --name ${name} --outdir test_output/${name} \
    --sampletable "${testdir}/tmt18_samples.txt" \
    --isobaric '0set-A:tmt18plex:131' \
    --tdb "${testdata}/tmt18_fa.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms Phospho \
    --targetpsms "${baseresults}/target_psmtable.txt" \
    --decoypsms "${baseresults}/decoy_psmtable.txt" \
    --ptmpsms "${baseresults}/ptm_psmtable.txt" \
    --targetpsmlookup "${baseresults}/target_psmlookup.sql" \
    --decoypsmlookup "${baseresults}/decoy_psmlookup.sql" \
    --pepconflvl 0.05 \
    --oldmzmldef ${testdir}/tmt18_mzmls.txt
