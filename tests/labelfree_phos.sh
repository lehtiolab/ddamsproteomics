#!/usr/bin/env bash
set -eu

echo Phospho, labelfree, one luciphor setB not enough PSMs
name=labelfree_phos
lfphos_dir=test_output/${name}
cat ${testdir}/lf_mzmls.txt | envsubst > test_output/mzmldef
$NXFCMD --name ${name} \
    --outdir ${lfphos_dir} \
    --input test_output/mzmldef \
    --genes \
    --tdb ${testdata}/lf.fa \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --minpsms_luciphor 3 \
    --report_seqmatch "${testdata}/seqmatch_lfphos.fa;${testdata}/seqmatch lfphos@2.fa" \
    --mods 'carbamidomethyl;oxidation' \
    --locptms 'Phospho'

# LF phos where we add replace a set with a "new" set
name=labelfree_phos_addset
cat <(head -n1 ${testdir}/lf_mzmls.txt) <(grep setA ${testdir}/lf_mzmls.txt) | envsubst > test_output/mzmldef
$NXFCMD --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --oldmzmldef ${testdir}/lf_mzmls.txt \
    --genes \
    --tdb ${testdata}/lf.fa \
    --psmconflvl 0.2 --pepconflvl 0.2 \
    --minpsms_luciphor 3 \
    --mods 'carbamidomethyl;oxidation' \
    --targetpsms ${lfphos_dir}/target_psmtable.txt \
    --decoypsms ${lfphos_dir}/decoy_psmtable.txt \
    --targetpsmlookup ${lfphos_dir}/target_psmlookup.sql \
    --decoypsmlookup ${lfphos_dir}/decoy_psmlookup.sql \
    --ptmpsms ${lfphos_dir}/ptm_psmtable.txt \
    --locptms 'Phospho'

# FIXME why are lookup outputs different between these two?

# Test for when there are no decoy PSMs: 
echo LF phos no decoy
name=labelfree_phos_nodecoy
cat <(head -n1 ${testdir}/lf_mzmls.txt) <(grep setA ${testdir}/lf_mzmls.txt) | envsubst > test_output/mzmldef
$NXFCMD --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --genes \
    --tdb ${testdata}/lf.fa \
    --mods 'carbamidomethyl;oxidation' \
    --psmconflvl 0.012 --pepconflvl 0.2 \
    --locptms 'Phospho'

# Test for when there are no target PSMs: crash in createPSMTable
set +eu
echo Run in which no target PSMs are found - crashes
name=labelfree_phos_notarget
cat <(head -n1 ${testdir}/lf_mzmls.txt) <(grep setA ${testdir}/lf_mzmls.txt) | envsubst > test_output/mzmldef
$NXFCMD --name ${name} \
    --outdir test_output/${name} \
    --input test_output/mzmldef \
    --psmconflvl 0.001 \
    --genes \
    --tdb ${testdata}/lf.fa \
    --mods 'carbamidomethyl;oxidation' \
    --locptms 'Phospho'
if [[ "$?" == 1 ]]
then
    exit 0
else
    exit 1
fi
