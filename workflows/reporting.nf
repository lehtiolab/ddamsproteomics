/* Reporting workflow for DDA pipeline */

include { listify; stripchars_infile; get_regex_specialchars } from '../modules.nf' 

process countMS2sPerPlate {
  label 'python'

// FIXME this could possibly go into the PSM QC R code?

  input:
  tuple val(setnames), file(mzmlfiles), val(platenames), path(speclookup), path(oldmzmls_fn), val(fractionation), val(regex_specialchars)

  output:
  tuple path('amount_spectra_files'), path('scans_per_plate'), emit: counted
  path('allplates'), emit: allplates

  script:
  splates = [setnames, platenames].transpose().collect() { "${it[0]}_${it[1]}" }
  """
  #!/usr/bin/env python
  import os
  import re
  import sqlite3
  # FIXME isnt a plate NA?
  # FIXME easier to just parse the --input file?
  con = sqlite3.Connection("$speclookup")
  cursor = con.execute("SELECT set_name, mzmlfilename, COUNT(*) FROM mzml JOIN mzmlfiles USING(mzmlfile_id) JOIN biosets USING(set_id) GROUP BY mzmlfile_id")
  if ${fractionation ? '1' : '0'}:
      platesets = [\"${splates.join('", "')}\"]
      plates = [\"${platenames.join('", "')}\"]
      setnames = [\"${setnames.join('", "')}\"]
      fileplates = {}
      for fn, setname, plate in zip([${mzmlfiles.collect() { "'${stripchars_infile(it)[1]}'" }.join(',')}], setnames, plates):
          try:
              fileplates[fn][setname] = plate
          except KeyError:
              fileplates[fn] = {setname: plate} 
      if ${oldmzmls_fn.name != 'NO__FILE' ? 1 : 0}:
          with open('$oldmzmls_fn') as oldmzfp:
              for line in oldmzfp:
                  if line.strip('\\n').split('\\t')[0] == 'mzmlfile':
                      continue
                  fpath, inst, setname, plate, fraction = line.strip('\\n').split('\\t')
                  # old mzmls also contain files that are no longer used (i.e. removed from set)
                  # filter by skipping any setname that is in the current new setnames
                  if setname in setnames:
                      continue
                  setplate = '{}_{}'.format(setname, plate)
                  fn = re.sub('${regex_specialchars.replaceAll("'", "\\\\'")}', '_', os.path.basename(fpath))
                  if setplate not in platesets:
                      platesets.append(setplate)
                  if fn not in fileplates:
                      fileplates[fn] = {setname: plate}
                  elif setname not in fileplates[fn]:
                      fileplates[fn][setname] = plate
      with open('allplates', 'w') as fp:
          fp.write('\\n'.join(platesets))
      platescans = {p: 0 for p in platesets}
      with open('amount_spectra_files', 'w') as fp:
          for setname, fn, scans in cursor:
              fp.write(f'{setname}\\t{fn}\\t{scans}\\n')
              setplate = f'{setname}_{fileplates[fn][setname]}'
              platescans[setplate] += int(scans)
      with open('scans_per_plate', 'w') as fp:
          for plate, scans in platescans.items():
              fp.write(f'{plate}\\t{scans}\\n')
  else:
      # make sure there's an empty file for scans_per_plate if not fractionated
      with open('amount_spectra_files', 'w') as fp, open('scans_per_plate', 'w') as _, open('allplates', 'w') as __:
          for setname, fn, scans in cursor:
              fp.write(f'{setname}\\t{fn}\\t{scans}\\n')
  """
}


process PSMQC {
    /* Need plotly from conda install 4.10.4, to get right libs and glue
    not available in biocontainers, so using own container.
    */
 
  container params.report_container

  input:
  tuple path('psms'), path('filescans'), path('platescans'), path('mzmldef'), path('oldmzmlfn'), val(fractionation), val(has_newmzmls), val(has_oldmzmls)

  output:
  tuple path('platescans'), path('amount_psms_files'), path("psmplothtml"), path('psmtable_summary.txt')

  script:
  newmzmls = has_newmzmls ? 'mzmldef' : 'FALSE'
  oldmzmls = has_oldmzmls ? 'oldmzmlfn' : 'FALSE'
  """
  sed '1s/\\#SpecFile/SpectraFile/' < psms > psms_clean
  qc_psms.R ${fractionation ? 'TRUE' : 'FALSE'} ${newmzmls} ${oldmzmls}
  mkdir psmplothtml
  mv *.html psmplothtml/
  """
}


process featQC {

  container params.report_container

  input:
  tuple val(acctype), path('feats'), path(normfacs), path(peptable), path('sampletable'), val(setnames), val(has_sampletable), val(conflvl)

  output:
  tuple path(htmldir), path("${acctype}__summary.txt"), path("${acctype}__overlap")

  script:
  parse_normfactors = normfacs.name != 'NO__FILE'
  htmldir = "${acctype}__plothtml"
  """
  ${parse_normfactors ? "cat ${normfacs} > allnormfacs" : ''}
  # Create QC plots and put them base64 into HTML, R also creates summary.txt
  qc_protein.R --sets ${listify(setnames).collect() { "'$it'" }.join(' ')} \
     --feattype ${acctype} --peptable $peptable \
     ${has_sampletable ? "--sampletable $sampletable" : ''} \
     --conflvl $conflvl \
     ${parse_normfactors ? '--normtable allnormfacs' : ''}

  mkdir ${htmldir}
  mv *.html ${htmldir}/
  """
}


process PTMQC {
  
  container params.report_container

  input:
  tuple path(ptmpsms), path(peptable)

  output:
  tuple path('ptmplothtml'), path('*.txt')

  script:
  """
  sed '1s/\\#SpecFile/SpectraFile/' < "${ptmpsms}" > psms_clean
  qc_ptms.R psms_clean "${peptable}"
  mkdir ptmplothtml
  mv *.html ptmplothtml/
  """
}


process summaryReport {
  
  container params.report_container

  input:
  tuple path(platescans), path(plotlibs), path('psmplots'), path(psm_summary), path(featplots), path(feat_summaries), path(feat_overlaps), path('ptmplots'), path(ptmfiles), path('warnings*')
  
  output:
  tuple path('report_groovy_template.html'), path('libs.js')
  
  script:
  has_plates = platescans.size()
  """
  # xargs removes trailing whitespace
  report_tables.py --version "${workflow.manifest.version}" --doi "${workflow.manifest.doi}" \
      --templatedir "$baseDir/assets" \
      ${has_plates ? "--plates \$(cut -f1 $platescans | sort -u | tr '\n' ' ' | xargs)" :''}
  """
}


workflow REPORTING {
  take:
  sets_mzmls_plates
  speclookup
  mzmldef
  oldmzmldef
  fractionation
  plot_psms
  pepprotgenes
  normfacs
  sampletable
  pepconflvl
  prot_gene_conflvl
  setnames
  ptms
  warnings

  main:
  nofile = "${baseDir}/assets/NO__FILE"

  sets_mzmls_plates
  | combine(speclookup)
  | combine(oldmzmldef)
  | map { it + [fractionation, get_regex_specialchars()] }
  | countMS2sPerPlate
  
  plot_psms
  | combine(countMS2sPerPlate.out.counted)
  | combine(mzmldef)
  | combine(oldmzmldef)
  | map { it + [fractionation, it[3].name != 'NO__FILE', it[4].name != 'NO__FILE'] }
  | PSMQC

  pepprotgenes
  | join(normfacs.groupTuple(), remainder: true)
  | map { [it[0], it[1], it[2] ?: nofile] }
  | combine(pepprotgenes.filter { it[0] == 'peptides' }.map { it[1] })
  | combine(sampletable.map { it[1] })
  | map { it + [setnames, it[4].name != 'NO__FILE', it[0] == 'peptides' ? pepconflvl : prot_gene_conflvl] }
  | featQC
  | toList
  | transpose
  | toList
  | set { feat_qc_ch }

  // Only pick one PTM table for reporting nr of sites etc
  ptms
  | filter { it[1].name.contains('_not_adjusted') }
  | PTMQC
  | ifEmpty([nofile, nofile])
  | set { ptmqc }

  PSMQC.out
  | combine(feat_qc_ch)
  | combine(ptmqc)
  | combine(warnings)
  | summaryReport
  
  emit:
  summaryReport.out
}
