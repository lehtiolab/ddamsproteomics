/* Reporting workflow for DDA pipeline */

include { listify; stripchars_infile; get_regex_specialchars } from '../modules.nf' 

process countMS2sPerPlate {
  container "python:3.12"

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
              setplate = '{}_{}'.format(setname, fileplates[fn][setname])
              platescans[setplate] += int(scans)
      with open('scans_per_plate', 'w') as fp:
          for plate, scans in platescans.items():
              fp.write('{}\\t{}\\n'.format(plate, scans))
  else:
      # make sure there's an empty file for scans_per_plate if not fractionated
      with open('amount_spectra_files', 'w') as fp, open('scans_per_plate', 'w') as _, open('allplates', 'w') as __:
          for line in cursor:
              fp.write('\\t'.join(line))
              fp.write('\\n')
  """
}


process PSMQC {
    /* Need plotly from conda install 4.10.4, to get right libs and glue
    not available in biocontainers, so using own container.
    */
  container 'lehtiolab/dda_report'

  input:
  tuple path('psms'), path('filescans'), path('platescans'), path('mzmlfrs'), path(oldmzmlfn), val(fractionation)

  output:
  tuple path(platescans), path("psmplothtml")

  script:
  """
  qc_psms.R ${fractionation ? 'TRUE' : 'FALSE'} ${oldmzmlfn}
  """
}


process featQC {
  container 'lehtiolab/dda_report'

  input:
  tuple val(acctype), path('feats'), path(normfacs), path(peptable), path('sampletable'), val(setnames), val(has_sampletable), val(conflvl)

  output:
  //tuple path(platescans), path("psmplothtml")

  script:
  parse_normfactors = normfacs.name != 'NO__FILE'
  """
  ${parse_normfactors ? "cat ${normfacs} > allnormfacs" : ''}
  # Create QC plots and put them base64 into HTML, R also creates summary.txt
  qc_protein.R --sets ${listify(setnames).collect() { "'$it'" }.join(' ')} \
     --feattype ${acctype} --peptable $peptable \
     ${has_sampletable ? "--sampletable $sampletable" : ''} \
     --conflvl $conflvl \
     ${parse_normfactors ? '--normtable allnormfacs' : ''}

  mkdir psmplothtml
  mv *.html psmplothtml/
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
  | map { it + fractionation }
  | PSMQC

  pepprotgenes
  | join(normfacs.groupTuple(), remainder: true)
  | map { [it[0], it[1], it[2] ?: nofile] }
  | combine(pepprotgenes.filter { it[0] == 'peptides' }.map { it[1] })
  | combine(sampletable.map { it[1] })
  | combine(setnames)
  | map { it + [it[4].name != 'NO__FILE', it[0] == 'peptides' ? pepconflvl : prot_gene_conflvl] }
  | featQC

}
