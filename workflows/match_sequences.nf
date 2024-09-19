include { get_field_nr; listify; get_regex_specialchars} from '../modules.nf' 

process createTrypticMatchDB {
  /* Create a sequence match SQLite database to be used for filtering
  PSM tables and percolator */

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple path(sequences), val(maxmiscleav), val(minpeplen)

  output:
  path(outfile)

  script:
  miscleav = maxmiscleav < 0 ? 2 : maxmiscleav
  outfile = sequences.baseName.replaceAll('[^a-zA-Z0-9_-]', '_')
  """
  msstitch storeseq -i ${sequences} -o ${outfile} \
    --minlen ${minpeplen} --cutproline \
    --nterm-meth-loss --map-accessions \
    ${miscleav ? "--miscleav ${miscleav}" : ''}
  """
}


process markPeptidesPresentInDB {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  /* Match peptide sequences to user-provided sequences in a storeseq SQLite DB.
  Mark peptides that match with a 1, else 0 in a new column.
  Field name for column will be user-provided 
  */

  input:
  tuple path(peptides), path(seqdbs)

  output:
  tuple val('peptides'), path(peptide_table_out), emit: peptides
  tuple path(jointable_out), emit: jointable

  script:
  peptide_table_out = 'peptide_table.txt'
  jointable_out = 'seqmatch_jointable'
  len_seqmatch = listify(seqdbs).size()
  """
  # Get bare peptide (works also with -f1)
  echo 'Peptide\tMSGFScore\tEValue' > pepseqs
  cut -f2 ${peptides} |  awk -v FS='\\t' -v OFS='\\t' '{print \$1, "NA", "NA"}' | tail -n+2 >> pepseqs
  for seqdb in ${listify(seqdbs).join(' ')}
    do msstitch seqmatch -i pepseqs -o tmppeps --dbfile \$seqdb --matchcolname \${seqdb}
    mv tmppeps pepseqs
  done
  paste <(cut -f4-${3+len_seqmatch} pepseqs) ${peptides} > ${jointable_out}
  paste ${peptides} <(cut -f4-${3+len_seqmatch} pepseqs) > ${peptide_table_out}
  """
}


process joinAnnotatedSeqmatchPeptides {

  tag 'sqlite'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(acctype), path('feats'), path('peptides'), path(seqdbs)

  output:
  tuple val(acctype), path(outfile)
  
  script:
  headfield = [proteins: "Protein(s)", ensg: "Gene ID(s)", genes: "Gene name(s)"][acctype]
  featfield = [proteins: "Protein ID", ensg: "Gene ID", genes: "Gene Name"][acctype]
  len_seqmatch = listify(seqdbs).size()
  seqdbs_clean = seqdbs.collect() { it.baseName.replaceAll('[^A-Za-z0-9_-]', '_') }
  outfile = "${acctype}_table.txt"
  """
  acc_col=${get_field_nr('peptides', headfield)}
  cut -f1-${len_seqmatch},\${acc_col} peptides > seqmatches
  echo "\$(head -n1 feats)\t${seqdbs_clean.join('\t')}" > joinedfeats
  sqlite3 test.db <<END_COMMANDS >> ${outfile}
.mode tabs
.import seqmatches seqmatch
.import feats feats

SELECT feats.*, ${listify(seqdbs_clean).collect() {"sqm_group.${it}"}.join(',')} FROM feats
JOIN (
  SELECT "${headfield}", ${listify(seqdbs_clean).collect() {"GROUP_CONCAT(DISTINCT $it) AS $it"}.join(',')}
  FROM seqmatch GROUP BY "${headfield}") AS sqm_group
ON sqm_group."${headfield}"=feats."${featfield}"
END_COMMANDS
  """
}


workflow MATCH_SEQUENCES {
  take:
  peptable
  otherfeattables
  sequence_fa_ch
  maxmiscleav
  minpeplen
  
  main:
  sequence_fa_ch
  | map { [it, maxmiscleav, minpeplen] }
  | createTrypticMatchDB

  peptable
  | combine(createTrypticMatchDB.out)
  | groupTuple
  | markPeptidesPresentInDB

  otherfeattables
  | combine(markPeptidesPresentInDB.out.jointable)
  | combine(sequence_fa_ch.toList().toList())
  | joinAnnotatedSeqmatchPeptides
  | concat(markPeptidesPresentInDB.out.peptides)
  | set { matched_out }

  emit:
  matched_out
}
