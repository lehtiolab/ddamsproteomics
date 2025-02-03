def get_regex_specialchars() {
  return '[&<>\'"]'
}

def get_field_nr(fn, fieldname) {
    return "\$(head -n1 ${fn} | tr '\\t' '\\n' | grep -wn '^${fieldname}\$' | cut -f 1 -d':')"
}

def get_field_nr_multi(fn, fieldnames) {
    /* return field nrs comma separated like: 1,2,5,9 */
    return "\$(head -n1 ${fn} | tr '\\t' '\\n' | grep -En '(${fieldnames.join('|')})' | cut -f 1 -d':' | tr '\\n' ',' | sed 's/\\,\$//')"
}

def get_complement_field_nr(fn, fieldname) {
    /* return field nrs NOT matching input, comma separated like: 1,2,5,9 */
    return "\$(head -n1 ${fn} | tr '\\t' '\\n' | grep -vwn '^${fieldname}\$' | cut -f 1 -d':' | tr '\\n' ',' | sed 's/\\,\$//')"
}

def parse_isotype(isobtype) {
  return ['tmt16plex', 'tmt18plex'].any { it == isobtype } ? 'tmtpro' : isobtype
}

def listify(it) {
  /* This function is useful when needing a list even when having a single item
  - Single items in channels get unpacked from a list
  - Processes expect lists. Even though it would be fine
  without a list, for single-item-lists any special characters are not escaped by NF
  in the script, which leads to errors. See:
  https://github.com/nextflow-io/nextflow/discussions/4240
  */
  return it instanceof java.util.List ? it : [it]
}


def stripchars_infile(x, return_oldfile=false) {
  // FIXME %, ?, * fn turns a basename into a list because they are wildcards
  // Replace special characters since they cause trouble in percolator XML output downstream
  // e.g. & is not allowed in XML apparently (percolator doesnt encode it)
  // and LXML then crashes on reading it.
  // Also NF doesnt quote e.g. semicolons it seems.
  def scriptinfile = "${x.baseName}.${x.extension}"
  def parsed_file = scriptinfile.replaceAll(get_regex_specialchars(), '_')
  if (return_oldfile) {
    return [parsed_file != scriptinfile, parsed_file, scriptinfile]
  } else {
    return [parsed_file != scriptinfile, parsed_file]
  }
}

def read_header(info_fn) {
  def header = []
  def info = file(info_fn).eachLine { line, ix ->
      if (ix == 1) {
        header = line.tokenize('\t')
      }
  }
  return header
}

def create_info_map(info_fn, possible_params) {
  /* From possible params this parses a tab separated input file with
  a header (which has some of those params names.
  It returns a map with params and their values, which are defaults for
  those not set in the input file.
  */
  def info = file(info_fn).readLines().collect { it.tokenize('\t') }
  def header = info.pop()

  def params_not_header = possible_params - header

  def info_map = [:]
  def empties = [null, '', 'none']

  info.findAll{ it[0][0] != '#' }.eachWithIndex { item, ix ->
    info_map[ix] = [:]
    header.eachWithIndex{ hfield, hix ->
      if (empties.any { it == item[hix] }) {
        info_map[ix][hfield] = ''
      } else {
        info_map[ix][hfield] = item[hix]
      }
    }
    info_map[ix].id = ix
    info_map[ix].mzmlfile = file(info_map[ix].mzmlfile)
    info_map[ix].filename = "${info_map[ix].mzmlfile.baseName}.${info_map[ix].mzmlfile.extension}"
    info_map[ix].fn_normalized_chars = info_map[ix].filename.replaceAll(get_regex_specialchars(), '_')
    info_map[ix].setname = info_map[ix].setname.replaceAll('[ ]+$', '').replaceAll('^[ ]+', '')
    info_map[ix].sample = info_map[ix].mzmlfile.baseName.replaceAll(get_regex_specialchars(), '_')
    params_not_header.each {
      info_map[ix][it] = params[it]
    }
  }
  return info_map
}


def msgf_info_map(info_fn) {
  expected_fields = ["mzmlfile", "setname", "plate", "fraction", "phospho", "activation",
        "prectol", "iso_err", "instrument", "enzyme", "frag", "terminicleaved"]
  def samples = create_info_map(info_fn, expected_fields)
  return samples
}


process createMods {

  tag 'python'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), val(isobtype), val(maxvarmods), val(search_engine), path(msgfmods), val(mods)

  output:
  tuple val(setname), path('mods.txt')

  script:
  """
  create_modfile.py $maxvarmods "${msgfmods}" $search_engine "${mods}${isobtype ? ";${parse_isotype(isobtype)}" : ''}"
  """
}
