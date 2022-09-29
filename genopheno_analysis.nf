nextflow.enable.dsl = 2

//include { EXOMISER; LIRICAL; SYMLINK_PREEXISTING_FILE } from './processes'
include { GENO_PHENO_ANALYSIS } from './main'

//================================================================================
// Read and derive file names and location from the params
//================================================================================

sample_vcf_params = file(params.input_fofn_preexisting_vcfs)


enable_exomiser_batchmode = false
enable_exomiser_batchmode = params.exomiser_batchmode

exomiser_path = file(params.exomiser_path)
exomiser_config = file(params.exomiser_config)
exomiser_yaml_template = file(params.exomiser_yaml_template)
exomiser_output_prefix = params.exomiser_output_prefix

lirical_path =file(params.lirical_path)
lirical_yaml_template = file(params.lirical_yaml_template)
lirical_output_prefix = params.lirical_output_prefix

//================================================================================
// Prepare channels
//================================================================================

iso8601_date = new Date().format("yyyy-MM-dd")

// Convert input manifest to a channel.
sample_params = channel.fromPath(sample_vcf_params)
  .splitText(keepHeader:true)
  .map { line ->
      cols = line.tokenize('\t')
      [
        cols[0], // samplename in vcf file
        cols[1], // samplename in redcap database
        cols[2], // hpo ids (list) ["value_1", "value_2", ...]
        cols[3], // cohort prefix
        cols[4], // sample specific vcf file path
        cols[5] // (main/cohort) output directory path
      ]
  }

sample_vcf_ch = sample_params
  .map { 
    it ->

    // Combine redcap and cohortvcf samplenames separated by "_" if not the same used in output
    if(it[0] == it[1]) { 
      samplename_combo = it[0] 
    } else { 
      samplename_combo = it[1] + '_' + it[0] 
    }

    // Add cohort_prefix infront of samplename_in_vcf if cohort_prefix is not empty
    if (it[3] == ""){
      samplename_in_vcf = it[0]
    }
    else {
      samplename_in_vcf = it[3] + it[0]
    }

    cohort_output_dir = it[5].trim()

    return [ 
            samplename_in_vcf,                                          // samplename in vcf
            file(cohort_output_dir + "/" + samplename_combo + "/" + iso8601_date), // (main/cohort) output directory path
            it[4],                                                      // sample vcf file path
            it[2].replaceAll("^\"|\"\$", "")                            // hpo ids (list)
           ] 
  }

//
//sample_vcf_ch.view()

//================================================================================
// LOGGING
//================================================================================

log.info """\
 G E N O T Y P E - P H E N O T Y P E
      A  N  A  L  Y  S  I  S
 -----------------------------------     
  on single (prepared: vt decomposed
    normalized, bgzipped, tabix 
        indexed) sample vcfs
      with Exomiser and LIRICAL
 ===================================

 """

//================================================================================
// Main workflow
//================================================================================

workflow {

  GENO_PHENO_ANALYSIS( sample_vcf_ch )

}
