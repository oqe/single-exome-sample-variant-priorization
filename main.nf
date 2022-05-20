nextflow.enable.dsl = 2

include { GATK_SELECTVARIANTS; VT_NORMALIZE; VT_DECOMPOSE; VCF_COMPRESS_INDEX; EXOMISER; LIRICAL; SYMLINK_PREEXISTING_FILE } from './processes'

//================================================================================
// Derive file names and location from the params.yaml
//================================================================================

// sample table
sample_params = file(params.input_fofn)

ref_fasta = file(params.fasta)
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta", ".dict"))

vt_path = file(params.vt_path)
vt_genome = file(params.vt_genome)

bgzip_path = file(params.bgzip_path)
tabix_path = file(params.tabix_path)

gatk_path = file(params.gatk_path)

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
sample_params_ch = channel.fromPath(sample_params)
      .splitText(keepHeader:true)
      .map { line ->
          cols = line.tokenize('\t')
          [
                cols[0], // samplename in cohort vcf file
                cols[1], // samplename in redcap database
                cols[2], // hpo ids (list) ["value_1", "value_2", ...]
                cols[3], // cohort prefix
                cols[4], // cohort vcf file path
                cols[5] // (main/cohort) output directory path
          ]
      
      }

// Further modify the channel
sample_params_ch_modified = sample_params_ch
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
      // cohort_vcf = it[4]

      return [ 
              samplename_in_vcf,      
              it[2].replaceAll("^\"|\"\$", ""),                           // hpo ids (list)
              it[4],                                                      // cohort vcf file path
              file("$cohort_output_dir/$samplename_combo/$iso8601_date")  // (main/cohort) output directory path
             ] 
    }

// create channel with only samplename in vcf and vcf path
samplename_and_vcf_ch = sample_params_ch_modified
  .map { 
    it ->
    return [
              it[0],  // samplename_in_vcf
              it[2],   // cohort vcf file path
              it[2] + ".tbi" // cohort vcf index file path
            ]
  }


//================================================================================
// Sub-workflow A
//================================================================================

workflow EXTRACT_PREPARE_SINGLE_SAMPLE_VCF {

  take: 
    samplename_and_vcf_ch
  main:
    GATK_SELECTVARIANTS(
    { samplename_and_vcf_ch },
//      cohort_prefix,
//      cohort_vcf,
//      cohort_vcf_index,
    gatk_path,
    ref_dict,
    ref_fasta,
    ref_fasta_fai
    )

    VT_DECOMPOSE(
        GATK_SELECTVARIANTS.out,
        vt_path,
        vt_genome
    )

    VT_NORMALIZE(
        VT_DECOMPOSE.out,
        vt_path,
        vt_genome
    )

    // add sample output directory to channel
    //  samplename_in_vcf, 
    //  output_dir, 
    //  vcf
    // =======================================
    VT_NORMALIZE.out 
      .join(sample_params_ch_modified, by: 0 )
      .map {
        it ->
        return[ 
                it[0], // samplename_in_vcf
                it[4], // output_dir
                it[1]  // vcf file
              ]
      }
      .set { temp_channel }

    // temp_channel.view()

    VCF_COMPRESS_INDEX(
        temp_channel,
        bgzip_path,
        tabix_path
    )

    // add hpos etc to channel
    //  samplename_in_vcf, 
    //  output_dir,
    //  vcf, 
    //  hpo ids
    // ============================================

    VT_NORMALIZE.out 
      .join(sample_params_ch_modified, by: 0 )
      .map {
        it ->
        return[ 
                it[0], // samplename_in_vcf
                it[4], // output_dir
                it[1], // vcf
                it[2]  // hpo_ids
              ] 
      }
      .set { sample_info }

  emit:
    sample_info 
}

//================================================================================
// Sub-workflow B
//================================================================================

workflow GENO_PHENO_ANALYSIS {

  take:
    sample_info

  main:

    // CREATE SYMLINK TO VCF.GZ and VCF.GZ.TBI TO CURRENT SAMPLE OUTPUT DIRECTORY
    // CREATE README.txt TO INFORM ABOUT THIS STEP but only if running genopheno_analysis.nf
    if( workflow.scriptName == 'genopheno_analysis.nf' ){
      
      // edit channel: remove last value (hpo ids)
      sample_info_nohpoids = sample_info 
        .map { 
          it -> return[ it[0], it[1], it[2], it[2] + ".tbi" ]
        }

      SYMLINK_PREEXISTING_FILE( sample_info_nohpoids )
    }

    EXOMISER(
        sample_info,
        exomiser_path,
        exomiser_config, 
        exomiser_output_prefix,
        exomiser_yaml_template
    )

    LIRICAL(
        sample_info,
        lirical_output_prefix,
        lirical_yaml_template,
        lirical_path
    )

}

//================================================================================
// Main workflow
//================================================================================

workflow {

  // WHOLE WORKFLOW - 1: EXTRACT SINGLE SAMPLE AND PREPARE
  //                - 2: GENOME/EXOME + PHENOTYPE ANALYSIS
  // ======================================================

  EXTRACT_PREPARE_SINGLE_SAMPLE_VCF( samplename_and_vcf_ch )
  GENO_PHENO_ANALYSIS( EXTRACT_PREPARE_SINGLE_SAMPLE_VCF.out )

}
