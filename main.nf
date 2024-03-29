nextflow.enable.dsl = 2

include { GATK_SELECTVARIANTS; VT_NORMALIZE; VT_DECOMPOSE; VCF_COMPRESS_INDEX; EXOMISER; LIRICAL; SYMLINK_PREEXISTING_FILE; EXOMISER_BATCH_ANALYSIS; GENERATE_EXOMISER_YAML; GENERATE_EXOMISER_BATCH_FILE; SAVE_BATCH_RESULTS_TO } from './processes'

//================================================================================
// Derive file names and location from the params.yaml
//================================================================================

// sample table
sample_params = file(params.input_fofn)

// switch for samplename combo
secondary_first = true
secondary_first = params.secondary_samplename_first

ref_fasta = file(params.fasta)
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta", ".dict"))

vt_path = file(params.vt_path)
vt_genome = file(params.vt_genome)

bgzip_path = file(params.bgzip_path)
tabix_path = file(params.tabix_path)

gatk_path = file(params.gatk_path)

// default is false for exomiser batch analysis mode
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
sample_params_ch = channel.fromPath(sample_params)
      .splitText(keepHeader:true)
      .map { line ->
          cols = line.tokenize('\t')
          [
                cols[0], // samplename in cohort vcf file
                cols[1], // samplename, secondary (for example from redcap database)
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

      // Combine secondary(for example, redcap) and cohortvcf samplenames separated by "_" if not the same used in output
      if(it[0] == it[1]) { 
        samplename_combo = it[0] 
      }
      else if(it[1]) { 
        if(secondary_first == 'true'){
          samplename_combo = it[1] + '_' + it[0] 
        } else {
          samplename_combo = it[0] + '_' + it[1] 
        }
      } else {
        samplename_combo = it[0] 
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
    // ... actually copy these file to designated sample output directory
    // CREATE README.txt TO INFORM ABOUT THIS STEP but only if running genopheno_analysis.nf
    if( workflow.scriptName == 'genopheno_analysis.nf' ){
      
      // edit channel: remove last value (hpo ids)
      sample_info_nohpoids = sample_info 
        .map { 
          it -> return[ 
            it[0], // samplename in vcf
            it[1], // sample output dir
            it[2], // sample vcf file (taken by process as path)
            it[2] + ".tbi", // sample vcf tbi file (taken by process as path)
            it[2] // sample vcf file (taken by process as val)
          ]
        }

      SYMLINK_PREEXISTING_FILE( sample_info_nohpoids )
    }

    //================================================================================
    // Exomiser BATCH ANALYSIS MODE enabled
    //================================================================================
    if( enable_exomiser_batchmode ){

      GENERATE_EXOMISER_YAML(
        sample_info,
        exomiser_yaml_template,
        exomiser_output_prefix
      )

      GENERATE_EXOMISER_BATCH_FILE(
        GENERATE_EXOMISER_YAML.out
      )

      EXOMISER_BATCH_ANALYSIS(

        // Gather all the yaml files
        GENERATE_EXOMISER_BATCH_FILE.out[1].collect(),

        // Gather all the vcf.gz files
        GENERATE_EXOMISER_BATCH_FILE.out[2].collect(),

        // Gather all yaml file (path writes) to single file
        GENERATE_EXOMISER_BATCH_FILE.out[0].collectFile(name: 'exomiser_batch.txt'),
        
        exomiser_path,
        exomiser_config
      )

      // Map outputs to samplenames (and sample_output_dir)
      // Use process publishDir to copy the result files from previous process
      //  to wanted output path

      EXOMISER_BATCH_ANALYSIS.out[0]
        .concat(EXOMISER_BATCH_ANALYSIS.out[1],
                EXOMISER_BATCH_ANALYSIS.out[2],
                EXOMISER_BATCH_ANALYSIS.out[3],
                EXOMISER_BATCH_ANALYSIS.out[4])
        .flatten()
        .map {
          it -> 
          return [ 
            it.baseName.replaceAll(/\..*$/,""),
            it
          ]
        }
        .groupTuple(by: 0)
        .set { batch_out }

      // Trim existing channel to only have
      //  samplename and sample output directory
      samples_dirs = sample_info_nohpoids
        .map { 
          it -> return[ 
            it[0], // samplename in vcf
            it[1], // sample output dir
          ]
        }

      // join channel by samplename
      // transpose so that each file (originally from previous process)
      //  has samplename and output directory
      save_samples = samples_dirs
        .join(batch_out, by: 0)
        .transpose()

      // save result files to previously defined sample output directory
      //   via publishDir
      SAVE_BATCH_RESULTS_TO(
        save_samples
      )


    }else{ // "Normal" non-batch mode

      EXOMISER(
        sample_info,
        exomiser_path,
        exomiser_config, 
        exomiser_output_prefix,
        exomiser_yaml_template
      )
    }

    LIRICAL(
        sample_info,
        lirical_output_prefix,
        lirical_yaml_template,
        lirical_path
    )
}

//================================================================================
// LOGGING
//================================================================================

log.info """\
            M A I N . N F
 G E N O T Y P E - P H E N O T Y P E
      A  N  A  L  Y  S  I  S
 ----------------------------------- 
  Workflow steps:
 -----------------------------------    
  1. extract single sample vcf
  2. vt decompose
  3. vt normalize
  4. compress (bgzip), index (tabix)
  5. Exomiser analysis
  6. LIRICAL analysis
 ===================================

 """

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
