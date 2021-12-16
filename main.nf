nextflow.enable.dsl = 2

//================================================================================
// Read and derive file names and location from the params
//================================================================================

// We want to use predefined files, we want to use files from previous workflow
ref_fasta = file(params.fasta)
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta", ".dict"))

output_dir = file(params.outdir)

cohort_prefix = params.cohort_prefix
cohort_vcf = file(params.cohort_vcf)
cohort_vcf_index = file(params.cohort_vcf_index)

conda_path = file(params.conda_path)

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
// Derive file names and location from the params.yaml
//================================================================================

sample_params = file(params.input_fofn)


//================================================================================
// Procesess of the workflow
//================================================================================

process GATK_SELECTVARIANTS {

  conda conda_path

  input:
    val(samplename_in_vcf)

    val(cohort_prefix)
    path(cohort_vcf)
    path(cohort_vcf_index)

    path(gatk_path)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

  output:
    tuple val(samplename_in_vcf),
      path("${samplename_in_vcf}.vcf")

  script:
  """
    set -euo pipefail

    ${gatk_path} --java-options "${params.gatk_java_opts}" \
      SelectVariants \
      -R ${ref_fasta} \
      -V ${cohort_vcf} \
      -sn ${samplename_in_vcf} \
      -O ${samplename_in_vcf}.vcf
  """

  stub:
  """
    touch "${sample_name_in_vcf}.vcf"
  """

}

process VT_NORMALIZE {

  conda conda_path

  input:
    tuple val(samplename_in_vcf),
      path(input_vcf)

    path(vt_path)
    path(vt_genome)

  output: 
    tuple val(samplename_in_vcf),
      path("${samplename_in_vcf}.vt_normalized.vcf")

  script:
  """
    set -euo pipefail

    ${vt_path} normalize \
    ${input_vcf} \
    -r ${vt_genome} \
    -o ${samplename_in_vcf}.vt_normalized.vcf
  """

  stub:
  """
    touch "${samplename_in_vcf}.vt_normalized.vcf"
  """
}

process VT_DECOMPOSE {

  conda conda_path

  input:
    tuple val(samplename_in_vcf),
      path(input_vcf)

    path(vt_path)
    path(vt_genome)

  output: 
    tuple val(samplename_in_vcf),
      path("${samplename_in_vcf}.vt_decomposed.vcf")

  script:
  """
    set -euo pipefail

    ${vt_path} decompose \
    -s \
    ${input_vcf} \
    -o ${samplename_in_vcf}.vt_decomposed.vcf
  """
  //    -r ${vt_genome} \

  stub:
  """
    touch "${samplename_in_vcf}.vt_decomposed.vcf"
  """
}

process VCF_COMPRESS_INDEX {

  conda conda_path
  publishDir "$sample_output_dir", mode: 'copy', overwrite: true

  input:
    tuple val(samplename_in_vcf),
      val(sample_output_dir),
      path(input_vcf)

    path(bgzip_path)
    path(tabix_path)

  output:
    tuple val(samplename_in_vcf),
      path("${samplename_in_vcf}.vcf.gz"),
      path("${samplename_in_vcf}.vcf.gz.tbi")

  script:
  """
    set -euo pipefail 

    ${bgzip_path} -c ${input_vcf} > ${samplename_in_vcf}.vcf.gz && \
    ${tabix_path} -p vcf ${samplename_in_vcf}.vcf.gz
  """

  stub:
  """
    touch "${samplename_in_vcf}.gz"
    touch "${samplename_in_vcf}.gz.tbi"
  """

}

process EXOMISER {

  conda conda_path
  publishDir "$sample_output_dir", mode: 'copy', overwrite: true

  input: 
    tuple val(samplename_in_vcf),
      val(sample_output_dir),
      path(input_vcf),
      val(hpo_ids)
    
    path(exomiser_path)
    path(exomiser_config)
    val(exomiser_output_prefix)
    path(exomiser_yaml_template)

  output:
    path("${samplename_in_vcf}.exomiser.yaml")
    path("${samplename_in_vcf}.${exomiser_output_prefix}.html")
    path("${samplename_in_vcf}.${exomiser_output_prefix}.json")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_AD.genes.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_AD.variants.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_AD.vcf")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_AR.genes.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_AR.variants.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_AR.vcf")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_XR.genes.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_XR.variants.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_XR.vcf")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_MT.genes.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_MT.variants.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_MT.vcf")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_XD.genes.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_XD.variants.tsv")
    path("${samplename_in_vcf}.${exomiser_output_prefix}_XD.vcf")

  script:

  hpo_ids = hpo_ids.trim()

  """
  set -euo pipefail

  python - <<END
  import yaml
  import os
  import ast

  # Exomiser YAML
  with open("${exomiser_yaml_template}") as f:
    list_doc = yaml.safe_load(f)

  hpos = "${hpo_ids}"
  hpos = ast.literal_eval(hpos)

  list_doc["analysis"]["vcf"] = "${input_vcf}"
  list_doc["analysis"]["hpoIds"] = hpos
  list_doc["analysis"]["proband"] = "${samplename_in_vcf}"
  list_doc["outputOptions"]["outputPrefix"] = "${samplename_in_vcf}.${exomiser_output_prefix}"

  with open("${samplename_in_vcf}.exomiser.yaml", "w") as f:
    yaml.dump(list_doc, f)
  END

  ${params.java_path} ${params.exomiser_java_opts} \
  -jar ${exomiser_path} \
  --analysis "${samplename_in_vcf}.exomiser.yaml" \
  --spring.config.location=${exomiser_config}
  """

  stub:
  """
    touch "${samplename_in_vcf}.exomiser.yaml"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}.html"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}.json"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_AD.genes.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_AD.variants.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_AD.vcf"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_AR.genes.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_AR.variants.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_AR.vcf"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_XR.genes.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_XR.variants.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_XR.vcf"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_MT.genes.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_MT.variants.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_MT.vcf"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_XD.genes.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_XD.variants.tsv"
    touch "${samplename_in_vcf}.${exomiser_output_prefix}_XD.vcf"
  """

}


process LIRICAL {
  
  conda conda_path
  publishDir "$sample_output_dir", mode: 'copy', overwrite: true

  input: 
    tuple val(samplename_in_vcf),
      val(sample_output_dir),
      path(input_vcf),
      val(hpo_ids)

    val(lirical_output_prefix)
    path(lirical_yaml_template)
    path(lirical_path)

  output:
    path("${samplename_in_vcf}.lirical.yaml")
    path("${samplename_in_vcf}.${lirical_output_prefix}.html")
    path("lirical.log")

  script:
  // trim line breaks
  hpo_ids = hpo_ids.trim()
  
  """
  set -euo pipefail

  python - <<END

  import yaml
  import ast
  import sys

  # LIRICAL YAML
  with open("${lirical_yaml_template}") as f:
    list_doc = yaml.safe_load(f)

  hpos = "${hpo_ids}"

  hpos = ast.literal_eval(hpos)

  list_doc["analysis"]["vcf"] = "${input_vcf}"
  list_doc["hpoIds"] = hpos
  list_doc["prefix"] = "${samplename_in_vcf}.${lirical_output_prefix}"
  
  with open("${samplename_in_vcf}.lirical.yaml", "w") as f:
    yaml.dump(list_doc, f)
  END

  #!/bin/bash

  ${params.java_path} ${params.lirical_java_opts} \
  -jar ${lirical_path} \
  yaml -y "${samplename_in_vcf}.lirical.yaml"

  """

  stub:
  """
    touch "${samplename_in_vcf}.lirical.yaml"
    touch "${samplename_in_vcf}.${lirical_output_prefix}"
    touch "lirical.log"
  """

}



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
                cols[2], // hpo ids (list)
          ]
      }
    .map { 
      it ->
      // Combine redcap and cohortvcf samplenames separated by "_" if not the same
      // ==========================================================================
      if(it[0] == it[1]) { 
        samplename_combo = it[0] 
      } else { 
        samplename_combo = it[1] + '_' + it[0] 
      }
      // Add cohort_prefix infront of samplename_in_vcf if cohort_prefix is not empty
      // =============================================================================
      if (cohort_prefix == ""){
        samplename_in_vcf = it[0]
      }
      else {
        samplename_in_vcf = cohort_prefix + it[0]
      }
      //            0,              1,     2,     3,                                4,                                                                       
      return [ samplename_in_vcf, it[0], it[1], it[2].replaceAll("^\"|\"\$", ""),  file("$output_dir/$samplename_combo/$iso8601_date") ] 
    }

// create channel with only samplename in vcf
vcf_samplenames_ch = sample_params_ch
  .map { 
    it ->
    return it[0]
  }

vcf_samplenames_ch.view()

//================================================================================
// Main workflow
//================================================================================

workflow {

  GATK_SELECTVARIANTS(
      { vcf_samplenames_ch },
      cohort_prefix,
      cohort_vcf,
      cohort_vcf_index,
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
  // samplename_in_vcf, output_dir, vcf
  // =======================================
  VT_NORMALIZE.out 
    .join(sample_params_ch, by: 0 )
    .map {
      it ->
      return[ it[0], it[5], it[1] ]
    }
    .set { temp_channel }

  //temp_channel.view()

  VCF_COMPRESS_INDEX(
      temp_channel,
      bgzip_path,
      tabix_path
  )

  // add hpos etc to channel
  // samplename_in_vcf, output_dir, vcf, hpo ids
  // ============================================
  VT_NORMALIZE.out 
    .join(sample_params_ch, by: 0 )
    .map {
      it ->
      return[ it[0], it[5], it[1], it[4]]
    }
    .set { sample_info }

  //sample_info.view()

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
