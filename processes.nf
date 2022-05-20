
conda_path = file(params.conda_path)

//================================================================================
// Procesess of the workflow
//================================================================================

process GATK_SELECTVARIANTS {

  conda conda_path

  input:
    tuple val(samplename_in_vcf),
      path(cohort_vcf),
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

process SYMLINK_PREEXISTING_FILE {

  publishDir "$sample_output_dir", mode: 'copy', overwrite: true

  input:
    tuple val(samplename_in_vcf),
          path(sample_output_dir),
          path(sample_vcf),
          path(sample_vcf_index)

  output:
    // tuple val(samplename_in_vcf),
    path("ln_${sample_vcf}")
    path("ln_${sample_vcf_index}")
    path("README_${samplename_in_vcf}.txt")

  script:

  // add index file ending
  //ln -s "${sample_vcf}" "ln_${sample_vcf}"
  //ln -s "${sample_vcf}.tbi" "ln_${sample_vcf}.tbi"

  """
    ln -s "${sample_vcf}" "ln_${sample_vcf}"
    ln -s "${sample_vcf_index}" "ln_${sample_vcf_index}"
    echo "Symbolic link created \$(date) to pre-existing file: ${sample_vcf}" > "README_${samplename_in_vcf}.txt"
    echo "Symbolic link created \$(date) to pre-existing file: ${sample_vcf_index}" >> "README_${samplename_in_vcf}.txt"
  """

}