# Single Exome Sample Variant Priorization

This is a nextflow[^nextflow] workflow for extracting (human) single sample from a larger exome sequencing multisample/cohort VCF and annote and prioritize exome variants of that sample. "Single sample" as in the context of the workflow stages per sample. You may execute the workflow for multiple samples in parallel.

![Analysis workflow figure](/extras/images/workflow_figure.png)
*<p class="text-center">Figure 1. Analysis Workflow simplified flowchart. Does not include samplename preparation (samplename combination, cohort prefix), nor all process input variables/channel values in the actual (nextflow) workflow.</p>*

## Workflow

**Here is a short description of the workflow stages (main.nf):**  

0. Generate samplename combination (depending on alternative samplename and possible cohort prefix) and date specifc output directories
1. Single sample vcf file is extracted from the cohort vcf with GATK SelectVariants (-sn) tool [^gatk].
2. Extracted vcf is decomposed   
3. ..and normalized by vt tool[^vt].
4. Decomposed and normalized vcf is compressed (.gz) using bgzip and indexed (.tbi) using tabix[^htslib].
5. Generate sample specific .yaml files[^python]<sup>,</sup>[^pyyaml] for Exomiser from input table information[^exomiser]. Annotation and priorization for the variants of the vcf with sample specific HPO terms is done with Exomiser [^hpo].
6. Generate sample specific .yaml files for LIRICAL from input table information and "Phenotype-driven priorization of candidate diseases and genes..." is done with LIRICAL (: LIkelihood Ratio Interpretation of Clinical AbnormaLities) [^lirical].  

**Alternative workflow for prepared sample vcf(s) workflow stages (genopheno_analysis.nf):**  
Stages 5. and 6. from previous list. VCFs are expected to be prepared(vt decompose, vt normalized, compressed and indexed).

### Why?

Streamlining individual patient genetic analysis from large cohort of (exome) sequenced patients.

Nextflow for fast deployment with good support for SLURM and cloud service providers[^nextflow]. Conda for easily staying up-to-date, easy management and deployment of available software (versions, for some of the software used)[^conda]. 

Exomiser and LIRICAL for applying phenotype information with genetic information for variant priorization and helping in diagnosis effort with human phenotype ontology terms (HPO).

Single sample vcf extraction is done for the purpose of speeding up analysis applications down the line. Applications don't therefore need to read large cohort vcf file(s) for each step.

**genopheno_analysis.nf** can be used for pre-existing sample vcfs, that might be from an earlier cohort vcf extraction workflow run or for other individual samples, for example from commercial sequencing services.

### 1. Requirements

<details><summary> Software </summary><blockquote>

+ \*nix based system (or WSL in Windows(10/11))  
	+ Nextflow  
	+ conda (conda enviroment .yaml is provided)  
		- java / Openjdk 11.0  
		- python  
			- pyyaml package  
		- gatk  
		- vt  
			- vt reference genome  
		- bgzip  
		- tabix  
	+ Exomiser 13.0.0 or newer  
		- with Exomiser data downloaded and configured to work  
		- application.properties configured  
	+ LIRICAL   
		- with LIRICAL data downloaded and configured to work  
</blockquote></details>

### 2. Input

Main workflow input is defined in **params.yaml**.  

Sample input table file (.tsv) for **main.nf** workflow (path given in params.yaml) is **input_fofn**.  
Sample input table file (.tsv) for **genopheno_analysis** (path given in params.yaml) is **input_fofn_preexisting_vcfs**.

See more details below.

#### 2.1. params.yaml details

Edit **params.yaml** for your inputs:

<details><summary> Sample input table(s) </summary><blockquote>
<details><summary> input_fofn </summary><blockquote>

Path to sample input table file. This is the main sample specific input file for **main.nf** workflow.

Table can have as many samples from any cohort you wish.

Sample input table where value (columns) are tab separated (.tsv), with following variables per sample per line. Columns:
1. required - unique id / sample name in cohort VCF file, with or without cohort/sub-cohort designation  
2. optional - secondary samplename (for example in REDCap/database where HPO terms were extracted), can be left empty
3. required - list of HPO terms per sample, in format: ['HPO:XXXXX', 'HPO:XXXXX']  
4. optional - cohort_prefix - if your samplename has cohort prefix in the vcf file you can add it here. you may also leave this empty but then the full samplename matching vcf samplename needs to be input in the first column  
5. required - cohort_vcf - path to cohort vcf file  
6. required - output_dir - base path where you want sample specific results to be stored (samplename combination and date of workflow execution are added to base path as subfolders)

**Example input_fofn table** 

| samplename_cohort_vcf | samplename_secondary | hpo_ids | cohort_prefix | cohort_vcf | output_dir |  
| --- | --- | --- | --- | --- | --- |  
| X1 | XRX200 | ['HPO:XXXX', 'HPO:YYYY', 'HPO:ZZZZZ'] | studyABC_ | /path/to/studyABC/data/studyABC.vcf | /path/to/studyABC/results/individual |  
| X2 | XR300 | ['HPO:X1000', 'HPO:00120', 'HPO:ZZZZZ'] | studyABC_ | /path/to/studyABC/data/studyABC_v2.vcf | /path/to/studyABC/results/individual |  
| Z001 | 2345 | ['HPO:00042', 'HPO:000124', 'HPO:00139'] |  | /path/to/studyXYZ/data/studyXYZ.vcf | /path/to/studyXYZ/results/ |  
| Z002 | Z002 | ['HPO:000252', 'HPO:00004'] |  | /path/to/studyXYZ/data/studyXYZ.vcf | /path/to/studyZ/results/individual |  
| COHORT_Y_Y02 |  | ['HPO:00042', 'HPO:000124', 'HPO:00139'] |  | /path/to/studyXYZ/data/studyXYZ.vcf | /path/to/studyY/results/patient |  

Note that header (first row) of the .tsv table is not used. Column names in the header are for reference when fillin in a table. Header row is required, otherwise table will be missing the first row of data.

**NAMING LOGIC**  
**Samplename combination**:
- 1st row: In the example input_fofn table 1st sample in the first row, X1 actual samplename in the vcf file is **studyABC_X1**. As the cohort prefix is defined as **studyABC_**. Samplename combination would in this case be **X1_XRX200**.  
- 4th row: Samplename combination is **Z002**
- 5th row: Samplename combination is **COHORT_Y_Y02**

**Samplename in vcf**:
- 1st row: Samplename in vcf is **studyABC_X1**
- 3rd row: Samplename in vcf is **Z001**
- 5th row: Samplename in vcf is **COHORT_Y_Y02**

**Sample output directory** (for results on a date 2022-05-24):
- 1st row: **/path/to/studyABC/results/individual/XRX200_X1/2022-05-24**
- 3rd row: **path/to/studyABC/results/2345_Z001/2022-05-24**
- 4th row: **/path/to/studyXYZ/results/individual/Z002/2022-05-24**
- 5th row: **/path/to/studyY/results/patient/COHORT_Y_Y02/2022-05-24**

As the samplenames can be different in a vcf file and in redcap(or other equivalent database) or sample might have secondary identifier there is a column for secondary (or REDCap specific) samplename and cohort-VCF file specific samplename. Cohort-VFC samplename can be full samplename and then the cohort prefix variable can be left empty (""). Cohort-VFC samplename can also be unique part of samplename if samplenames have cohort specific prefix name. In this case you need to fill in cohort prefix variable so that the cohort_prefix + samplename_cohort_vcf create a valid samplename that is in the vcf file.
</blockquote></details>

<details><summary> input_fofn_preexisting_vcfs </summary><blockquote>

Path to sample input table file. This is the sample specific input file for **genopheno_analysis.nf** workflow.

**input_fofn_preexisting_vcfs**

See above, input_fofn for more details. Same logic applies to input_fofn_preexisting_vcfs.

| samplename_vcf | samplename_secondary | hpo_ids | cohort_prefix | sample_vcf | output_dir |  
| --- | --- | --- | --- | --- | --- |  
| X1 | XRX200 | ['HPO:XXXX', 'HPO:YYYY', 'HPO:ZZZZZ'] | studyABC_ | /path/to/studyABC/results/individual/studyABC_X1_XRX200.vcf | /path/to/studyABC/results/individual |  
</blockquote></details>
</blockquote></details>

<details><summary> secondary_samplename_first </summary><blockquote>

Switch true/false whether to have secondary samplename (input_fofn... 2nd column) first in samplename combination (samplename_combo) subfolder name in output path for results.

</blockquote></details>

<details><summary> Applications, settings, enviroment and settings related input </summary><blockquote>
<details><summary> Exomiser related input </summary><blockquote>
<details><summary> exomiser_path </summary><blockquote>

Path to exomiser jar file  
</blockquote></details>

<details><summary> exomiser_config </summary><blockquote>

Path to **application.properties** file of exomiser.  
</blockquote></details>

<details><summary> exomiser_yaml_template </summary><blockquote>

Path to exomiser yaml analysis file.  

Template is provided in /templates/exomiser_13.0.0_hg38_template.yaml, edit this file with your values.
</blockquote></details>

<details><summary> exomiser_output_prefix </summary><blockquote>

String output prefix for exomiser result file.
</blockquote></details>
</blockquote></details>

<details><summary> LIRICAL related input </summary><blockquote>
<details><summary> lirical_path </summary><blockquote>

Path to LIRICAL.jar file.
</blockquote></details>

<details><summary> lirical_yaml_template </summary><blockquote>
Path to lirical yaml analysis file.  
  
Template is provided in /templates/lirical_1.3.4_template.yaml, edit this file with your values and fill in exomiser data directory path in the template.
</blockquote></details>

<details><summary> lirical_output_prefix </summary><blockquote>

String output prefix for rirical result file.
</blockquote></details>
</blockquote></details>

<details><summary> CONDA related input </summary><blockquote>
<details><summary> conda_path </summary><blockquote>
Path to your conda enviroment.  
You can import the conda enviroment from /templates/conda_env_variant-tools.yml  
</blockquote></details>

<details><summary> vt_path </summary><blockquote>

Path to vt, default provided as per conda env: vt  
</blockquote></details>

<details><summary> vt_genome </summary><blockquote>
Path to vt generated reference genome.
</blockquote></details>

<details><summary> bgzip_path </summary><blockquote>

Path to bgzip, default is provided as per conda env: bgzip 
</blockquote></details>

<details><summary> tabix_path </summary><blockquote>

Path to bgzip, default is provided as per conda env: tabix  
</blockquote></details>

<details><summary> gatk_path </summary><blockquote>

Path to gatk, default is provided as per conda env: gatk
</blockquote></details>
</blockquote></details>

<details><summary> Reference (genome) input </summary><blockquote>
<details><summary> fasta </summary><blockquote>

Path to reference genome fasta file.  
  
!Note .fai and .dict files are also required, but their paths are generated based on fasta file path.
</blockquote></details>
</blockquote></details>

</blockquote></details>

#### 2.2. nextflow.config

Edit your **nextflow.config** according to how and where you plan to execute the workflow. For example fill in SLURM server parameters.

### 3. Execution

Note that **main.nf** requires **input_fofn** to be defined. **input_fofn_preexisting_vcfs** is ignored.

```bash
nextflow run main.nf \\
-params-file params.main.2022-05-24.yaml \\
-profile local \\
-c nextflow.config
```

Note that **genopheno_analysis.nf** requires **input_fofn_preexisting_vcfs** to be defined. **input_fofn** is ignored.

```bash
nextflow run genopheno_analysis.nf \\
-params-file params.genopheno.2022-05-24.yaml \\
-profile local \\
-c nextflow.config
```

Alternatively you can give any variable defined in **params.yaml** as specified parameter with nextflow run...

```bash
nextflow run main.nf \\
-params-file params.2022-05-24.yaml \\
-profile local \\
-c nextflow.config
--input_fofn /my/research/project/cohort_abc/data/nf_run_2022-05-24.tsv
```

This way you might reuse your **params.yaml** file if there are no other changes to variables needed.

## 4. Output

Output file base folder path is given in **input_fofn** or **input_fofn_preexisting_vcfs** table files per each (row) sample in the 6th column. Samplename combination of input_fofn file's 1st column (samplename_in_vcf) and 2nd column (samplename_redcap) is designated as subfolder and (workflow execution) date (yyyy-MM-dd) subfolder are added to base output path.

	<input_main_output_dir> / <samplename_combo > / <date_in_yyyy-MM-dd>

Results (of specified process) outputs are saved to path described above.

### 4.1. main.nf

After the workflow has completed the following files should be in each date specific folder inside sample name folder:  

- (vt decomposed and normalized) bgzip compressed .vcf.gz and tabix indexed .vcf.gz.tbi  
- exomiser yaml analysis file  
- exomiser output files (many files)   
- lirical yaml analysis file  
- lirical output file  

### 4.2. genopheno_analysis.nf

- symbolic link to used sample vcf.gz and vcf.gz.tbi file
- exomiser yaml analysis file  
- exomiser output files (many files)  
- lirical yaml analysis file  
- lirical output file  

## 5. Extras

**_**

**wdl workflow**

**input_fofn_importer.R**  

HPO term extraction per samplename is separated as an external R Script[^r]<sup>,</sup>[^r-redcap]. R script is for sample(name) exploration/checking from REDCap database as well as extracting HPO terms and generating sample input table required as input in this workflow. Expectation is that in your use case your institution/or your affiliates are hosting REDCap server and have HPO terms saved there[^redcap1]<sup>,</sup>[^redcap2]. Inputting the HPO terms per patient, manually into *input_fofn* table is of course possible.

## TO-DO

- [ ] Document zero-to-finish example of this workflow to a separate file  
- [ ] Finalize extra R script (REDCap HPO-terms extraction)  
	- [ ] Write up documentation  
- [ ] Add original .wdl workflow in extras/wdl folder
- [ ] Add VEP annotation?
- [ ] Clean-up, modulize code?  
- [ ] Inbedding sample specific .bam/.cram file to exomiser output .html with IGV.js  
- [ ] Docker/Apptainer formerly known as singularity) container for programs ???  
- [x] create simple graphic/graph about the workflow  

- [x] Replace samplename_in_redcap with samplename_secondary (input_fofn, 2nd column)
	- [x] Explain that you can leave secondary samplename empty, check and modify the source code for this
	- [x] ? add a switch to to params.yaml to select if samplename_in_vcf or samplename_secondary comes first in the samplename_combo
- [x] Update conda enviroment due to Log4 vulnerability + check and test openjdk  
~~- [ ] Add checking if sample specific .vcf (.gz and .tbi) files are found in the sample output directory~~  
	~~- [ ] Add switch to utilize found .vcf file and skip first step (GATK_SELECTVARIANTS) to save time~~  
- [x] Change/add sample specific parametres to **input_fofn** table...  
	- [x] output directory
	- [x] cohort-vcf  
	*enables executing sample variant priorization for samples from different cohort or cohort-vcf versions in with the same input_fofn table file*
- [x] Gather and write-out references to this README.md and add references to the end  
	- [x] nextflow  
	- [x] conda  
	- [x] gatk  
	- [x] htslib: bgzip,tabix  
	- [x] vt  
	- [x] exomiser  
		- [x] HPO  
	- [x] lirical  
	- [x] R  
		- [x] REDCap package  
	- [x] python  
		- [x] yaml package / pyyaml

## Tool links

* **GATK SelectVariants**: [https://gatk.broadinstitute.org/hc/en-us/articles/4414594350619-SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/4414594350619-SelectVariants)  
* **Vt - normalize, decompose**: [https://genome.sph.umich.edu/wiki/](https://genome.sph.umich.edu/wiki/)  
* **HTSlib - _bgzip, tabix_**:  
	- [http://www.htslib.org/doc/bgzip.html](http://www.htslib.org/doc/bgzip.html)  
	- [http://www.htslib.org/doc/tabix.html](http://www.htslib.org/doc/tabix.html)  
* **HPO**: [https://hpo.jax.org/app/](https://hpo.jax.org/app/)  
* **Exomiser**: [https://github.com/exomiser/Exomiser](https://github.com/exomiser/Exomiser)  
* **LIRICAL**:   
	- [https://github.com/TheJacksonLaboratory/LIRICAL](https://github.com/TheJacksonLaboratory/LIRICAL)  
	- [https://lirical.readthedocs.io/en/latest/](https://lirical.readthedocs.io/en/latest/)  
* **nextflow**: [https://www.nextflow.io/](https://www.nextflow.io/)  
* **conda**: [https://docs.conda.io/en/latest/](https://docs.conda.io/en/latest/)  
* **PyYAML**: [https://pyyaml.org/wiki/PyYAML](https://pyyaml.org/wiki/PyYAML)  

## References

[^gatk]: Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.  
[^vt]: Tan, A., Abecasis, G. R., & Kang, H. M. (2015). Unified representation of genetic variants. Bioinformatics (Oxford, England), 31(13), 2202–2204. [https://doi.org/10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112)  
[^htslib]: BBonfield, J. K., Marshall, J., Danecek, P., Li, H., Ohan, V., Whitwham, A., Keane, T., & Davies, R. M. (2021). HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience, 10(2), giab007. [https://doi.org/10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007)
[^python]: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.
[^pyyaml]: Simonov Kirill, YAML and Python communities (2020). PyYAML. Retrieved from [https://pyyaml.org/wiki/PyYAML](https://pyyaml.org/wiki/PyYAML)
[^exomiser]: Robinson PN, Köhler S, Oellrich A, et al. Improved exome prioritization of disease genes through cross-species phenotype comparison. Genome Research. 2014 Feb;24(2):340-348. [DOI: 10.1101/gr.160325.113](http://doi.org/10.1101/gr.160325.113). PMID: 24162188; PMCID: PMC3912424. 
[^hpo]: Köhler, S., Gargano, M., Matentzoglu, N., Carmody, L. C., Lewis-Smith, D., Vasilevsky, N. A., Danis, D., Balagura, G., Baynam, G., Brower, A. M., Callahan, T. J., Chute, C. G., Est, J. L., Galer, P. D., Ganesan, S., Griese, M., Haimel, M., Pazmandi, J., Hanauer, M., Harris, N. L., … Robinson, P. N. (2021). The Human Phenotype Ontology in 2021. Nucleic acids research, 49(D1), D1207–D1217. [https://doi.org/10.1093/nar/gkaa1043](https://doi.org/10.1093/nar/gkaa1043)  
[^lirical]: Robinson, P. N., Ravanmehr, V., Jacobsen, J., Danis, D., Zhang, X. A., Carmody, L. C., Gargano, M. A., Thaxton, C. L., UNC Biocuration Core, Karlebach, G., Reese, J., Holtgrewe, M., Köhler, S., McMurry, J. A., Haendel, M. A., & Smedley, D. (2020). Interpretable Clinical Genomics with a Likelihood Ratio Paradigm. American journal of human genetics, 107(3), 403–417. [https://doi.org/10.1016/j.ajhg.2020.06.021](https://doi.org/10.1016/j.ajhg.2020.06.021)  
[^r]: R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL [https://www.R-project.org/](https://www.R-project.org/).  
[^r-redcap]: Will Beasley (2021). REDCapR: Interaction Between R and REDCap. R package version 1.0.0. [https://CRAN.R-project.org/package=REDCapR](https://CRAN.R-project.org/package=REDCapR)  
[^redcap1]: Harris, P. A., Taylor, R., Thielke, R., Payne, J., Gonzalez, N., & Conde, J. G. (2009). Research electronic data capture (REDCap)--a metadata-driven methodology and workflow process for providing translational research informatics support. Journal of biomedical informatics, 42(2), 377–381. [https://doi.org/10.1016/j.jbi.2008.08.010](https://doi.org/10.1016/j.jbi.2008.08.010)  
[^redcap2]: Harris, P. A., Taylor, R., Minor, B. L., Elliott, V., Fernandez, M., O'Neal, L., McLeod, L., Delacqua, G., Delacqua, F., Kirby, J., Duda, S. N., & REDCap Consortium (2019). The REDCap consortium: Building an international community of software platform partners. Journal of biomedical informatics, 95, 103208. [https://doi.org/10.1016/j.jbi.2019.103208](https://doi.org/10.1016/j.jbi.2019.103208)  
[^nextflow]: Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. [doi:10.1038/nbt.3820](doi:10.1038/nbt.3820)  
[^conda]: Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from [https://docs.anaconda.com/](https://docs.anaconda.com/)  
[^openjdk_java]: [https://openjdk.java.net/](https://openjdk.java.net/)
