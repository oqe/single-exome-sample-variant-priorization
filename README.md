# Single Exome Sample Variant Priorization

This is a nextflow[^nextflow] workflow for extracting (human) single sample from a larger exome sequencing multisample/cohort VCF and annote and prioritize exome variants of that sample. "Single sample" as in the context of the workflow stages per sample. You may execute the workflow for multiple samples in parallel.

**Here is a short description of the workflow stages:**  
1. Generate sample and date specifc output directories
2. Single sample vcf file is extracted from the cohort vcf with GATK SelectVariants (-sn) tool [^gatk].
3. Extracted vcf is decomposed and normalized by vt tool[^vt].
4. Decomposed and normalized vcf is compressed (.gz) using bgzip and indexed (.tbi) using tabix[^htslib].
5. Generate sample specific .yaml files[^python]<sup>,</sup>[^pyyaml] for Exomiser from input table information[^exomiser]. Annotation and priorization for the variants of the vcf with sample specific HPO terms is done with Exomiser [^hpo].
6. Generate sample specific .yaml files for LIRICAL from input table information and "Phenotype-driven priorization of candidate diseases and genes..." is done with LIRICAL (: LIkelihood Ratio Interpretation of Clinical AbnormaLities) [^lirical].  

## Why?

Streamlining individual patient genetic analysis from large cohort of (exome) sequenced patients.

Nextflow for fast deployment with good support for SLURM and cloud service providers[^nextflow]. Conda for easily staying up-to-date, easy management and deployment of available software (versions)[^conda]. 

Exomiser and LIRICAL for applying phenotype information with genetic information for variant priorization and helping in diagnosis effort with human phenotype ontology terms.

HPO term extraction per samplename is separated as an external R Script[^r]<sup>,</sup>[^r-redcap]. R script is for sample(name) exploration/checking from REDCap database as well as extracting HPO terms and generating sample input table required as input in this workflow. Expectation is that in your use case your institution/or your affiliates is/are hosting REDCap and have HPO terms saved there[^redcap1]<sup>,</sup>[^redcap2]. Alternatively you can input the HPO terms per patient manually into *input_fofn* table.

## Software/data requirements
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
		- with Exomiser data download and configured to work
		- application.properties configured
	+ LIRICAL 
		- with LIRICAL data downloaded and configured to work

+ Reference genome files
	- .fasta, .fai, .dict
+ Input
	- multi-patient vcf (and index)
	- samplenames matching to vcf
	- hpo-terms per patient OR REDCap API access and REDCap specific samplenames (matching to vcf samplenames)

## Input

### params.yaml

Edit **params.yaml** and for your inputs:
<details><summary>DETAILS</summary>
<p>
- input_fofn file path:  
	- file should be sample input table where (column) values should be tab separated (.tsv), with following variables per sample per line  
		1. unique id / sample name in cohort VCF file, with or without cohort/sub-cohort designation  
		2. samplename in REDCap/database where HPO terms were extracted  
		3. list of HPO terms per sample ['HPO:XXX', 'HPO:XXX']  
- cohort_prefix  
	- if your samplename has cohort prefix in the vcf file you can add it here  
	- you may also leave this empty but then the full samplename needs to be input in the first column of input_fofn file  
- cohort_vcf  
	- path to cohort vcf file  
- cohort_vcf_index  
	- path to cohort vcf index file  
- conda_path  
	- path to your conda enviroment  
	- you can import the conda enviroment from /templates/conda_env_variant-tools.yml  
- vt_path  
	- path to vt, default provided as per conda env: vt  
- vt_genome   
	- vt generated reference genome  
- bgzip_path  
	- path to bgzip, default is provided as per conda env: bgzip  
- tabix_path  
	- path to bgzip, default is provided as per conda env: tabix  
- exomiser_path  
	- path to exomiser jar file  
- exomiser_config  
	- path to application.properties file of exomiser  
- exomiser_yaml_template  
	- path to exomiser yaml analysis file  
	- template is provided in /templates/exomiser_13.0.0_hg38_template.yaml, edit this file with your values  
- exomiser_outputprefix (string output prefix)  
	- exomiser result file output prefix string  
- lirical_path  
	- path to LIRICAL.jar file  
- lirical_yaml_template  
	- path to lirical yaml analysis file  
	- template is provided in /templates/lirical_1.3.4_template.yaml, edit this file with your values and fill in exomiser data directory path in the template  
- lirical_outputprefix  
	- lirical result file output prefix string  
- outdir  
	- main output directory path  
- fasta  
	- reference genome fasta file path  
	- !Note .fai and .dict files are also required  
</p>
</details>  
  
  
**input_fofn** sample input file's format is tab separated values (.tsv) text file. First field is samplename in cohort vcf file, separated by tab, then samplename in REDCap database and tab separted by HPO terms field with (possibly multiple) HPO terms inside, where each HPO term is quoted with single quotes and everything is enclosed with angle brackets for example: **['HPO:XXXXX', 'HPO:ZZZZZ', 'HPO:YYYYYY']**. 

| samplename_cohort_vcf | samplename_redcap | hpo_ids |
| --- | --- | --- |
| X1 | XRX200 | ['HPO:XXXX', 'HPO:YYYY', 'HPO:ZZZZZ'] |

As the samplenames can be different in a vcf file and in redcap(or other equivalen database) there is a column for REDCap specific samplename and cohort-VCF file specific samplename. Cohort-VFC samplename can be full samplename and then the cohort prefix variable can be left empty (""). Cohort-VFC samplename can also be unique part of samplename if samplenames have cohort specific prefix name. In this case you need to fill in cohort prefix variable that so that the cohort_prefix + samplename_cohort_vcf create a valid samplename that is in the vcf file.

Exomiser and LIRICAL outputprefix variables (string) are for filename output prefixes. These output prefixes are later added in the whole output file path to yaml files when sample specific yaml files are generated.

+ Reference genome files (Note! be sure to use same version as Exomiser/LIRICAL)

You can download reference files...


### nextflow.config

Edit your **nextflow.config** according to how you plan to execute the workflow, for example local execution, SLURM etc.


## Output

In the beginning sample specific output directories are created based on given main output directory parameter, samplename and date (yyyy-MM-dd). Each stage will save specific files to this sample specific directory.

	<input_main_output_dir> / <sample_name> / <date_in_yyyy-MM-dd>

After the workflow has completed the following files should be in each date specific folder inside sample name folder:  
- exomiser yaml analysis file  
- lirical yaml analysis file  
- (vt decomposed and normalized) bgzip compressed .vcf.gz and tabix indexed .vcf.gz.tbi  
- exomiser output files (many files)  
- lirical output file  

## TO-DO

- [ ] Update conda enviroment due to Log4 vulnerability
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

- [ ] Change/add sample specific parametres to **input_fofn** table...  
	- [ ] output directory
	- [ ] cohort-vcf  
	*enables executing sample variant priorization for samples from different cohort or cohort-vcf versions in with the same input_fofn table file*
- [ ] Document zero-to-finish example of this workflow to a separate file  
- [ ] Finalize extra R script (REDCap HPO-terms extraction)  
	- [ ] Write up documentation  
- [ ] Add original .wdl workflow in extras/wdl folder
- [ ] Add VEP annotation?
- [ ] Clean-up, modulize code?  
- [ ] Add checking if sample specific .vcf (.gz and .tbi) files are found in the sample output directory  
	- [ ] Add switch to utilize found .vcf file and skip first step (GATK_SELECTVARIANTS) to save time  
- [ ] Inbedding sample specific .bam/.cram file to exomiser output .html with IGV.js  
- [ ] Docker/Apptainer formerly known as singularity) container for programs ???  
- [ ] create simple graphic/graph about the workflow  

## Tool specific manual/download/homepage links

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
