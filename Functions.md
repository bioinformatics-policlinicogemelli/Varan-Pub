
# In-depth functions

### Block One: Create Study Folder ex-Novo

<p align="justify">Below are reported and explained the functions called by the main `varan.py`. These functions can also be called separetely by giving their specific input.

* **walk.py**<br><p align="justify">The `walk.py`  function is a Python script that processes Variant Call Format (VCF) files containing genetic information about Single Nucleotide Polymorphism (SNP) and Copy Number Variation (CNV) variants. Once the `walk` is executed, you will have an output containing the Minor Allele Frequencies (MAF).


* **filter_clinvar.py**<br><p align="justify">The `filter_clinvar.py` function is a Python script to filter mutations according to ClinVar annotation and/or mutation consequence and according to VAF and gnomAD frequency.

	-  **VUS option**:

		Through this filter all mutations that **aren't annotated as *benign* or *likely benign*** will be maintained.

  

	-  **NO_VUS option**:

		The first step is the elimination of all mutations annotated as *benign* or *likely benign*.

		Next, the script retaines variants annotated in **Clinvar** such as

		-  *risk_factor*

		-  *pathogenic*

		-  *likely_pathogenic*

		-  *drug_response*.

		
		or mutations that have terms in the Consequence column such as

  

		-  *splice_region_variant*,

		-  *frameshift_variant*,

		-  *splice_donor_variant*,

		-  *stop_gained*,

		-  *splice_acceptor_variant*,

		-  *start_lost*,

		-  *inframe_insertion*,

		-  *inframe_deletion*

	Based on frequency, are maintened mutations that have a **VAF** (Variant Allele Frequency) > 0.05 and mutations that have a **gnomAD** frequency < 0.98

  
* **Concatenate.py**<p align="justify"> The `concatenate.py` is a Python script that can be used to merge the MAF files of single 	patients, into a unique file<br> 

  
* **Make_meta_and_cases.py**<br><p align="justify"> The `Make_meta_and_cases.py` function is a Python script that creates all the required metafiles
	* meta_study.txt
	* meta_clinical_patient.txt
	* meta_clinical_sample.txt
	* data_mutations_extended.txt
	* meta_sv.txt
	* meta_cna.txt
	* meta_cna_hg19_seg.txt

	It also create the case list files:
	* cases_sequenced.txt
	* cases_cna.txt
	* cases_sv.txt		


* **ValidateFolder.py**<br><p align="justify"> The `ValidateFolder.py` checks the contents of a study folder against a set of required files for different categories (e.g., Patient, Study, CNA, Fusion, SNV) that are necessary for uploading data to cBioPortal. It logs any missing files and provides a success message if all required files are present.<br>

#
### Block Two: Modify Existing Study Folder

 * **Update_script.py**<br><p align="justify">Th function `Update_script.py` covers the update samples process. In particular, after giving an existing study folder and a newer version of this study folder containing novel or updated information, *Update_script.py* concatenates the newer info and update the oldest with the new ones. These modifications are applied to all of the txt files contained inside the study folder<p align="justify">The output will be a new study folder saved as the study folder original name tagged with a version number to  keep track of changes. This new folder will be an updated version of the original one.<br>
 <br>


* **Delete_script.py**<br> <p align="justify">Th covers the delete samples process. In particular, after giving an existing study folder and a list of samples to remove, *Delete_script.py* deletes the info regarding the samples of interest. These modifications are applied to all of the txt files contained inside the study folder.<br>

<br>

* **ExtractSamples_script.py**<br><p align="justify">This function covers the extract samples process. In particular, after giving an existing study folder and a list of samples to extract, *ExtractSample_script.py* extracts the info regarding the samples of interest. These modifications are applied to all of the txt files contained inside the study folder.<br><br>