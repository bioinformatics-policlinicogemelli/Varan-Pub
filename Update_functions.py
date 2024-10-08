import os
import pandas as pd
import re
import numpy as np
from loguru import logger
import shutil

def update_clinical_samples(oldfile_path, newfile_path, output_folder):
    """ 
    Update clinical sample information inside data_clinical_sample.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the sample info founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'data_clinical_sample.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original data_clinical_sample.
        newfile_path (str): Path to the new data_clinical_sample.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_clinical_samples('old_data_clinical_sample.txt', 'new_data_clinical_sample.txt', 'output_folder/')
    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype = str, index_col=0)
    new = pd.read_csv(newfile_path, sep="\t", dtype = str, index_col=0) #, skiprows=4) #, names=old.columns, skiprows=5)

    only_old_col = list(set(old.columns) - set(new.columns))
    # TODO scrivere queste colonne nel summary indicando che sono state rimosse perche presenti nell'old e assenti nel new

    common_samples = old.index.intersection(new.index)

    old_updated = old.drop(only_old_col, axis=1)

    old_updated = old.drop(index=common_samples)
    updated = pd.concat([new, old_updated], axis=0)
    updated.to_csv(os.path.join(output_folder, "data_clinical_sample.txt"), index=True, sep="\t", index_label=old.index.name)


def update_clinical_patient(oldfile_path, newfile_path, output_folder):
    """
    Update clinical patient informations inside data_clinical_patient.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the patients info founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'data_clinical_patient.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original data_clinical_patient.
        newfile_path (str): Path to the new data_clinical_patient.
        output_folder (str): Path to the output folder where the updated file will be saved.

    Example:
      >>>  update_clinical_patient('data_clinical_patient.txt', 'new_data_clinical_patient.txt', 'output_folder/')
    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype = str)
    old.set_index(old.columns[0], inplace=True)
    new = pd.read_csv(newfile_path, sep="\t", dtype = str) #, skiprows=4)
    new.set_index(new.columns[0], inplace=True)

    only_old_col = set(old.columns) - set(new.columns)
    # TODO scrivere queste colonne nel summary indicando che sono state rimosse perche presenti nell'old e assenti nel new

    common_patients = old.index.intersection(new.index)

    old_updated = old.drop(list(only_old_col), axis=1)
    old_updated = old_updated.drop(index=common_patients)

    updated = pd.concat([new, old_updated], axis=0).fillna(value=np.nan)
    updated.to_csv(os.path.join(output_folder, "data_clinical_patient.txt"), index=True, sep="\t", index_label=old.index.name)
    

def update_cna_hg19(oldfile_path, newfile_path, output_folder):
    """
    Update sample inside a copy number alteration (CNA) data file in hg19 genome format.
    This function reads the original tab separated CNA data from the given 'oldfile_path',
    insert new rows with the sample CNA data founded inside the new file from the given 'newfile_path' 
    and save the updated file named 'data_cna_hg19.seg' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original CNA data file.
        newfile_path (str): Path to the new CNA data file.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_cna_hg19('old_data_cna_hg19.seg', 'new_data_cna_hg19.seg', 'output_folder/')
    """
    old = pd.read_csv(oldfile_path, sep="\t")
    new = pd.read_csv(newfile_path, sep="\t")

    updated = pd.concat([old, new])
    updated = updated.drop_duplicates(subset=["ID","chrom","loc.start","loc.end"], keep='last')
    updated.to_csv(os.path.join(output_folder, "data_cna_hg19.seg"), index=False, sep="\t")
        
    

def update_cna(oldfile_path, newfile_path, output_folder):
    """
    Update sample inside a copy number alteration (CNA) data file.
    This function reads the original tab separated CNA data from the given 'oldfile_path',
    insert new rows with the sample CNA data founded inside the new file from the given 'newfile_path' 
    and save the updated file named 'data_cna.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original CNA data file.
        newfile_path (str): Path to the new CNA data file.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_cna('old_data_cna.txt', 'new_data_cna.txt', 'output_folder/')
    """

    old = pd.read_csv(oldfile_path, sep="\t", index_col=0)
    new = pd.read_csv(newfile_path, sep="\t", index_col=0)
    
    sample_old = old.columns 
    sample_new = new.columns
    to_remove = sample_old.intersection(sample_new)

    updated = pd.concat([old.drop(columns=to_remove), new], axis=1).replace(np.nan, 0).astype(int)
    updated.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")
        
    
def update_mutations(oldfile_path, newfile_path, output_folder):
    """
    Update samples' mutation data inside data_mutations_extended.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the samples' mutation data founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'data_mutations_extended.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original data_mutations_extended.
        newfile_path (str): Path to the new data_mutations_extended.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_mutations('old_data_mutations_extended.txt', 'new_data_mutations_extended.txt', 'output_folder/')
    """
    old = pd.read_csv(oldfile_path, sep="\t")
    new = pd.read_csv(newfile_path, sep="\t")

    updated = pd.concat([old, new])
    updated.drop_duplicates(keep='last', inplace=True)
    updated.to_csv(os.path.join(output_folder, "data_mutations_extended.txt"), index=False, sep="\t")
    
def update_sv(oldfile_path, newfile_path, output_folder):
    """
    Update samples' structural variation (SV) data inside data_sv.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the samples' SV data founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'data_sv.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original data_sv.
        newfile_path (str): Path to the new data_sv.
        output_folder (str): Path to the output folder where the updated file will be saved.

    Example:
      >>>  update_sv('old_data_sv.txt', 'new_data_sv.txt', 'output_folder/')
    """
    with open(oldfile_path,"r") as old_file:
        with open(newfile_path,"r") as new_file:
            with open(os.path.join(output_folder, "data_sv.txt"), "w") as of:

                for line in old_file:
                    list_split = line.split("\t\t")
                    list_strip = [elem.strip() for elem in list_split]
                    new_row = "\t".join(list_strip) + '\n'
                    of.write(new_row)
                for linenew in old_file:
                    if linenew.startswith("Sample_ID"):
                        new_row = linenew.replace("Sample_ID", "Sample_Id")
                        of.write(new_row)
                    if not linenew.startswith("Sample") :
                        list_split = linenew.split("\t\t")
                        list_strip = [elem.strip() for elem in list_split]
                        new_row = "\t".join(list_strip) + '\n'
                        of.write(new_row)
    data_sv = pd.read_csv(os.path.join(output_folder, "data_sv.txt"), sep="\t")
    data_sv = data_sv.drop_duplicates(subset=["Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol", "Normal_Paired_End_Read_Count"], keep='last')
    data_sv.to_csv(os.path.join(output_folder,"data_sv.txt"), sep="\t")
    

def update_caselist_cna(oldfile_path, newfile_path, output_folder):
    """
    Update cases' CNA data inside cases_cna.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the cases' CNA data founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'cases_cna.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original cases_cna.
        newfile_path (str): Path to the new cases_cna.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_caselist_cna('old_cases_cna.txt', 'new_cases_cna.txt', 'output_folder/')
    """
    with open(oldfile_path, "r") as old:
        with open(newfile_path, "r") as new:
            for line in new:
                line = line.strip()
                if line.startswith("case_list_ids"):
                    new_samples = line.split(":")[1]
                    len_new_sample = len(new_samples.split("\t"))
             
        with open(os.path.join(output_folder,"cases_cna.txt"),"w") as updated:           
            for line in old:
                if line.startswith("case_list_description"):
                    n_old_samples = re.findall(r'\d+', line)[0]
                    line = line.replace(n_old_samples, str(int(n_old_samples) + len_new_sample))
                if line.startswith("case_list_ids"):
                    new_samples_filtered = [sample for sample in new_samples.split("\t") if sample not in line]
                    line = "\t".join([line, "\t".join(new_samples_filtered)])
                updated.write(line)
            

def update_caselist_sequenced(oldfile_path, newfile_path, output_folder):
    """
    Update cases' sequenced data inside cases_sequenced.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the cases' sequenced data founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'cases_sequenced.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original cases_sequenced.
        newfile_path (str): Path to the new cases_sequenced.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_caselist_sequenced('old_cases_sequenced.txt', 'new_cases_sequenced.txt', 'output_folder/')
    """
    with open(oldfile_path, "r") as old:
        with open(newfile_path, "r") as new:
            for line in new:
                line = line.strip()
                if line.startswith("case_list_ids"):
                    new_samples = line.split(":")[1]
                    len_new_sample = len(new_samples.split("\t"))
             
        with open(os.path.join(output_folder, "cases_sequenced.txt"), "w") as updated:           
            for line in old:
                if line.startswith("case_list_description"):
                    n_old_samples = re.findall(r'\d+', line)[0]
                    line = line.replace(n_old_samples, str(int(n_old_samples) + len_new_sample))
                if line.startswith("case_list_ids"):
                    new_samples_filtered = [sample for sample in new_samples.split("\t") if sample not in line]
                    line = "\t".join([line, "\t".join(new_samples_filtered)])
                updated.write(line)
     
def update_caselist_sv(oldfile_path, newfile_path, output_folder):
    """
    Update cases' structural variation (SV) data inside cases_sv.txt file.
    This function reads the original tab separated version txt file from the given 'oldfile_path',
    insert new rows with the cases' SV data founded inside the new txt file from the given 'newfile_path' 
    and save the updated file named 'cases_sv.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original cases_sv.
        newfile_path (str): Path to the new cases_sv.
        output_folder (str): Path to the output folder where the updated file will be saved.
    
    Example:
      >>>  update_caselist_sv('old_cases_sv.txt', 'new_cases_sv.txt', 'output_folder/')
    """
    with open(oldfile_path, "r") as old:
        with open(newfile_path, "r") as new:
            for line in new:
                line = line.strip()
                if line.startswith("case_list_ids"):
                    new_samples = line.split(":")[1]
                    len_new_sample = len(new_samples.split("\t"))
             
        with open(os.path.join(output_folder, "cases_sv.txt"), "w") as updated:           
            for line in old:
                if line.startswith("case_list_description"):
                    n_old_samples = re.findall(r'\d+', line)[0]
                    line = line.replace(n_old_samples, str(int(n_old_samples) + len_new_sample))
                if line.startswith("case_list_ids"):
                    new_samples_filtered = [sample for sample in new_samples.split("\t") if sample not in line]
                    line = "\t".join([line, "\t".join(new_samples_filtered)])
                updated.write(line)


def check_files(oldpath, newpath, output, file_name):
    o_data = os.path.join(oldpath, file_name)
    n_data = os.path.join(newpath, file_name)
    if os.path.exists(o_data) and os.path.exists(n_data):
        if file_name == "data_clinical_sample.txt":
            update_clinical_samples(o_data,n_data,output)
        elif file_name == "data_clinical_patient.txt":
            update_clinical_patient(o_data,n_data,output)
        elif file_name == "data_cna_hg19.seg":
            update_cna_hg19(o_data,n_data,output)
        elif file_name == "data_cna.txt":
            update_cna(o_data,n_data,output)
        elif file_name == "data_mutations_extended.txt":
            update_mutations(o_data,n_data,output)
        elif file_name == "data_sv.txt":
            update_sv(o_data,n_data,output)
    elif os.path.exists(o_data) and not os.path.exists(n_data):
        logger.warning(f"{file_name} was not found in path 2. The file is being copied from path 1.")
        shutil.copy(o_data, output)
    elif not os.path.exists(o_data) and os.path.exists(n_data):
        logger.warning(f"{file_name} was not found in path 1. The file is being copied from path 2.")  
        shutil.copy(n_data, output)
    else:
        logger.warning(f"{file_name} not found in current folders. Skipping")


def check_files_cases(oldpath, newpath, output_caseslists, file_name):
    o_data = os.path.join(oldpath,"case_lists",file_name)
    n_data = os.path.join(newpath,"case_lists",file_name)
    if os.path.exists(o_data) and os.path.exists(n_data):
        if file_name == "cases_cna.txt":
            update_caselist_cna(o_data,n_data,output_caseslists)
        elif file_name == "cases_sequenced.txt":
            update_caselist_sequenced(o_data,n_data,output_caseslists)
        elif file_name == "cases_sv.txt":
            update_caselist_sv(o_data,n_data,update_caselist_sv)
    elif os.path.exists(o_data) and not os.path.exists(n_data):
        logger.warning(f"{file_name} was not found in path 2 case_lists folder. The file is being copied from path 1.")  
        shutil.copy(o_data, output_caseslists)
    elif not os.path.exists(o_data) and os.path.exists(n_data):
        logger.warning(f"{file_name} was not found in path 1 case_lists folder. The file is being copied from path 2.") 
        shutil.copy(n_data, output_caseslists)
    else:
        logger.warning(f"{file_name} not found in 'case_lists' folders. Skipping")
    