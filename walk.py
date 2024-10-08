#####################################
# NAME: walk.py
# Date: 10/01/2023
version = "1.0"
# ===================================

from operator import index
import os
import sys
import ast
import subprocess
import vcf2tab_cnv
import vcf_filter
import tsv
import pandas as pd
from configparser import ConfigParser
import shutil
from loguru import logger 
import random
import string
import numpy as np
from filter_clinvar import filter_OncoKB
from versioning import get_newest_version, get_version_list
import pandas as pd
from datetime import datetime

config = ConfigParser()
configFile = config.read("conf.ini")

VCF2MAF = config.get('Paths', 'VCF2MAF')
REF_FASTA = config.get('Paths', 'REF_FASTA')
VEP_PATH = config.get('Paths', 'VEP_PATH')
VEP_DATA = config.get('Paths', 'VEP_DATA')
CLINV = config.get('Paths', 'CLINV')
CNA = ast.literal_eval(config.get('Cna', 'HEADER_CNV'))
PLOIDY = int(config.get('Cna', 'PLOIDY'))

output_filtered = "snv_filtered"
tmp = "scratch"

def create_random_name_folder():
    nome_cartella = ''.join(random.choices(string.ascii_lowercase + string.digits, k=10))
    temporary = os.path.join(tmp, nome_cartella)
    try:
        os.mkdir(temporary)
    except FileNotFoundError:
        logger.critical(f"Scratch folder '{tmp}' not found!")
        raise(FileNotFoundError("Error in create_random_name_folder: exiting from walk script!"))
    except Exception:
        logger.critical("Something went wrong while creating the vep tmp folder")
        raise(Exception("Error in create_random_name_folder: exiting from walk script!"))
    return(temporary)


def clear_scratch():
    for root, dirs, files in os.walk(tmp):
        for dir in dirs:
            shutil.rmtree(os.path.join(root,dir))
        

def clear_temp(folder):
    shutil.rmtree(os.path.join(folder, "temp"))


def get_cnv_from_folder(inputFolderCNV):
    files = os.listdir(inputFolderCNV)
    cnv_vcf_files = [file for file in files if file.endswith("vcf")]
    return cnv_vcf_files


def get_sampleID_from_cnv(cnv_vcf):
    if "_CopyNumberVariants.vcf" in cnv_vcf:
        sample=cnv_vcf.replace("_CopyNumberVariants.vcf", ".bam")
    else:
        sample=cnv_vcf.replace("vcf", "bam")
    return sample


def reshape_cna(input, cna_df_path, cancer, output_dir):
   
    if not os.path.isfile(input):
        input_file = pd.read_csv(os.path.join(input, "sample.tsv"), sep="\t")    
    else:
        input_file = pd.read_csv(input, sep="\t")

    cna_df = pd.read_csv(cna_df_path, sep="\t")
    
    cna_df.rename({"ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"}, inplace=True, axis=1)
    input_file.rename({"SampleID":"Tumor_Sample_Barcode"}, inplace=True, axis=1)

    if not "ONCOTREE_CODE" in input_file.columns:
        input_file["ONCOTREE_CODE"] = cancer
  
    input_file["Tumor_Sample_Barcode"] = input_file["Tumor_Sample_Barcode"] + ".cnv.bam"
    
    annotate = pd.merge(cna_df[["Tumor_Sample_Barcode", "Hugo_Symbol", "discrete", "Copy_Number_Alteration"]], \
                        input_file[["Tumor_Sample_Barcode", "ONCOTREE_CODE"]], on="Tumor_Sample_Barcode")
   
    
    annotate.to_csv(os.path.join(output_dir, "temp_cna.txt"), index=False, sep="\t")
    
    return os.path.join(output_dir, "temp_cna.txt")


def annotate_cna(path_cna, output_folder):
    
    out = path_cna.replace(".txt", "2.txt")
    os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {path_cna}\
                        -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")
                
    cna = pd.read_csv(out, sep="\t", dtype={"Copy_Number_Alteration":int})
    cna = cna[cna["ONCOGENIC"].isin(["Oncogenic", "Likely Oncogenic"])]
        
    data_cna = cna.pivot_table(index="Hugo_Symbol", columns="Tumor_Sample_Barcode", values="Copy_Number_Alteration", fill_value=0)

    data_cna.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")


def cnv_type_from_folder(input, cnv_vcf_files, output_folder, oncokb, cancer, multiple):
    
    c = 0
    sID_path = dict()
    for case_folder in cnv_vcf_files:
        if os.path.exists('data_cna_hg19.seg'):
            MODE = 'a'
        else:
            MODE = 'w'
        try:
            cnv_vcf = case_folder 
            sampleID = get_sampleID_from_cnv(case_folder)

            if sampleID in sID_path:
                dup = open('sampleID_dup'.log, 'w')
                dup.write(sampleID + '\t' + 'cnv_vcf')
                dup.close()
            else:
                if multiple:
                    sID_path[sampleID] = os.path.join(os.path.join(input, "CNV", "single_sample_vcf"), cnv_vcf)
                else:
                    sID_path[sampleID] = os.path.join(os.path.join(input, "CNV"), cnv_vcf)
                
                vcf2tab_cnv.vcf_to_table(sID_path[sampleID], os.path.join(output_folder, 'data_cna_hg19.seg'), sampleID, MODE)
                vcf2tab_cnv.vcf_to_table_fc(sID_path[sampleID], os.path.join(output_folder, 'data_cna_hg19.seg.fc.txt'), sampleID, MODE)
               
        except Exception:
            log_noparsed = open('noParsed_cnv.log', 'a')
            log_noparsed.write('[WARNING] ' + case_folder + '\n')
            log_noparsed.close()
        
        c = c + 1
    logger.info("Writing data_cna_hg19.seg succefully completed!")
    logger.info("Writing data_cna_hg19.seg.fc.txt succefully completed!")
    
    ############################
    ### MANAGE DISCRETE TABLE ##
    ############################
    
    logger.info("Starting cna evaluation (this step could take a while)...")
    df_table = pd.read_csv(os.path.join(output_folder, 'data_cna_hg19.seg.fc.txt'), sep="\t", header=0)
    result = table_to_dict(df_table)
        
    df_table.rename(columns={"discrete":"Copy_Number_Alteration", "ID":"Tumor_Sample_Barcode", "gene":"Hugo_Symbol"}, inplace=True)
    df_table_filt = df_table[df_table["Copy_Number_Alteration"].isin([-2,2])]
        
    if oncokb:
        
        ################################# OPZIONE 1 #########################################################
        # if not os.path.isfile(input):
        #     tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
        #     input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
        
        # else:
        #     input_file=pd.read_csv(input,sep="\t")
        
        # for index, row in df_table.iterrows():
        #     tc= int(input_file[input_file["SampleID"]+".cnv.bam"==row["ID"]]["TC"])
        #     df_table.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc
            
        
        # df_table_filtered=df_table#[(df_table["discrete"]>=3)|(df_table["discrete"]<=0.8)] 
        
        #df_table_filtered["Copy_Number_Alteration"]=0
        # df_table_filtered.loc[(df_table_filtered["discrete"]>3)&(df_table_filtered["discrete"]<5), "Copy_Number_Alteration"]=1
        # df_table_filtered.loc[df_table_filtered["discrete"]>5, "Copy_Number_Alteration"]=2
        # df_table_filtered.loc[(df_table_filtered["discrete"]>0)&(df_table_filtered["discrete"]<0.8), "Copy_Number_Alteration"]=-1
        # df_table_filtered.loc[df_table_filtered["discrete"]<=0, "Copy_Number_Alteration"]=-2
        
       
        #df_table_filtered["Copy_Number_Alteration"]=df_table_filtered["discrete"]
        #df_table_filtered.to_csv(os.path.join(output_folder,"data_cna_hg19.seg.filtered.txt"),sep="\t",index=False)
        #temp_cna=reshape_cna(input,os.path.join(output_folder,"data_cna_hg19.seg.filtered.txt"),cancer,output_folder)
        #annotate_cna(temp_cna,output_folder)
   
        ################################# OPZIONE 2 #########################################################
   
        # rename discrete 
        
        #df_table.rename(columns={"discrete":"Copy_Number_Alteration","ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"},inplace=True)
        #df_table_filt=df_table[df_table["Copy_Number_Alteration"].isin([-2,2])]
    
        if not os.path.isfile(input):
            input_file = pd.read_csv(os.path.join(input, "sample.tsv"), sep="\t")
        else:
            input_file = pd.read_csv(input, sep="\t")
        
        #TODO inserire filtro per TC (mettere specifiche su valori di TC da considerare)
        
        #check for nan TC values
        try:
            input_file["TC"]
        except:
            input_file["TC"] = np.nan
        
        if len(input_file[input_file["TC"].isna()])>0:
            nan_sbj = input_file[input_file["TC"].isna()]
            nan_sbj = list(nan_sbj["SAMPLE_ID"])
            logger.warning(f"Some subject have Nan TC in tsv input file: {nan_sbj}!")
        
        if not "ONCOTREE_CODE" in input_file.columns:
            input_file["ONCOTREE_CODE"] = cancer

        input_file["Tumor_Sample_Barcode"] = input_file["SAMPLE_ID"] #+ ".cnv.bam"
        
        annotate = pd.merge(df_table_filt[["Tumor_Sample_Barcode","Hugo_Symbol", "seg.mean", "Copy_Number_Alteration"]], \
                           input_file[["Tumor_Sample_Barcode","ONCOTREE_CODE","TC"]], on="Tumor_Sample_Barcode")
        temppath = os.path.join(output_folder, "temp_cna_toannotate.txt")
        annotate.to_csv(temppath, index=False, sep="\t")

        out=temppath.replace("toannotate.txt", "annotated.txt")
        os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {temppath}\
                        -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")
        
           
        cna = pd.read_csv(out,sep="\t", dtype={"Copy_Number_Alteration":int})
        cna = cna[cna["ONCOGENIC"].isin(["Oncogenic", "Likely Oncogenic"])]
        cna["ESCAT"]="Unmatched"
        df_table["ESCAT"]="Unmatched"
        
        logger.info("Analyzing cna sample(s)")
        for _, row in cna.iterrows():
            
            #logger.info("Analyzing cna sample " + row["Tumor_Sample_Barcode"])            
            try:
                tc = int(input_file[input_file["Tumor_Sample_Barcode"] == row["Tumor_Sample_Barcode"]]["TC"])
            except ValueError:
                continue
            except Exception:
                raise(Exception(f"Something went wrong while reading TC!"))

            #decomm se cutoff su copy number discreto
            #cna.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc 
            
            purity = tc / 100
            copy_nums = np.arange(6)
            c = 2 ** (np.log2((1 - purity) + purity * (copy_nums + .5) / PLOIDY ))
            
            # ESCAT classification
            escat_class(df_table, input_file, row)
                    
        #nuovi filtri (filtro cnvkit - conservativo)
        
        if not cna.empty:
            cna["Copy_Number_Alteration"]=0
            cna.loc[(cna["seg.mean"]<c[0]), "Copy_Number_Alteration"]=-2
            cna.loc[(cna["seg.mean"]>=c[0])&(cna["seg.mean"]<c[1]), "Copy_Number_Alteration"]=-1
            cna.loc[(cna["seg.mean"]>=c[1])&(cna["seg.mean"]<c[3]), "Copy_Number_Alteration"]=0
            cna.loc[(cna["seg.mean"]>=c[3])&(cna["seg.mean"]<c[5]), "Copy_Number_Alteration"]=1
            cna.loc[cna["seg.mean"]>=c[5], "Copy_Number_Alteration"]=2
        
        else:
            logger.warning("After filtering cna dataframe is empty!")
               
        #vecchi filtri (decomm se filtro copy number discreto)
        
        #decomm se cutoff su copy number discreto
        #cna.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc 
            
        # cna["Copy_Number_Alteration"]=0
        # cna.loc[(cna["discrete"]>3)&(cna["discrete"]<5), "Copy_Number_Alteration"]=1
        # cna.loc[cna["discrete"]>5, "Copy_Number_Alteration"]=2
        # cna.loc[(cna["discrete"]>0.25)&(cna["discrete"]<1), "Copy_Number_Alteration"]=-1
        # cna.loc[cna["discrete"]<=0.25, "Copy_Number_Alteration"]=-2

        df_table.to_csv(os.path.join(output_folder, "data_cna_hg19_escat.seg.fc.txt"), index=True, sep="\t")
        
        cna.to_csv(os.path.join(output_folder, "CNA_ANNOTATI_NDISCRETO.txt"), index=True, sep="\t")
        
        cna["Tumor_Sample_Barcode"] = cna["Tumor_Sample_Barcode"].str.replace(".cnv.bam", "")
        data_cna=cna.pivot_table(index="Hugo_Symbol", columns="Tumor_Sample_Barcode", values="Copy_Number_Alteration", fill_value=0)
        data_cna.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")
        
    else:
        
        df_table_filt["Tumor_Sample_Barcode"] = df_table_filt["Tumor_Sample_Barcode"].str.replace(".cnv.bam", "")
        
        data_cna = df_table_filt.pivot_table(index="Hugo_Symbol", columns="Tumor_Sample_Barcode", values="Copy_Number_Alteration", fill_value=0)
        data_cna.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")
    
    return sID_path

def escat_class(df_table, input_file, row):
    oncocode = input_file[input_file["Tumor_Sample_Barcode"] == row["Tumor_Sample_Barcode"]]["ONCOTREE_CODE"].values
    gene = df_table[(df_table["Tumor_Sample_Barcode"] == row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"] == row["Hugo_Symbol"])]

    if oncocode =="NSCLC":  
        if (row["Hugo_Symbol"]=="MET" and row["seg.mean"]>1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
        elif (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
            
    elif oncocode=="BREAST":
        if (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
            row["ESCAT"]="IA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"
                    
    elif oncocode=="BOWEL":
        if (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
                    
    elif oncocode=="PROSTATE":
        if (row["Hugo_Symbol"]=="BRCA1" and row["seg.mean"]<1):
            row["ESCAT"]="IA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"
        elif (row["Hugo_Symbol"]=="BRCA2" and row["seg.mean"]<1):
            row["ESCAT"]="IA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"
        elif (row["Hugo_Symbol"]=="PTEN" and row["seg.mean"]<1):
            row["ESCAT"]="IIA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIA"
        elif (row["Hugo_Symbol"]=="ATM" and row["seg.mean"]<1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
        elif (row["Hugo_Symbol"]=="PALB2" and row["seg.mean"]<1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"


# def table_to_dict(df):
#     result = {}
#     for _, row in df.iterrows():
#         row_values = (row['chrom'], row['loc.start'], row['loc.end'], row['num.mark'], row['seg.mean'], row['gene'], row['discrete'])
#         if row['ID'] not in result:
#             result[row['ID']] = []
#         result[row['ID']].append(row_values)
#     return result

def table_to_dict(df):
    result = {}
    for row in df.itertuples(index=False):
        row_values = (row.chrom, row._2, row._3, row._4, row._5, row.gene, row.discrete)
        if row.ID not in result:
            result[row.ID] = []
        result[row.ID].append(row_values)
    return result

def get_snv_from_folder(inputFolderSNV):
    files = os.listdir(inputFolderSNV)
    
    snv_vcf_files = [file for file in files if file.endswith("vcf")]
    return snv_vcf_files


def get_sampleID_from_snv(snv_vcf):
    if "MergedSmallVariants.genome.vcf" in snv_vcf:
        sample = snv_vcf.replace("_MergedSmallVariants.genome.vcf", ".bam")
    else:
        sample = snv_vcf.replace("vcf", "bam")
    return sample


def snv_type_from_folder(input, snv_vcf_files):
    c = 0
    sID_path = dict()
    for case_folder in snv_vcf_files:
        try:
            snv_vcf= case_folder
            sampleID = get_sampleID_from_snv(case_folder)
            if sampleID in sID_path:
                dup = open('sampleID_dup'.log, 'w')
                dup.write(sampleID + '\t' + 'snv_vcf')
                dup.close()
            else:
                sID_path[sampleID] = os.path.join(input, snv_vcf)
        except Exception:
            log_noparsed = open('noParsed_snv.log', 'a')
            log_noparsed.write('[WARNING]' + case_folder + '\n')
            log_noparsed.close()
        c = c + 1
    return sID_path


def vcf_filtering(sID_path, output_folder, output_filtered):
    sID_path_filtered = dict()
    if output_filtered.strip() == "":
        output_filtered = "snv_filtered"
    os.makedirs(os.path.join(output_folder, output_filtered), exist_ok=True)
    for k, v in sID_path.items():
        _, vcf_file = os.path.split(v)
        out_filt = os.path.join(output_folder, output_filtered) #TEST
        vcf_filtered = os.path.join(out_filt, vcf_file.replace('.vcf','') + '.FILTERED.vcf')
        #logger.info(f'[FILTERING] {v}')
        vcf_filter.main(v, vcf_filtered)
        #logger.info(f'filtered file {vcf_filtered} created!')
        sID_path_filtered[k] = vcf_filtered
    return sID_path_filtered


def vcf2maf_constructor(k, v, temporary, output_folder):

    CACHE=config.get('Paths', 'CACHE')
    #cmd = "grep $'\t'" + k.split(".")[0] + " " + v
    cmd = "vcf-query -l "+v
    try:
        tum_id = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
        #tum_id = [i for i in tum_id.split("\t") if k.split(".")[0] in i][0]
    except Exception:
        tum_id = ""
    
    if VCF2MAF == "" or REF_FASTA == "" or VEP_PATH == "" or VEP_DATA == "":
        logger.critical("[Paths] section in conf.ini is not correctly compiled. Please check again!")
        raise Exception ("Input error")
    
    cl = ['perl']
    cl.append(VCF2MAF)
    cl.append('--input-vcf')
    cl.append(v)
    _, file_vcf = os.path.split(v)

    out_file = os.path.join(output_folder, os.path.join("maf", file_vcf + '.maf'))
    if not CLINV.strip() == "":
        cl.append('--vep-custom')
        cl.append(CLINV)
    else:
        logger.warning(f"CLINV section in [Paths] in conf.ini is not compiled. This step will be skipped")
    cl.append('--output-maf')
    cl.append(out_file)
    cl.append('--ref-fasta')
    cl.append(REF_FASTA)
    cl.append('--tmp-dir')
    cl.append(temporary)
    cl.append('--retain-fmt')
    cl.append('GT,GQ,AD,DP,VF,AF')
    cl.append('--vep-path')
    cl.append(VEP_PATH)
    cl.append('--vep-data')
    cl.append(VEP_DATA)
    cl.append('--tumor-id') 
    cl.append(tum_id)
    cl.append('--cache-version')
    cl.append(CACHE)
    
    return cl

def run_vcf2maf(cl):
    logger.info('Starting vcf2maf conversion...')
    logger.info(f'args={cl}')
    sout = subprocess.run(cl, capture_output=True)
 
    if sout.stderr != None:
        if 'ERROR' not in sout.stderr.decode('ascii'):
            logger.warning(sout.stderr.decode('ascii').replace('ERROR: ',''))
        else:
            logger.error(sout.stderr.decode('ascii').replace('ERROR: ',''))


def create_folder(output_folder, overwrite_output, resume):
    
    output_list = get_version_list(output_folder)
    output_list=list(map(lambda x: os.path.join(os.path.dirname(output_folder),x),output_list))
    
    if output_list != [] and os.path.exists(output_list[-1]):
        output = output_list[-1]
        logger.warning(f"It seems that a version of the folder '{output_folder}' already exists.")
        if overwrite_output:
            logger.info(f"Overwrite option set. Start removing folder")
            shutil.rmtree(output)
        elif resume:
            _,current_version = get_newest_version(output_folder)
            return output_folder + current_version

    if output_list == []:
        version = "_v1"
        output_folder_version = output_folder + version  
    else:
        output_folder_version,_ = get_newest_version(output_folder)

    logger.info(f"Creating the output folder '{output_folder_version}' in {os.path.dirname(output_folder)}...")
    
    os.makedirs(output_folder_version, exist_ok=True)
    maf_path = os.path.join(output_folder_version, 'maf')
    os.makedirs(maf_path, exist_ok=True)
    logger.info(f"The folder '{output_folder_version}' was correctly created!")
    
    return output_folder_version    


def get_table_from_folder(tsvpath):
    table_dict = dict()
    file = pd.read_csv(tsvpath, sep="\t", index_col=False, dtype=str)
    for _, row in file.iterrows():
        sampleID = str(row["SAMPLE_ID"])
        if ".bam" in sampleID:
           sampleID = sampleID.replace(".bam", "")
        if sampleID not in table_dict.keys():
            table_dict[sampleID] = [str(row["PATIENT_ID"])]  
    return table_dict


def flatten(nested_list):
    flat_list = []
    for sublist in nested_list:
        for item in sublist:
            flat_list.append(item)
    return flat_list


def check_multiple_file(input_file, multiple):
    
    conf_snv = config.get('Multiple', 'SNV')
    conf_cnv = config.get('Multiple', 'CNV')
    file_paths = pd.read_csv(input_file, sep="\t", header=0)
    snv_file_path = file_paths['snv_path'].isnull().all() # True mancano tutti i valori
    cnv_file_path = file_paths['cnv_path'].isnull().all() # False se almeno uno è pieno
    
    #TODO ridiscutere CASE 4
    #CASE 1: multiple = True; neither snv or cnv are filled in sample.tsv; multiple snv or cnv are filled in conf.ini
    if multiple and not (snv_file_path or cnv_file_path) and not (conf_snv == "" or conf_cnv == ""):
        logger.critical("-m was selected and Muliple section in conf.ini was filled but the file doesn't looks like a multiVCF.")
        raise Exception ("Input error")
    #CASE 2: multiple = True; neither snv or cnv are filled in sample.tsv; neither multiple snv or cnv are filled in conf.ini
    elif multiple and not (snv_file_path or cnv_file_path) and (conf_snv == "" or conf_cnv == ""):
        logger.critical("-m was selected but Muliple section in conf.ini wasn't filled and the file doesn't looks like a multiVCF.")
        raise Exception ("Input error")
    #CASE 3: multiple = False; snv or cnv are filled in sample.tsv; multiple snv or cnv are filled in conf.ini
    elif not multiple and (snv_file_path or cnv_file_path) and not (conf_snv == "" or conf_cnv == ""):
        logger.critical("-m was not selected but both sample.tsv and Muliple section in conf.ini were filled.")
        raise Exception ("Input error")
    #CASE 4: multiple = False; snv or cnv are filled in sample.tsv; neither multiple snv or cnv are filled in conf.ini
    elif not multiple and (snv_file_path or cnv_file_path) and (conf_snv == "" or conf_cnv == ""):
        logger.warning("SNV and/or CNV columns in sample.tsv were not filled.")
        # raise Exception ("Input error")


def check_multiple_folder(input_dir, multiple):
    
    snv_mulitple, snv_single, cnv_mulitple, cnv_single = False, False, False, False

    snv_folder = os.path.join(input_dir, "SNV")
    multiple_vcf_snv = [file for file in os.listdir(snv_folder) if file.endswith(".vcf")][0]
    snv_file = os.path.join(input_dir, "sample_id_snv.txt")
    cmd_snv = "vcf-query -l " + os.path.join(snv_folder, multiple_vcf_snv) + " > " + snv_file
    os.system(cmd_snv)
    with open (snv_file, "r") as file:
        lines = file.readlines()
        if len(lines) >= 2 and not multiple:
            snv_mulitple = True
            logger.error("-m option was not selected but the SNV.vcf file is multiple!")
        elif len(lines) < 2 and multiple:
            snv_single = True
            logger.error("-m option was selected but the SNV.vcf file is not multiple!")        
    os.remove(snv_file)

    cnv_folder = os.path.join(input_dir, "CNV")
    multiple_vcf_cnv = [file for file in os.listdir(cnv_folder) if file.endswith(".vcf")][0]
    cnv_file = os.path.join(input_dir, "sample_id_cnv.txt")
    cmd_cnv = "vcf-query -l " +  os.path.join(cnv_folder, multiple_vcf_cnv) + " > " + cnv_file
    os.system(cmd_cnv)
    with open (cnv_file, "r") as file:
        lines = file.readlines()
        if len(lines) >= 2 and not multiple:
            cnv_mulitple = True
            logger.error("-m option was not selected but the CNV.vcf file is multiple!")
        elif len(lines) < 2 and multiple:
            cnv_single = True
            logger.error("-m option was selected but the CNV.vcf file is not multiple!")   
    os.remove(cnv_file)

    if cnv_mulitple or snv_mulitple:
        raise Exception ("Input error")
    if cnv_single or snv_single:
        raise Exception ("Input error")


def write_clinical_sample(clin_samp_path, output_folder, table_dict):

    logger.info("Writing data_clinical_sample.txt file...")

    conf_header_short = config.get('ClinicalSample', 'HEADER_SAMPLE_SHORT')
    conf_header_long = config.get('ClinicalSample', 'HEADER_SAMPLE_LONG')
    conf_header_type = config.get('ClinicalSample', 'HEADER_SAMPLE_TYPE')

    data_clin_samp = pd.read_csv(clin_samp_path, sep="\t", header=0, dtype=str)

    data_clin_samp.drop(['snv_path', 'cnv_path', 'comb_path'], axis=1, inplace=True)
    
    data_clin_samp.columns = data_clin_samp.columns.str.upper()

    combout_df = pd.DataFrame.from_dict(table_dict).transpose().reset_index()
    combout_df = combout_df.rename(columns={"index": "SAMPLE_ID", 0: "PATIENT_ID"})
    
    final_data_sample = data_clin_samp
    
    if len(combout_df.columns) > 2:
        data_clin_samp.drop(columns=["MSI_THR", "TMB_THR"], inplace=True)
        combout_df = combout_df.rename(columns={1: "MSI", 2: "TMB", 3: "MSI_THR", 4:"TMB_THR"})
        if any(data_clin_samp["MSI"] != combout_df["MSI"]) or any(data_clin_samp["TMB"] != combout_df["TMB"]):
            logger.warning("MSI and/or TMB values are reported in sample.tsv and CombinedOutput but they mismatch! CombinedOutput values were selected by default")
            data_clin_samp.drop(columns=["MSI", "TMB"], inplace=True)
      
        final_data_sample = pd.merge(data_clin_samp, combout_df, on=["PATIENT_ID", "SAMPLE_ID"])

    new_cols = ["SAMPLE_ID", "PATIENT_ID", "MSI", "TMB", "MSI_THR", "TMB_THR", "ONCOTREE_CODE"]+final_data_sample.columns[3:-4].to_list()
    final_data_sample = final_data_sample[new_cols]
    
    dataclin_columns = list(final_data_sample.columns)

    # Add header's fifth row
    default_row = pd.DataFrame([dataclin_columns], columns=dataclin_columns)
    final_data_sample = pd.concat([default_row, final_data_sample], ignore_index=True) 

    # Add header's fourth row (SERIES OF 1s)
    header_numbers = pd.DataFrame([[1] * len(final_data_sample.columns)], columns=dataclin_columns)
    final_data_sample = pd.concat([header_numbers, final_data_sample], ignore_index=True)
    
    # Add header's third row (HEADER_SAMPLE_TYPE)
    if not conf_header_type:
        header_row = ["STRING"] * len(dataclin_columns)
        header_row[dataclin_columns.index('MSI')] = "NUMBER"
        header_row[dataclin_columns.index('TMB')] = "NUMBER"
        sample_header_type = pd.DataFrame([header_row], columns=dataclin_columns)
    else:
        types_list = conf_header_type.split(',')
        types_list = list(map(lambda x: x.strip(), types_list))
        for type in types_list:
            if type.upper() not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{type} is not a valid type. Please check the given input in conf.ini. Valid types: STRING, NUMBER, BOOLEAN")
                raise(NameError("The type is not valid: exiting from walk script!"))
        try:
            types_list = list(map(lambda x: x.upper(), types_list))
            #types_list = types_list + ["NUMBER", "NUMBER", "STRING", "STRING"]
            
            sample_header_type = pd.DataFrame([types_list], columns=dataclin_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(types_list)}) in HEADER_SAMPLE_TYPE is different from the effective number of columns ({len(final_data_sample.columns)}).") 
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_sample = pd.concat([sample_header_type, final_data_sample], ignore_index=True)   

    # Add header's second row (HEADER_SAMPLE_LONG)
    if not conf_header_long:
        sample_header_long = default_row
    else:
        try:
            combined_headers = conf_header_long.split(",")
            sample_header_long = pd.DataFrame([combined_headers], columns=dataclin_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(combined_headers)}) in HEADER_SAMPLE_LONG in conf.ini is different from the effective number of columns ({len(final_data_sample.columns)}).") 
            raise(NameError("Different number of columns: exiting from walk script!"))    
    final_data_sample = pd.concat([sample_header_long, final_data_sample], ignore_index=True)

    # Add header's first row (HEADER_SAMPLE_SHORT)
    if not conf_header_short:
        sample_header_short = default_row
    else:
        try:
            combined_headers = conf_header_short.split(",")
            sample_header_short = pd.DataFrame([combined_headers], columns=dataclin_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(combined_headers)}) in HEADER_SAMPLE_SHORT in conf.ini is different from the effective number of columns ({len(final_data_sample.columns)}).") 
            raise(NameError("Different number of columns: exiting from walk script!")) 
    final_data_sample = pd.concat([sample_header_short, final_data_sample], ignore_index=True)
    
    #final_data_sample.loc[4].replace({'SAMPLEID': 'SAMPLE_ID', 'PATIENTID': 'PATIENT_ID'}, inplace=True)
    #final_data_sample.replace({'SAMPLEID': 'SAMPLE_ID', 'PATIENTID': 'PATIENT_ID'}, inplace=True)
    final_data_sample.loc[0:3, 'SAMPLE_ID'] = final_data_sample.loc[0:3, 'SAMPLE_ID'].apply(lambda x: f'#{x}')

    data_clin_txt = os.path.join(output_folder, 'data_clinical_sample.txt')
    final_data_sample.to_csv(data_clin_txt, sep="\t", index=False, header=False)


def write_default_clinical_patient(output_folder, table_dict):
    logger.info("Writing data_clinical_patient.txt file...")
    data_clin_samp = os.path.join(output_folder, 'data_clinical_patient.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tAge\tGender\n')
    cil_sample.write('#Patient identifier\tAge\tGender\n')
    cil_sample.write('#STRING\tNUMBER\tSTRING\n')
    cil_sample.write('#1\t1\t1\n')
    cil_sample.write('PATIENT_ID\tAGE\tGENDER\n')

    nested_list = list(table_dict.values())
    list_patients = flatten(nested_list)
    list_patients = set(list_patients)
    for v in list_patients:
        cil_sample.write(v + "\tNaN\tNaN\n")
    cil_sample.close()

def add_header_patient_type(sample_pzt, datapat_columns, conf_header_type, final_data_pat):
    if not conf_header_type:
        def_type = ["STRING", "NUMBER", "STRING"]
        header_type = def_type + ["STRING"] * (len(datapat_columns)-3)
        header_type_df = pd.DataFrame([header_type], columns=datapat_columns)
    else:
        types_list = conf_header_type.split(',')
        types_list = list(map(lambda x: x.strip().upper(), types_list))
        for type in types_list:
            if type not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{type} is not a valid type. Please check the given input in conf.ini. " +
                                    "Valid types: STRING, NUMBER, BOOLEAN")
                raise(NameError("The type is not valid: exiting from walk script!"))
        try:
            header_type_df = pd.DataFrame([types_list], columns=datapat_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(types_list)}) in HEADER_PATIENT_TYPE " +
                                f"is different from the effective number of columns ({len(datapat_columns)}) in {sample_pzt}.") 
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_pat = pd.concat([header_type_df, final_data_pat], ignore_index=True)
    return final_data_pat

def add_header_patient_short(sample_pzt, datapat_columns, conf_header_short, default_row, final_data_pat):
    if not conf_header_short:
        pat_header_short = default_row
    else:
        try:
            pat_header_short = pd.DataFrame([conf_header_short.split(', ')], columns=datapat_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(conf_header_short.split(', '))}) in HEADER_PATIENT_SHORT in conf.ini " +
                                f"is different from the effective number of columns ({len(datapat_columns)}) in {sample_pzt}.") 
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_pat = pd.concat([pat_header_short, final_data_pat], ignore_index=True)
    return final_data_pat

def add_header_patient_long(sample_pzt, datapat_columns, conf_header_long, default_row, final_data_pat):
    if not conf_header_long:
        pat_header_long = default_row
    else:
        try:
            pat_header_long = pd.DataFrame([conf_header_long.split(', ')], columns=datapat_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(conf_header_long.split(', '))}) in HEADER_PATIENT_LONG " +
                                f"in conf.ini is different from the effective number of columns ({len(datapat_columns)}) in {sample_pzt}.") 
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_pat = pd.concat([pat_header_long, final_data_pat], ignore_index=True)
    return final_data_pat

def write_clinical_sample_empty(output_folder, table_dict):
    logger.info("Writing data_clinical_sample.txt file...")
    data_clin_samp = os.path.join(output_folder, 'data_clinical_sample.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tSample Identifier\n')
    cil_sample.write('#Patient identifier\tSample Identifier\n')
    cil_sample.write('#STRING\tSTRING\n')
    cil_sample.write('#1\t1\n')
    cil_sample.write('PATIENT_ID\tSAMPLE_ID\n')
    for k, v in table_dict.items():
        cil_sample.write(v[0] + '\t' + k + '\n')
    cil_sample.close()

# def check_cna_vcf(file, inputFolderCNV, multivcf):
#     vcf=pd.read_csv(os.path.join(inputFolderCNV, file), comment="#", sep="\t", \
#                     names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"])
        
#     nsample = os.popen(f'bcftools query -l {os.path.join(inputFolderCNV,file)}').read().split("\n")
#     nsample = [sample for sample in nsample if not sample == ""]
#     if len(nsample)>1 and not multivcf:
#         logger.critical("VCF contains multiple samples")
#         exit(1)

#     if vcf.loc[0]["FORMAT"] == "FC":
#         return True
#     else:
#         return False
    
# def check_snv_vcf(file,inputFolderSNV, multivcf):
#     vcf = pd.read_csv(os.path.join(inputFolderSNV, file), comment="#", sep="\t", \
#                     names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"])
    
#     nsample = os.popen(f'bcftools query -l {os.path.join(inputFolderSNV, file)}').read().split("\n")
#     nsample = [sample for sample in nsample if not sample == ""]
#     if len(nsample)>1 and not multivcf:
#         logger.critical("VCF contains multiple samples")
#         exit(1)
    
#     if vcf.loc[0]["FORMAT"].startswith("GT"):
#         return True
#     else:
#         return False
   

def extract_multiple_cnv(multiple_vcf, input_dir):
    single_sample_vcf_dir = os.path.join(input_dir, "single_sample_vcf")  
    if not os.path.exists(os.path.join(input_dir, "single_sample_vcf")):
        os.mkdir(os.path.join(input_dir, "single_sample_vcf"))
    cmd = "vcf-query -l " + multiple_vcf + " > " + os.path.join(input_dir, "sample_id.txt")
    os.system(cmd)

    with open (os.path.join(input_dir, "sample_id.txt"), "r") as file:
        lines = file.readlines()
        for sample in lines:
            sample = sample.strip()
            single_sample_vcf = os.path.join(single_sample_vcf_dir, f"{sample}.vcf")
            cmd = f"vcftools --vcf {multiple_vcf} --indv {sample} --recode --recode-INFO-all --stdout > {single_sample_vcf}"
            os.system(cmd)
    
def extract_multiple_snv(multiple_vcf,input_dir):    
    if not os.path.exists(os.path.join(input_dir, "single_sample_vcf")):
        os.mkdir(os.path.join(input_dir, "single_sample_vcf"))

    cmd = "vcf-query -l " + multiple_vcf + " > " + os.path.join(input_dir, "sample_id.txt")
    os.system(cmd)

    with open (os.path.join(input_dir, "sample_id.txt"), "r") as file:
        lines = file.readlines()
        for sample in lines:
            print(sample.strip())
            cmd="vcf-subset --exclude-ref -c " + sample.strip() + " " + multiple_vcf + " > " + os.path.join(os.path.join(input_dir, "single_sample_vcf"), sample.strip() + ".vcf")
            os.system(cmd)

def check_field_tsv(row, name):
    try:
        field=str(row[name])
    except KeyError as e: 
        logger.critical(f"KeyError: {e} not found! Check if column name is correctly spelled or if there are tabs/spaces before or after the coloumn key: \n{row.index}. \nThis error may also occur if the table columns have not been separated by tabs!")
        sys.exit()
        #raise(KeyError("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))
    return field

def get_combinedVariantOutput_from_folder(inputFolder, file, isinputfile):
    combined_dict = dict()

    for _,row in file.iterrows():
        
        #patientID = check_field_tsv(row, "PatientID")
        # combined_file = patientID+COMBOUT #da verificare
        # combined_path = os.path.join(inputFolder,"CombinedOutput",combined_file)
        sampleID = check_field_tsv(row, "SAMPLE_ID")
        patientID = check_field_tsv(row, "PATIENT_ID")
        if isinputfile:
            combined_path = check_field_tsv(row, "comb_path")
        else: # TODO rivedere se mantenere così o meno
            combined_path = os.path.join(inputFolder, "CombinedOutput", patientID + "_CombinedVariantOutput.tsv")
        
        if os.path.exists(combined_path):
            pass 
        else:
            logger.warning(f"{combined_path} not exists")
        combined_dict[sampleID] = combined_path
    return combined_dict

def check_input_file(output_folder, file, copy_to):
    if os.path.exists(file):
        os.system("cp " + file + " " + os.path.join(output_folder, "temp", copy_to))
    elif not file:
        logger.warning(f"No final_path set in conf.ini!")
    else:
        logger.warning(f"{file} not found")
        
def check_folders(output_folder, snv_path, cnv_path, combout):
    check_input_file(output_folder, snv_path, "SNV")
    check_input_file(output_folder, cnv_path, "CNV")
    check_input_file(output_folder, combout, "CombinedOutput")
        
def transform_input(tsv, clin_pzt, fusion_tsv, output_folder, multiple):
    os.makedirs(os.path.join(output_folder, "temp"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "SNV"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "CNV"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "CombinedOutput"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "FUSIONS"), exist_ok=True)
   
    os.system("cp " + tsv + " " + os.path.join(output_folder, "temp", "sample.tsv"))
    if clin_pzt!="":
        os.system("cp " + clin_pzt + " " + os.path.join(output_folder, "temp", "patient.tsv"))
    if fusion_tsv!="":
        os.system("cp " + fusion_tsv + " " + os.path.join(output_folder, "temp", "FUSIONS", "fusion.tsv"))

    if multiple:
        snv_path = config.get('Multiple', 'SNV')
        cnv_path = config.get('Multiple', 'CNV')
        combout = config.get('Multiple', 'COMBOUT')
        
        check_folders(output_folder, snv_path, cnv_path, combout)
    
    else:   
        tsv_file = pd.read_csv(tsv, sep="\t", dtype="string", keep_default_na=False)
            
        for _,row in tsv_file.iterrows():

            # Decommentare se si fornsicono path come colonne dell'input
            snv_path = row["snv_path"]
            cnv_path = row["cnv_path"]
            combout = row["comb_path"]

            # Decommentare se deve costruire i path in base alle posizioni delle cartelle
            # #res_folder = "/data/data_storage/novaseq_results"
            # res_folder = INPUT_PATH
            # snv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+SNV)     #"_MergedSmallVariants.genome.vcf")
            # #snv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+"_MergedSmallVariants.genome.vcf")
            # cnv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+CNV) # "_CopyNumberVariants.vcf")
            # combout=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["PatientID"]+COMBOUT)

            check_folders(output_folder, snv_path, cnv_path, combout)
                
    return os.path.join(output_folder, "temp")


def fill_fusion_from_temp(input, fusion_table_file, clin_file, fusion_files):

    nfusion=len(fusion_files)
    logger.info(f"Found {nfusion} Fusion file(s) ")

    fusion_input = os.path.join(input, "FUSIONS", fusion_files[0])
    with open(fusion_input) as template:
        header = template.readline()
            
    with open(fusion_table_file, "w") as fusion_table:

        fusion_table.write(header)
                
        for fusion_file in fusion_files:
            ff = pd.read_csv(fusion_input, sep="\t")
                    
            if not set(["Sample_Id","SV_Status","Site1_Hugo_Symbol","Site2_Hugo_Symbol"]).issubset(ff.columns):
                logger.warning(f"{fusion_file} does not contain required columns")
                continue
                
            if ff.shape[0]==0:
                logger.info(f"No Fusions found in {fusion_file}")
                continue
            else:
                logger.info(f"Fusions found in {fusion_file}")          

            for fus in ff.itertuples(index=False):
                if fus.Sample_Id in clin_file["SAMPLE_ID"].values:
                    if int(fus.Normal_Paired_End_Read_Count) >= 15:
                        fusion_table.write('\t'.join(map(str, fus)) + '\n')

            #TODO controllare possibili campi (FUSION, INVERSION, DELETION)

def annotate_fusion(cancer, fusion_table_file, data_sv, input_file):
    
    if "ONCOTREE_CODE" in input_file.columns:
        input_file["SAMPLE_ID"] = input_file["SAMPLE_ID"] + ".bam"
        fusion_table_df = data_sv.merge(input_file, how="inner", left_on="Sample_Id", right_on="SAMPLE_ID")
        fusion_table_df.to_csv(fusion_table_file, sep="\t", index=False)
        fusion_table_file_out = fusion_table_file.replace(".txt", "ann.txt")
        os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                        -o {fusion_table_file_out} -b {config.get('OncoKB', 'ONCOKB')}")   

    else:   
        fusion_table_file_out = fusion_table_file.replace(".txt", "ann.txt")       
        os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                            -o {fusion_table_file_out} -t {cancer.upper()}  -b {config.get('OncoKB', 'ONCOKB')}")
                    
    return fusion_table_file_out

def fill_fusion_from_combined(fusion_table_file, combined_dict, THR_FUS):
    logger.info(f"Creating data_sv.txt file...")
    
    with open(fusion_table_file, "w") as fusion_table:

        header = 'Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\tNormal_Paired_End_Read_Count\tEvent_Info\tRNA_Support\n'
        fusion_table.write(header)

        for k, v in combined_dict.items():
            fusions=[]
            #logger.info(f"Reading Fusion info in CombinedOutput file {v}...")
            try:
                fusions = tsv.get_fusions(v)
            except Exception as e:
                logger.error(f"Something went wrong while reading Fusion section of file {v}")
            if len(fusions) == 0:
                #logger.info(f"No Fusions found in {v}")
                continue
            # else:
            #     logger.info(f"Fusions found in {v}")

            for fus in fusions:
                if len(fusions) > 0:
                    Site1_Hugo_Symbol = fus['Site1_Hugo_Symbol']
                    Site2_Hugo_Symbol = fus['Site2_Hugo_Symbol']
                    if Site2_Hugo_Symbol == 'CASC1':
                        Site2_Hugo_Symbol = 'DNAI7'
                    Site1_Chromosome = fus['Site1_Chromosome']
                    Site2_Chromosome = fus['Site2_Chromosome']
                    Site1_Position = fus['Site1_Position']
                    Site2_Position = fus['Site2_Position']

                    if eval("int(fus['Normal_Paired_End_Read_Count'])" + THR_FUS):
                        fusion_table.write(k+'\tSOMATIC\tFUSION\t'+\
            str(Site1_Hugo_Symbol)+'\t'+str(Site2_Hugo_Symbol)+'\t'+fus['Normal_Paired_End_Read_Count']+\
            '\t'+fus['Event_Info']+' Fusion\t'+'Yes\n')

def fill_from_file(table_dict_patient, fileinputclinical, MSI_THR, TMB):

    #logger.info(f"Reading Tumor clinical parameters info in sample.tsv...")
    for k, m, t in zip(fileinputclinical["SAMPLE_ID"], fileinputclinical["MSI"], fileinputclinical["TMB"]):
        table_dict_patient[k].append(m)
        table_dict_patient[k].append(t)
        
        if np.isnan(float(m)):
            table_dict_patient[k].append("NA")
        else:
            if eval("float(m)" + MSI_THR):
                table_dict_patient[k].append("Stable")   
            else:
                table_dict_patient[k].append('Unstable')
        
        if not np.isnan(float(t)):
            for _k, _v in TMB.items():
                if eval("float(t)" + _v):
                    table_dict_patient[k].append(_k)
                    break
        else:
            table_dict_patient[k].append("NA")
    return table_dict_patient

def fill_from_combined(combined_dict, table_dict_patient, MSI_SITES_THR, MSI_THR, TMB):
    for k, v in combined_dict.items():
        #logger.info(f"Reading Tumor clinical parameters info in CombinedOutput file {v}...")
        try:
            tmv_msi = tsv.get_msi_tmb(v)
        except Exception as e:
            logger.error(f"Something went wrong!")
        #logger.info(f"Tumor clinical parameters Values found: {tmv_msi}")

        if not tmv_msi['MSI'][0][1]=="NA" and eval("float(tmv_msi['MSI'][0][1])" + MSI_SITES_THR):    
            table_dict_patient[k].append(tmv_msi['MSI'][1][1])   
        else:
            table_dict_patient[k].append('NA')
        table_dict_patient[k].append(tmv_msi['TMB_Total'])
        if not tmv_msi['MSI'][0][1]=="NA":
            if not tmv_msi['MSI'][1][1] =="NA":
                if eval("float(tmv_msi['MSI'][1][1])" + MSI_THR):
                    table_dict_patient[k].append("Stable")   
                else:
                    table_dict_patient[k].append('Unstable')

                found = False
            else:
                table_dict_patient[k].append('NA')
        else:
            table_dict_patient[k].append('NA')
        found = False

        if not tmv_msi["TMB_Total"]=="NA":
            for _k, _v in TMB.items():
                if eval(tmv_msi["TMB_Total"] + _v):
                    table_dict_patient[k].append(_k)
                    found = True
                    break
        else:
            table_dict_patient[k].append("NA")
        # if found==False:
        #     table_dict_patient[k].append("NA")
        #     table_dict_patient[k].append(list(TMB.keys())[-1])
    return table_dict_patient

def input_extraction(input):
    sample_tsv=input[0]
    patient_tsv, fusion_tsv = "",""
    
    if len(input)>1:
        patient_tsv=input[1].strip()
        if len(input)>2:
            fusion_tsv=input[2].strip()

    return sample_tsv, patient_tsv, fusion_tsv

def validate_input(oncokb, vcf_type, filters, cancer):
    
    #check that oncokb key is filled in conf.ini when oncokb annotation is selected
    if oncokb:
        assert config.get('OncoKB', 'ONCOKB')!="", \
               f"oncokb option was set but ONCOKB field in conf.ini is empty!"
    
    #check that vep info in conf.ini are set when snv analysisi is request
    if vcf_type==None or "snv" in vcf_type:
        assert VEP_PATH !="" and VEP_DATA !="", \
               f"vep_path and/or vep_data field in conf.ini is empty!"
        assert REF_FASTA !="", \
               f"re_fasta field in conf.ini is empty!"
    
    #verify that oncokb filter function only when oncokb annotation is set 
    if "o" in filters and oncokb==False:
        logger.warning("OncoKB filter was selected in filters options but -k option was not set. This filtering will be ignored.")
    
    #check if cancer id is compatible with cbioportal
    cancer_cbio=pd.read_csv("cancer_list.txt", sep="\t")
    cancer_cbio=cancer_cbio["TYPE_OF_CANCER_ID"].values.tolist()
    
    if cancer not in cancer_cbio:
        logger.critical(f"cancer_id '{cancer}' is not recognize by cbioportal. Check in cancer_list.txt to find the correct cancer id")
        sys.exit()

def write_filters_in_report(output_folder):
    now = datetime.now()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")
    report_file_path = os.path.join(output_folder, "report.txt")  
    sections_to_include = {
            "Filters": ["BENIGN", "CLIN_SIG", "CONSEQUENCES", "ONCOKB_FILTER", 
                        "t_VAF_min", "t_VAF_min_novel", "t_VAF_max", 
                        "AF", "POLYPHEN", "IMPACT", "SIFT"],
            "Cna": ["HEADER_CNV", "PLOIDY"],
            "TMB": ["THRESHOLD"],
            "MSI": ["THRESHOLD_SITES", "THRESHOLD"],
            "FUSION": ["THRESHOLD"]
        }

    conf_content = []

    for section, keys in sections_to_include.items():
        if section in config:
            conf_content.append(f"[{section}]")  
            for key in keys:
                if key in config[section]:  
                    value = config[section][key]
                    conf_content.append(f"{key.upper()} = {value}")  

    conf_content = "\n".join(conf_content) + "\n\n"

    with open(report_file_path, "r") as file:
        val_report = file.readlines()

    with open(report_file_path, "w") as file:
        file.write(f"Varan run - {date}\n\nThe following configuration and filters have been used:\n")
        file.write(conf_content)
        file.write("This is the report from cBioPortal Validator:\n")
        file.writelines(val_report)



def walk_folder(input, multiple, output_folder, oncokb, cancer, overwrite_output=False, resume=False, vcf_type=None, filters=""):
    
    logger.info("Starting walk_folder script:")
    logger.info(f"walk_folder args [input:{input}, output_folder:{output_folder}, overwrite:{overwrite_output}, resume:{resume}, vcf_type:{vcf_type}, filters:{filters}, multiple:{multiple}]")
 
    config.read("conf.ini")
    
    input, sample_pzt, fusion_tsv = input_extraction(input)
    assert os.path.exists(input), \
            f"No valid file/folder {input} found. Check your input path"

    validate_input(oncokb, vcf_type, filters, cancer)
    
    
    
    ###############################
    ###      OUTPUT FOLDER      ###
    ###############################
    
    if not resume or not os.path.exists(os.path.join(output_folder, "temp")):
        output_folder = create_folder(output_folder, overwrite_output, resume)
    
    global isinputfile 
    if os.path.isdir(input):
        isinputfile = False
        input_folder = input
        check_multiple_folder(input_folder, multiple)
   
    elif os.path.isfile(input):
        isinputfile = True
        check_multiple_file(input, multiple)
        input_folder = transform_input(input, sample_pzt, fusion_tsv, output_folder, multiple)
    #input=os.path.join(output_folder, "temp")

    else:
        logger.critical(f"The input {input} isn't a file nor a folder")
        raise(FileNotFoundError("Exiting from walk script!"))

    inputFolderSNV = os.path.abspath(os.path.join(input_folder, "SNV"))
    inputFolderCNV = os.path.abspath(os.path.join(input_folder, "CNV"))
    inputFolderCombOut = os.path.abspath(os.path.join(input_folder, "CombinedOutput"))
    inputFolderFusion = os.path.abspath(os.path.join(input_folder, "FUSIONS"))

    
    assert (len(os.listdir(inputFolderCNV))>0 or len(os.listdir(inputFolderCNV))>0 or len(os.listdir(inputFolderCombOut))>0), \
        "No valid input file was found for neither SNV, CNV or CombinedOutput! Check your input file/folder and input options."
       
    if os.path.exists(inputFolderCNV) and not vcf_type in ["snv", "fus", "tab"]:
        if multiple:
            multivcf = [i for i in os.listdir(inputFolderCNV) if i.endswith('.vcf')][0]
            extract_multiple_cnv(os.path.join(inputFolderCNV, multivcf), inputFolderCNV)
            inputFolderCNV= os.path.join(inputFolderCNV, "single_sample_vcf")
        logger.info("Check CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(inputFolderCNV)
        logger.info("Everything ok!")

    if os.path.exists(inputFolderSNV) and not vcf_type in ["cnv", "fus", "tab"]:
        os.makedirs(tmp, exist_ok=True)
        if multiple:
            multivcf = [i for i in os.listdir(inputFolderSNV) if i.endswith('.vcf')][0]
            extract_multiple_snv(os.path.join(inputFolderSNV, multivcf), inputFolderSNV)
            inputFolderSNV = os.path.join(inputFolderSNV, "single_sample_vcf")
        logger.info("Check SNV files...")
    case_folder_arr = get_snv_from_folder(inputFolderSNV)
    logger.info("Everything ok!")

    ###############################
    ###       SNV AND CNV       ###
    ###############################
    
    logger.info("Managing CNV files...")
    if os.path.exists(inputFolderCNV) and not vcf_type in ["snv", "fus", "tab"]:
        sID_path_cnv = cnv_type_from_folder(input_folder, case_folder_arr_cnv, output_folder, oncokb, cancer, multiple)
    
    logger.info("Managing SNV files...")
    if os.path.exists(inputFolderSNV) and not vcf_type in ["cnv", "fus", "tab"]:
        sID_path_snv = snv_type_from_folder(inputFolderSNV, case_folder_arr)
        
        logger.info("Check maf folder...")
        maf_path = os.path.join(output_folder, "maf")
        if os.path.isdir(maf_path) and len([i for i in os.listdir(maf_path) if i.endswith('.maf')])>0:
            logger.info("a non empty maf folder already exists!")
        
        if not resume:
            logger.info("Starting vcf2maf conversion...")
            if "d" in filters:
                logger.info("filtering out vcfs with dots in ALT column")
                sID_path_snv = vcf_filtering(sID_path_snv, output_folder, output_filtered)
            temporary = create_random_name_folder()
            for k, v in sID_path_snv.items():
                cl = vcf2maf_constructor(k, v, temporary, output_folder)
                run_vcf2maf(cl)
    
    logger.info("Clearing scratch folder...")
    clear_scratch()

   
    ###############################
    ###       GET FUSION        ###
    ###############################

    clin_sample_path = os.path.join(input_folder, "sample.tsv")
    fusion_table_file = os.path.join(output_folder, 'data_sv.txt')
    
    try:
        clin_file = pd.read_csv(clin_sample_path, sep="\t", dtype=str)
    except Exception as e:
        logger.critical(f"Something went wrong while reading {clin_sample_path}!")
        raise(Exception("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))

    if os.path.exists(os.path.join(input_folder, "CombinedOutput")) and len(os.listdir(os.path.join(input_folder, "CombinedOutput")))>0 and not type in ["cnv","snv","tab"]:
        THR_FUS = config.get('FUSION', 'THRESHOLD')
        combined_dict = get_combinedVariantOutput_from_folder(input_folder, clin_file, isinputfile)       
        fill_fusion_from_combined(fusion_table_file, combined_dict, THR_FUS)
    
    elif os.path.exists(os.path.abspath(os.path.join(input_folder,"FUSIONS"))) and not type in ["cnv","snv","tab"]:    
        fusion_files=[file for file in os.listdir(os.path.join(input_folder,"FUSIONS")) if "tsv" in file]
        fill_fusion_from_temp(input_folder, fusion_table_file, clin_file, fusion_files)  

    if oncokb and os.path.exists(fusion_table_file):
        data_sv = pd.read_csv(fusion_table_file, sep="\t")
        input_file = pd.read_csv(clin_sample_path, sep="\t")
        fusion_table_file_out = annotate_fusion(cancer, fusion_table_file, data_sv, input_file)
        
        if "o" in filters:
            fus_file = pd.read_csv(fusion_table_file_out, sep="\t")
            fus_file = filter_OncoKB(fus_file)
            fus_file.to_csv(fusion_table_file_out, index=False, sep="\t")
        os.system(f"mv {fusion_table_file_out}  {fusion_table_file}") 
            
      
    ##############################
    ##       MAKES TABLE       ###
    ##############################    

    table_dict_patient = get_table_from_folder(clin_sample_path)
    
    logger.info("Writing clinical files...")
    
    if sample_pzt.strip()!="":
        logger.info("Writing data_clinical_patient.txt file...")
        
        input_file_path = os.path.join(input_folder, "patient.tsv")#os.path.basename(sample_pzt))
        data_clin_pat = pd.read_csv(input_file_path, sep="\t", header=0, dtype=str)
        
        data_clin_pat.columns = data_clin_pat.columns.str.upper()
        datapat_columns = list(data_clin_pat.columns)
        
        # Get headers from conf.ini if they're present
        conf_header_short = config.get('ClinicalPatient', 'HEADER_PATIENT_SHORT')
        conf_header_long = config.get('ClinicalPatient', 'HEADER_PATIENT_LONG')
        conf_header_type = config.get('ClinicalPatient', 'HEADER_PATIENT_TYPE')

        # Add header's fifth row
        default_row = pd.DataFrame([datapat_columns], columns=datapat_columns)
        final_data_pat = pd.concat([default_row, data_clin_pat], ignore_index=True) 

        # Add header's fourth row (1s)
        header_numbers = pd.DataFrame([[1] * len(datapat_columns)], columns=datapat_columns)
        final_data_pat = pd.concat([header_numbers, final_data_pat], ignore_index=True)

        # Add header's third row (HEADER_PATIENT_TYPE)
        final_data_pat = add_header_patient_type(sample_pzt, datapat_columns, conf_header_type, final_data_pat)

        # Add header's second row (HEADER_PATIENT_LONG)
        final_data_pat = add_header_patient_long(sample_pzt, datapat_columns, conf_header_long, default_row, final_data_pat)

        # Add header's first row (HEADER_PATIENT_SHORT)
        final_data_pat = add_header_patient_short(sample_pzt, datapat_columns, conf_header_short, default_row, final_data_pat)
        
        #final_data_pat.replace({'PATIENT ID': 'PATIENT_ID'}, inplace=True)
        final_data_pat.loc[0:3, 'PATIENT_ID'] = final_data_pat.loc[0:3, 'PATIENT_ID'].apply(lambda x: f'#{x}')

        data_clin_txt = os.path.join(output_folder, 'data_clinical_patient.txt')
        final_data_pat.to_csv(data_clin_txt, sep="\t", index=False, header=False)

    else:
        write_default_clinical_patient(output_folder, table_dict_patient)

    fileinputclinical = pd.read_csv(os.path.join(input_folder, "sample.tsv"), sep="\t", index_col=False, dtype=str)
    
    MSI_THR = config.get('MSI', 'THRESHOLD')
    TMB_THR = ast.literal_eval(config.get('TMB', 'THRESHOLD'))
    
    if os.path.exists(os.path.join(input_folder, "CombinedOutput")) and \
        len(os.listdir(os.path.join(input_folder, "CombinedOutput"))) > 0:
        MSI_SITES_THR = config.get('MSI', 'THRESHOLD_SITES')

        combined_dict = get_combinedVariantOutput_from_folder(input_folder, clin_file, isinputfile)        
        table_dict_patient = fill_from_combined(combined_dict, table_dict_patient, MSI_SITES_THR, MSI_THR, TMB_THR)
    else:
        table_dict_patient = fill_from_file(table_dict_patient, fileinputclinical, MSI_THR, TMB_THR)

    write_clinical_sample(clin_sample_path, output_folder, table_dict_patient)

    #DECOMMENTARE UNA VOLTA FINITI TEST
    # if isinputfile:
    #     try:
    #         logger.info("Deleting temp folder")
    #         shutil.rmtree(os.path.join(output_folder,"temp"))
    #     except Exception as e:
    #         print("No combined output found")
    #         write_clinical_sample_empty(output_folder, table_dict_patient)
 
    # logger.info("Deleting temp folder")
    # clear_temp(output_folder)
    
    logger.success("Walk script completed!\n")

    return output_folder, input, fusion_tsv
        
    
        
    