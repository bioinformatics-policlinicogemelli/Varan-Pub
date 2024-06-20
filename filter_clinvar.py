import os
import ast
import argparse
import concatenate
import pandas as pd
from configparser import ConfigParser
from loguru import logger
import sys
import shutil

config = ConfigParser()
configFile = config.read("conf.ini")


def print_unique_clin_sig(df):
    unique_clin_sig = df['CLIN_SIG'].unique()
    print(unique_clin_sig)


def filter_OncoKB(df):
    oncokb_filter=ast.literal_eval(config.get('Filters', 'ONCOKB_FILTER'))
    
    df_filtered=df[df["ONCOGENIC"].isin(oncokb_filter)]
    return df_filtered



def filter_benign(df):
    benign_filter = ~df['CLIN_SIG'].str.contains(config.get('Filters', 'BENIGN')
            , case=False
            , na=False
            , regex=True)
    return df[benign_filter]


def keep_risk_factors(df):

    benign_filter = ~df['CLIN_SIG'].str.contains(config.get('Filters', 'BENIGN')
        , case=False
        , na=False
        , regex=True)
    df=df[benign_filter]
    df = df[df.apply(check_CLIN_SIG,axis=1)|df.apply(check_consequences,axis=1)]
    return df


def filter_PASS(df):
    pass_filter=ast.literal_eval(config.get('Filters', 'FILTER'))
    
    df_filtered=df[df["FILTER"].isin(pass_filter)]
    return df_filtered


  
def check_CLIN_SIG(row):
    clin_sig=ast.literal_eval(config.get('Filters', 'CLIN_SIG'))
    output=[]
    for _e in str(row["CLIN_SIG"]).split(","):
        if _e in clin_sig:
            output.append(True)
        else:
            output.append(False)
    return any(output)
        
def check_consequences(row):
    consequences=ast.literal_eval(config.get('Filters', 'CONSEQUENCES'))
    output=[]
    for _e in str(row['Consequence']).split(","):
        if _e in consequences:
            output.append(True)
        else:
            output.append(False)
    return any(output)

def check_polyphen(row):
    consequences=ast.literal_eval(config.get('Filters', 'POLYPHEN'))
    output=[]
    for _e in str(row['PolyPhen']).split(","):
        if any(_e in s for s in consequences):
            output.append(True)
        else:
            output.append(False)
    return any(output)


def check_sift(row):
    consequences=ast.literal_eval(config.get('Filters', 'SIFT'))
    output=[]
    for _e in str(row['SIFT']).split(","):
        if any(_e in s for s in consequences):
            output.append(True)
        else:
            output.append(False)
    return any(output)

def filter_vf(df):
    t_vaf_min=float(config.get('Filters', 't_VAF_min'))
    t_vaf_max=float(config.get('Filters', 't_VAF_max'))
    
    df = df[(df['t_VF'] > t_vaf_min) | (df['t_VF'].isnull()) | (df['t_VF'] <= t_vaf_max)]
    return df 


def filter_main(input,folder, output_folder,filters,oncokb,cancer,overwrite=False):
    
    logger.info("Starting filter_main script:")
    logger.info(f"filter_main args [maf_folder:{folder}, output_folder:{output_folder},  overwrite:{overwrite}]")

    if os.path.exists(os.path.join(output_folder,'MAF_OncoKB')) and len(os.listdir(os.path.join(output_folder,'MAF_OncoKB')))>0:
        if overwrite:
            logger.warning(f"It seems that the folder 'MAF_OncoKB' already exists. Start removing process...")        
            shutil.rmtree(os.path.join(output_folder,'MAF_OncoKB'))
        else:
            logger.critical(f"The folder 'MAF_OncoKB' already exists. To overwrite an existing folder add the -w option!")
            logger.critical(f"Exit without completing the task!")
            exit()

    file_list = concatenate.get_files_by_ext(folder, 'maf')
    if len(file_list)==0:
        logger.warning(f"The maf folder {os.path.join(folder, 'maf')} seems to be empty! Filtering cannot be done.")
        logger.critical("Empty maf folder: Filter script exited before completing!")
        
    else:
   
        # ANNOTATION    
        if oncokb:
            output_onco=os.path.join(output_folder, 'MAF_OncoKB')
            os.mkdir(output_onco)
            extension="_OncoAnnotated.maf"
        
            tsv_file=[file for file in os.listdir(input) if "tsv" in file][0]
            input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
            for f in file_list:
                root, file = os.path.split(f)
                file_No = file.replace('.maf','') + extension
                file_path = os.path.join(output_onco, file_No)
                if "ONCOTREE_CODE" in input_file.columns:
                    for _ ,row in input_file.iterrows():
                        if row["SampleID"] in file_No:
                            cancer_onco=row["ONCOTREE_CODE"]
                            os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f}\
                                    -o {file_path} -t {cancer_onco.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
                else:               
                    os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f}\
                                -o {file_path} -t {cancer.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
            
        
        # FILTER 

        if oncokb:   
            file_list = concatenate.get_files_by_ext(output_onco, 'maf')
        else:
            file_list=concatenate.get_files_by_ext(folder, 'maf')
            
        
        if not filters==None: 
            os.mkdir(os.path.join(output_folder, 'MAF_Filtered'))
            
            for file in file_list:
                logger.info(f"Filtering file {file}")
                file_to_filter=pd.read_csv(file,sep="\t")
                
                if "f" in filters:
                    file_to_filter[file_to_filter["FILTER"].isin(ast.literal_eval(config.get('Filters', 'FILTER')))]
                
                if "b" in filters:
                    benign_filter = ~file_to_filter['CLIN_SIG'].str.contains(config.get('Filters', 'BENIGN')
                    , case=False
                    , na=False
                    , regex=True)   
                    file_to_filter=file_to_filter[benign_filter]

                if "k" in filters:
                    file_to_filter=filter_OncoKB(file_to_filter)
    
                if "v" in filters:
                    t_vaf_min=float(config.get('Filters', 't_VAF_min'))
                    t_vaf_max=float(config.get('Filters', 't_VAF_max'))
            
                    file_to_filter = file_to_filter[(file_to_filter['t_VF'] > t_vaf_min) | (file_to_filter['t_VF'].isnull()) | (file_to_filter['t_VF'] <= t_vaf_max)]
                
                if "g" in filters:     
                    gnomAD=float(config.get('Filters', 'gnomAD'))
                    file_to_filter =file_to_filter[(file_to_filter['AF'] <gnomAD) | (file_to_filter['AF'].isnull())]
                    
                if "c" in filters:
                    file_to_filter = file_to_filter[file_to_filter.apply(check_CLIN_SIG,axis=1)]
                    
                if "i" in filters:
                    file_to_filter= file_to_filter[file_to_filter["IMPACT"].isin(ast.literal_eval(config.get('Filters',"IMPACT")))]    
                    
                if "q" in filters:
                    file_to_filter = file_to_filter[file_to_filter.apply(check_consequences,axis=1)]
                    
                if "y" in filters:
                    file_to_filter = file_to_filter[file_to_filter.apply(check_polyphen,axis=1)]
                
                if "s" in filters:
                    file_to_filter=file_to_filter[file_to_filter.apply(check_sift,axis=1)]
                    
                logger.info(f"Filtered file: {file}")
                file_to_filter.to_csv(os.path.join(output_folder,"MAF_Filtered",file.split("/")[-1]),sep="\t",index=False)

                
        logger.success("Filter script completed!\n")
