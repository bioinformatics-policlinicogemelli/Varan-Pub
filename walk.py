#####################################
# NAME: walk.py
# Date: 10/01/2023
version = "1.0"
# ===================================

import os
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
import traceback
from versioning import get_newest_version
from filter_clinvar import filter_OncoKB

config = ConfigParser()

configFile = config.read("conf.ini")
OUTPUT_FILTERED = config.get('Paths', 'OUTPUT_FILTERED')
OUTPUT_MAF = config.get('Paths', 'OUTPUT_MAF')
VCF2MAF = config.get('Paths', 'VCF2MAF')
REF_FASTA = config.get('Paths', 'REF_FASTA')
TMP = config.get('Paths', 'TMP')
VEP_PATH = config.get('Paths', 'VEP_PATH')
VEP_DATA = config.get('Paths', 'VEP_DATA')
CLINV = config.get('Paths', 'CLINV')


def create_folder(output_folder,type):
    version="_v1"
    output_folder_version=output_folder+version
    
    if os.path.exists(output_folder_version):
        logger.warning(f"It seems that the folder '{output_folder_version}' already exists.")
        output_folder_version=get_newest_version(output_folder)
    logger.info(f"Creating the output folder '{output_folder_version}' in {os.getcwd()}...")
    os.mkdir(output_folder_version)
    if not type in ["cnv","fus","tab"]:
        maf_path = os.path.join(output_folder_version, 'maf')
        os.mkdir(maf_path)
        filtered_path = os.path.join(output_folder_version, 'snv_filtered')
        os.mkdir(filtered_path)
        logger.info(f"The folder '{output_folder_version}' was correctly created!")
    return output_folder_version

## GET CNV

def check_cna_vcf(file,inputFolderCNV,multivcf):
    
    vcf=pd.read_csv(os.path.join(inputFolderCNV,file),comment="#",sep="\t",header=None)
    
    nsample=os.popen(f'bcftools query -l {os.path.join(inputFolderCNV,file)}').read().split("\n")
   
    nsample=[sample for sample in nsample if not sample ==""]
    if len(nsample)>1 and not multivcf:
        logger.critical("VCF contains multiple samples")
        exit(1)
      
    ### TO CHECK WITH NEW VCF
    if all(vcf.apply(lambda row: row.astype(str).str.contains('FC').any(), axis=1)) or all(vcf.apply(lambda row: row.astype(str).str.contains('SM').any(), axis=1)):
        return True
    else:
        return False


def get_cnv_from_folder(inputFolderCNV,multivcf):
    files= os.listdir(inputFolderCNV)
    cnv_vcf_files=[file for file in files if file.endswith("vcf")]
    check=list(map(lambda x: check_cna_vcf(x,inputFolderCNV,multivcf),cnv_vcf_files))
    incorrect_files=[]
    for i, check_res in enumerate(check):
        if not check_res:
            incorrect_files.append(cnv_vcf_files[i])
    if len(incorrect_files)!=0:
        logger.critical(f"It seems that the files \n{incorrect_files} \nare not CNV! Please check your CNV input data and try again.")
        exit()
    logger.info(f"#{len(cnv_vcf_files)} vcf files found in CNV folder")
    return cnv_vcf_files

def cnv_type_from_folder(input,tsv_input, cnv_vcf_files,oncokb,cancer,output_folder,Filter):
    c = 0
    sID_path = dict()
    for case_folder in cnv_vcf_files:
        if os.path.exists('data_cna_hg19.seg'):
            MODE = 'a'
        else:
            MODE = 'w'
        
        cnv_vcf=case_folder 
        
        tsv_input_file=pd.read_csv(tsv_input,sep="\t")
        sampleID = tsv_input_file[tsv_input_file['SampleID'].apply(lambda x: x in cnv_vcf)]["SampleID"].values.tolist()[0]
        
        if sampleID in sID_path:
            dup = open('sampleID_dup'.log, 'w')
            dup.write(sampleID+'\t'+'cnv_vcf')
            dup.close()
        else:
            sID_path[sampleID] = os.path.join(input,cnv_vcf)
            vcf2tab_cnv.vcf_to_table(sID_path[sampleID], os.path.join(output_folder,'data_cna_hg19.seg'), sampleID, MODE)
            vcf2tab_cnv.vcf_to_table_fc(sID_path[sampleID], os.path.join(output_folder,'data_cna_hg19.seg.fc.txt'), sampleID, MODE)

        c = c +1
    logger.info("Writing data_cna_hg19.seg succefully completed!")
    logger.info("Writing data_cna_hg19.seg.fc.txt succefully completed!")
    
    
    ############################
    ### MANAGE DISCRETE TABLE ##
    ############################

    df_table = pd.read_csv(os.path.join(output_folder,'data_cna_hg19.seg.fc.txt'),sep="\t",header=0)
    result = tabella_to_dict(df_table)

    df = pd.DataFrame()
    for key in result.keys():
        if df.empty:
            for elem in result[key]: 
                genedata = {"Hugo_Symbol":elem[-2],key:elem[-1]}
                temp = pd.DataFrame(genedata,index=[0])
                df = pd.concat([df, temp])
        else:
            valuedf = pd.DataFrame()
            for elem in result[key]: 
                genedata = {key:elem[-1]}
                temp = pd.DataFrame(genedata,index=[0])
                valuedf = pd.concat([valuedf, temp])
            df = pd.concat([df, valuedf], axis=1)

    df.to_csv(os.path.join(output_folder,'data_cna.txt'),  sep='\t', index=False)

    
    if oncokb:   
             
        
        df_table_temp=df_table
        df_table_temp.rename({"discrete":"Copy_Number_Alteration"},inplace=True,axis=1)

        df_table_filtered=df_table_temp[df_table_temp["Copy_Number_Alteration"].isin([-2,2])]
        df_table_filtered.to_csv(os.path.join(output_folder,"data_cna_hg19.seg.to_annotate.txt"),sep="\t",index=False)
        temp_cna=reshape_cna(tsv_input,os.path.join(output_folder,"data_cna_hg19.seg.to_annotate.txt"),cancer,output_folder)
        annotate_cna(temp_cna,output_folder,tsv_input,Filter)
        

    return sID_path

def tabella_to_dict(df):
    result = {}
    for index, row in df.iterrows():
        row_values = (row['chrom'], row['loc.start'], row['loc.end'], row['num.mark'], row['seg.mean'], row['gene'], row['discrete'])
        if row['ID'] not in result:
            result[row['ID']] = []
        result[row['ID']].append(row_values)
    return result


def annotate_cna(path_cna,output_folder,tsv_input,Filter):


    out=path_cna.replace(".txt","_annotate.txt")
    os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {path_cna}\
                        -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")
                
                
    cna=pd.read_csv(out,sep="\t",dtype={"Copy_Number_Alteration":int})
    ##
    
    if "k" in Filter:

        logger.info("Filtering CNV based on OncoKB Annotation")
        cna=filter_OncoKB(cna)

    input_file=pd.read_csv(tsv_input,sep="\t")
    if "TC" in input_file.columns:
        for index, row in cna.iterrows():
            tc= int(input_file[input_file["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]]["TC"])
            cna.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc

        cna["Copy_Number_Alteration"]=0
        cna.loc[(cna["discrete"]>3)&(cna["discrete"]<5), "Copy_Number_Alteration"]=1
        cna.loc[cna["discrete"]>5, "Copy_Number_Alteration"]=2
        cna.loc[(cna["discrete"]>0)&(cna["discrete"]<0.8), "Copy_Number_Alteration"]=-1
        cna.loc[cna["discrete"]<=0, "Copy_Number_Alteration"]=-2

    data_cna=cna.pivot_table(index="Hugo_Symbol",columns="Tumor_Sample_Barcode",values="Copy_Number_Alteration",fill_value=0)
    data_cna.to_csv(os.path.join(output_folder,"data_cna.txt"),index=True,sep="\t")
    

def reshape_cna(input,cna_df_path,cancer,output_dir):
   
    input_file=pd.read_csv(input,sep="\t")    

    cna_df=pd.read_csv(cna_df_path,sep="\t")
    
    cna_df.rename({"ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"},inplace=True,axis=1)
    input_file.rename({"SampleID":"Tumor_Sample_Barcode"}, inplace=True,axis=1)

    if not "ONCOTREE_CODE" in input_file.columns:
        input_file["ONCOTREE_CODE"]=cancer
    
    annotate= pd.merge(cna_df[["Tumor_Sample_Barcode","Hugo_Symbol","seg.mean","Copy_Number_Alteration"]] ,input_file[["Tumor_Sample_Barcode","ONCOTREE_CODE"]],on="Tumor_Sample_Barcode")
    annotate.to_csv(os.path.join(output_dir,"temp_cna.txt"),index=False,sep="\t")
    
    return os.path.join(output_dir,"temp_cna.txt")


#  GET SNV

def check_snv_vcf(file,inputFolderSNV,multivcf):
    vcf=pd.read_csv(os.path.join(inputFolderSNV,file),comment="#",sep="\t",header=None)
    nsample=os.popen(f'bcftools query -l {os.path.join(inputFolderSNV,file)}').read().split("\n")
    nsample=[sample for sample in nsample if not sample ==""]
    if len(nsample)>1 and not multivcf:
        logger.critical("VCF contains multiple samples")
        exit(1)
     
    if all(vcf.apply(lambda row: row.astype(str).str.contains('GT').any(), axis=1)):
        return True
    else:
        return False

def get_snv_from_folder(inputFolderSNV,multivcf):
    files= os.listdir(inputFolderSNV)
    
    snv_vcf_files=[file for file in files if file.endswith("vcf")]
    try:
        check=list(map(lambda x: check_snv_vcf(x,inputFolderSNV,multivcf),snv_vcf_files))
        incorrect_files=[]
        for i, check_res in enumerate(check):
            if not check_res:
                incorrect_files.append(snv_vcf_files[i])
        if len(incorrect_files)!=0:
            logger.critical(f"It seems that the files \n{incorrect_files} \nare not SNV! Please check your SNV input data and try again.")
            exit()
    except:
        logger.warning("Could not check VCF type")
    logger.info(f"#{len(snv_vcf_files)} vcf files found in SNV folder")
    return snv_vcf_files

def snv_type_from_folder(input,tsv_input,snv_vcf_files):
    c = 0
    sID_path = dict()
    for case_folder in snv_vcf_files:
        try:
            snv_vcf= case_folder
            tsv_input_file=pd.read_csv(tsv_input,sep="\t")
            sampleID = tsv_input_file[tsv_input_file['SampleID'].apply(lambda x: x in snv_vcf)]["SampleID"].values.tolist()[0]
            if sampleID in sID_path:
                dup = open('sampleID_dup'.log, 'w')
                dup.write(sampleID+'\t'+'snv_vcf')
                dup.close()
            else:
                sID_path[sampleID] = os.path.join(input,snv_vcf)
        except Exception:
            log_noparsed = open('noParsed_snv.log', 'a')
            log_noparsed.write('[WARNING]'+case_folder+'\n')
            log_noparsed.close()
        c = c +1
    return sID_path



def extract_multiple_snv(multiple_vcf,input_dir):    

    if not os.path.exists(os.path.join(input_dir,"single_sample_vcf")):
        os.mkdir(os.path.join(input_dir,"single_sample_vcf"))
    os.system("bcftools query -l "+multiple_vcf+ " > "+input_dir+"/sample_id.txt")
    os.system("while read sample; do vcf-subset --exclude-ref -c $sample " + multiple_vcf +" > "+input_dir+"/single_sample_vcf/${sample}.vcf ; done <"+ input_dir+"/sample_id.txt")

    
def extract_multiple_cnv(multiple_vcf,input_dir):    

    if not os.path.exists(os.path.join(input_dir,"single_sample_vcf")):
        os.mkdir(os.path.join(input_dir,"single_sample_vcf"))
    os.system("bcftools query -l "+multiple_vcf+ " > "+input_dir+"/sample_id.txt")
    os.system("while read sample; do bcftools view -s $sample " + multiple_vcf +" > "+input_dir+"/single_sample_vcf/${sample}.vcf ; done <"+ input_dir+"/sample_id.txt")

    

# SCRATCH
def create_random_name_folder():
    nome_cartella = ''.join(random.choices(string.ascii_lowercase + string.digits, k=10))
    temporary = os.path.join(TMP, nome_cartella)
    try:
        os.mkdir(temporary)
    except FileNotFoundError:
        logger.critical(f"Scratch folder '{TMP}' not found! Check TMP field in conf.ini")
        exit()
    except Exception:
        logger.critical("Something went wrong while creating the vep tmp folder")
        exit()
    return(temporary)

def clear_scratch():
    for root, dirs, files in os.walk(TMP):
        for dir in dirs:
            shutil.rmtree(os.path.join(root,dir))
        

def vcf_filtering(sID_path,output_folder):
    sID_path_filtered = dict()
    for k, v in sID_path.items():
        root, vcf_file = os.path.split(v)
        
        out_filt=os.path.join(output_folder,OUTPUT_FILTERED) #TEST
        vcf_filtered = os.path.join(out_filt, vcf_file.replace('.vcf','')+'.FILTERED.vcf')
        logger.info(f'[FILTERING] {v}')
        vcf_filter.main(v, vcf_filtered)
        logger.info(f'filtered file {vcf_filtered} created!')
        sID_path_filtered[k] = vcf_filtered
    return sID_path_filtered


def vcf2maf_constructor(k, v, temporary,output_folder):
    cl = ['perl']
    cl.append(VCF2MAF)
    cl.append('--input-vcf')
    cl.append(v)
    root, file_vcf = os.path.split(v)
    out_file = os.path.join(output_folder,os.path.join(OUTPUT_MAF, file_vcf+'.maf'))
    cl.append('--output-maf')
    cl.append(out_file)
    if not CLINV =="":
        cl.append('--vep-custom')
        cl.append(CLINV)
    cl.append('--ref-fasta')
    cl.append(REF_FASTA)
    cl.append('--tmp-dir')
    cl.append(temporary)
    cl.append('--retain-fmt')
    cl.append('GT,GQ,AD,DP,VF')
    cl.append('--vep-path')
    cl.append(VEP_PATH)
    cl.append('--vep-data')
    cl.append(VEP_DATA)
    cl.append('--tumor-id')
    cl.append(k)
    return cl

def run_vcf2maf(cl):
    logger.info('Starting vcf2maf conversion...')
    logger.info(f'args={cl}')
    sout = subprocess.run(cl, capture_output=True)
 
    if sout.stderr!=None:
        if 'ERROR' not in sout.stderr.decode('ascii'):
            logger.warning(sout.stderr.decode('ascii').replace('ERROR: ',''))
        else:
            logger.error(sout.stderr.decode('ascii').replace('ERROR: ',''))
    
def get_table_from_folder(tsvpath):
    table_dict = dict()
    file=pd.read_csv(tsvpath,sep="\t",index_col=False, dtype=str)
    for _, row in file.iterrows():
        sampleID=str(row["SampleID"])
        if ".bam" in sampleID:
           sampleID=sampleID.replace(".bam","")
        if sampleID not in table_dict.keys():
            if "PatientID" in file.columns:
                table_dict[sampleID]=[str(row["PatientID"])]  
            else:
                table_dict[sampleID]=[str(row["SampleID"])]
    return table_dict
    
def flatten(nested_list):
    flat_list = []
    for sublist in nested_list:
        for item in sublist:
            flat_list.append(item)
    return flat_list

def write_clinical_patient(output_folder, table_dict):
    logger.info("Writing data_clinical_patient.txt file...")
    data_clin_samp = os.path.join(output_folder,'data_clinical_patient.txt')
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
        cil_sample.write(v+"\tNan\tNa\n")
    cil_sample.close()

def write_clinical_sample(output_folder, table_dict):
    logger.info("Writing data_clinical_sample.txt file...")
    data_clin_samp = os.path.join(output_folder, 'data_clinical_sample.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tSample Identifier\tMSI\tTMB\tMSI_THR\tTMB_THR\n')
    cil_sample.write('#Patient identifier\tSample Identifier\tMicro Satellite Instability\tMutational Tumor Burden\tMSI_THR\tTMB_THR\n')
    cil_sample.write('#STRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\n')
    cil_sample.write('#1\t1\t1\t1\t1\t1\n')
    cil_sample.write('PATIENT_ID\tSAMPLE_ID\tMSI\tTMB\tMSI_THR\tTMB_THR\n')
    
    for k, v in table_dict.items():
        cil_sample.write(str(v[0])+'\t'+k+'\t'+str(v[1])+'\t'+str(v[3])+'\t'+str(v[2])+'\t'+str(v[4])+'\n')
    cil_sample.close()
    
def write_clinical_sample_only(output_folder, table_dict):
    logger.info("Writing data_clinical_sample.txt file...")
    data_clin_samp = os.path.join(output_folder, 'data_clinical_sample.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tSample Identifier\n')
    cil_sample.write('#Patient identifier\tSample Identifier\n')
    cil_sample.write('#STRING\tSTRING\n')
    cil_sample.write('#1\t1\n')
    cil_sample.write('PATIENT_ID\tSAMPLE_ID\n')
    for k, v in table_dict.items():
        cil_sample.write(v[0]+'\t'+k+'\n')
    cil_sample.close()


def walk_folder(input, multiplevcf,oncokb, cancer,output_folder, type=None,Filter=""):

    logger.info("Starting walk_folder script:")
    logger.info(f"walk_folder args [input:{input}, output_folder:{output_folder},type:{type},Filter:{Filter}]")
    config.read('conf.ini')


    try:
        tsvfile=[file for file in os.listdir(input) if file.endswith("tsv")][0]
    except IndexError:
        logger.critical(f"It seems that no tsv file is in your folder!")
        exit()
    tsvpath=os.path.join(input,tsvfile)    
    inputFolderSNV=os.path.abspath(os.path.join(input,"SNV"))
    inputFolderCNV=os.path.abspath(os.path.join(input,"CNV"))
    inputFolderFusions=os.path.abspath(os.path.join(input,"FUSIONS"))
   
    ###############################
    ###       OUTPUT FOLDER     ###
    ###############################
 
    output_folder=create_folder(output_folder,type)  
 

    if os.path.exists(inputFolderCNV) and not type in ["snv","fus","tab"]:
        if multiplevcf:
            multivcf = os.listdir(inputFolderCNV)[0]
            extract_multiple_cnv(os.path.join(inputFolderCNV,multivcf),inputFolderCNV)
            inputFolderCNV= os.path.join(inputFolderCNV,"single_sample_vcf")
        
        logger.info("Check CNV files...")
        
        case_folder_arr_cnv = get_cnv_from_folder(inputFolderCNV,multiplevcf)
        
        logger.info("Everything ok!")

    if os.path.exists(inputFolderSNV) and not type in ["cnv","fus","tab"]:
        
        if multiplevcf:
            multivcf = os.listdir(inputFolderSNV)[0]
            extract_multiple_snv(os.path.join(inputFolderSNV,multivcf),inputFolderSNV)
            inputFolderSNV= os.path.join(inputFolderSNV,"single_sample_vcf")
        
       
        logger.info("Check SNV files...")
        case_folder_arr = get_snv_from_folder(inputFolderSNV,multiplevcf)
        logger.info("Everything ok!")

    ###############################
    ###       SNV AND CNV       ###
    ###############################
 
    if not type in ["snv","fus","tab"] and os.path.exists(inputFolderCNV) :
        sID_path_cnv = cnv_type_from_folder(inputFolderCNV,tsvpath,case_folder_arr_cnv,oncokb,cancer,output_folder,Filter)
    
    if not type in ["cnv","fus","tab"] and os.path.exists(inputFolderSNV):
        sID_path_snv = snv_type_from_folder(inputFolderSNV,tsvpath,case_folder_arr)
        
        logger.info("Starting vcf2maf conversion...")
       
        temporary = create_random_name_folder()
        sID_path_filtered = vcf_filtering(sID_path_snv,output_folder)
        for k, v in sID_path_filtered.items():
           
            cl = vcf2maf_constructor(k, v, temporary,output_folder)
            run_vcf2maf(cl)
    
    
    logger.info("Clearing scratch folder...")
    clear_scratch()
    
    ###############################
    ###       GET FUSION        ###
    ###############################
    
    if os.path.exists(inputFolderFusions) and not type in ["cnv","snv","tab"]:
        
        fusion_table_file = os.path.join(output_folder,'data_sv.txt')
        tsv_file=[file for file in os.listdir(input) if file.endswith("tsv")][0]
        tsvpath=os.path.join(input,tsv_file)
        patient_file=pd.read_csv(tsvpath,sep="\t",index_col=False, dtype=str)
        if len(os.listdir(os.path.join(input,"FUSIONS")))==0 :
            logger.warning("No fusions found")
        else:    
            fusion_files=[file for file in os.listdir(os.path.join(input,"FUSIONS")) if "tsv" in file]
            
            nfusion=len(fusion_files) 
            logger.info(f"Found {nfusion} Fusion files ")

            for fusion_file in fusion_files :
                
                ff = pd.read_csv(os.path.join(input,"FUSIONS",fusion_file),sep="\t")
                
                if not set(["Sample_Id","SV_Status","Site1_Hugo_Symbol","Site2_Hugo_Symbol"]).issubset(ff.columns):
                    logger.warning(f"{fusion_file} does not contain required columns")
                    continue
            
                if ff.shape[0]==0:
                    logger.info(f"No Fusions found in {fusion_file}")
                    continue
                else:
                    logger.info(f"Fusions found in {fusion_file}")
                if not os.path.exists(fusion_table_file):
                    logger.info(f"Creating data_sv.txt file...")
                    fusion_table = open(fusion_table_file, 'w')
                    header = 'Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\n'
                    fusion_table.write(header)
                else:
                    fusion_table = open(fusion_table_file, 'a')
                for _, fus in ff.iterrows():
                    
                    if fus["Sample_Id"] in patient_file["SampleID"].values:
                        
                        fusion_table.write(str(fus["Sample_Id"])+'\t'+fus["SV_Status"]+'\tFUSION'+'\t'+str(fus['Site1_Hugo_Symbol'])+'\t'+str(fus['Site2_Hugo_Symbol'])+'\n') 
            fusion_table.close()
    
        if oncokb and  os.path.exists(fusion_table_file):
            
            data_sv=pd.read_csv(fusion_table_file,sep="\t")
            
            tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
            input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
            
            if "ONCOTREE_CODE" in input_file.columns:
            
                input_file["SampleID"]=input_file["SampleID"]
                fusion_table_df=data_sv.merge(input_file,how="inner",left_on="Sample_Id",right_on="SampleID")
                fusion_table_df.to_csv(fusion_table_file,sep="\t",index=False)
                fusion_table_file_out=fusion_table_file.replace(".txt","ann.txt")
                os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                        -o {fusion_table_file_out} -b {config.get('OncoKB', 'ONCOKB')}")
                
            else:   
                    
                fusion_table_file_out=fusion_table_file.replace(".txt","ann.txt")       
                os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                            -o {fusion_table_file_out} -t {cancer.upper()}  -b {config.get('OncoKB', 'ONCOKB')}")

            if "k" in Filter:
                fus_file=pd.read_csv(fusion_table_file_out,sep="\t")
                fus_file=filter_OncoKB(fus_file)
                fus_file.to_csv(fusion_table_file_out,index=False,sep="\t")
            os.system(f"mv {fusion_table_file_out}  {fusion_table_file}")     
        

    ###############################
    ###       MAKES TABLE       ###
    ###############################
    
    if not type in ["snv","cnv","fus"]:
        tsv_file=[file for file in os.listdir(input) if file.endswith("tsv")][0]
        tsvpath=os.path.join(input,tsv_file)    
        table_dict_patient = get_table_from_folder(tsvpath)
        logger.info("Writing clinical files...")
        write_clinical_patient(output_folder, table_dict_patient)

        tmb=True
        if not os.path.exists(os.path.join(input,"TMB")) or len(os.listdir(os.path.join(input,"TMB")))==0  :
            logger.warning("Data clinical sample will be written without TMB information")
            tmb=False
        else:
            tmb_file=pd.read_csv(os.path.join(input,"TMB","TMB.tsv"),sep="\t")

        msi=True   
        if not os.path.exists(os.path.join(input,"MSI")) or len(os.listdir(os.path.join(input,"MSI")))==0  :
            logger.warning("Data clinical sample will be written without MSI information")
            msi=False
        else:
            msi_file=pd.read_csv(os.path.join(input,"MSI","MSI.tsv"),sep="\t")

        if any([tmb,msi]):
            
            MSI_THR=config.get('MSI', 'THRESHOLD')
            TMB=ast.literal_eval(config.get('TMB', 'THRESHOLD'))

            for k, v in table_dict_patient.items():
                logger.info(f"Reading Tumor clinical parameters info for {v}...")
                
                if msi:
                    msi_value=msi_file[msi_file["SampleID"]==k]["MSI"]
                    if msi_value.empty:
                        table_dict_patient[k].append("NA")
                    else:
                        table_dict_patient[k].append(msi_file[msi_file["SampleID"]==k]["MSI"].values[0])

                    ## MSI Threshold
                    try:
                        if float(msi_file[msi_file["SampleID"]==k]["MSI"].values[0]) < float(MSI_THR):
                            table_dict_patient[k].append("Stable")   
                        else:
                            table_dict_patient[k].append('Unstable')
                    except:
                        table_dict_patient[k].append('NI')
                else:
                    table_dict_patient[k].append('NA')
                    table_dict_patient[k].append('NI')


                if tmb:
                    tmb_value=tmb_file[tmb_file["SampleID"]==k]["TMB"]
                    if tmb_value.empty:
                        table_dict_patient[k].append("NA")
                    else:
                        table_dict_patient[k].append(tmb_file[tmb_file["SampleID"]==k]["TMB"].values[0])

                    try:
                        found = False
                        for _k, _v in TMB.items():
                            if float(tmb_file[tmb_file["SampleID"]==k]["TMB"].values[0])<=float(_v):
                                    table_dict_patient[k].append(_k)
                                    found=True
                                    break
                            else:
                                table_dict_patient[k].append(list(TMB.keys())[-1])
                        if found==False:
                            table_dict_patient[k].append(list(TMB.keys())[-1])
                    except:
                        table_dict_patient[k].append('NI')

                else:
                    table_dict_patient[k].append('NA')
                    table_dict_patient[k].append('NI')

            write_clinical_sample(output_folder, table_dict_patient)
        else:
            write_clinical_sample_only(output_folder, table_dict_patient)

        logger.success("Walk script completed!\n")
    return output_folder
