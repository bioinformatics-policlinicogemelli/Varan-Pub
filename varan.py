import os 
import sys
import argparse
from loguru import logger 
import shutil

@logger.catch()
def varan(args):
    logger.info(f"Varan args [input:{args.input}, output_folder:{args.output_folder},  cancer:{args.Cancer}, \
                            type:{args.type}], \
                            update:{args.Update}, extract:{args.Extract}, remove:{args.Remove}")
    

    if not any([args.Update ,args.Extract , args.Remove]) :       

        if not args.Metafile:
            if not args.FilterOnly:

                ###########################
                #        1.  WALK         #
                ###########################
                from walk import walk_folder
                logger.info("Starting preparation study folder")
                if not os.path.exists("./scratch"):
                    logger.info("Creating scratch dir")
                    os.mkdir("./scratch")
                output_folder=walk_folder(args.input,args.multiple,args.oncoKB, args.Cancer, args.output_folder, args.type,args.Filter)
          
            ###########################
            #       2. FILTER         #
            ###########################
            from filter_clinvar import filter_main
            try:
                output_folder
            except NameError:
                output_folder=args.input
            if args.type =="snv" or (args.type==None and os.path.exists(os.path.join(args.input,"SNV"))):
                logger.info("Starting filter") 
                filter_main(args.input,output_folder, output_folder,args.Filter, args.oncoKB,args.Cancer) 

            ############################
            #      3. CONCATENATE      #
            ############################
            from concatenate import concatenate_main
            if args.type =="snv" or (args.type==None and os.path.exists(os.path.join(args.input,"SNV"))):
                logger.info("Concatenate mutation file")
                folders=[]
                if os.path.exists(os.path.join(output_folder,"MAF_Filtered")):
                    folders=["MAF_Filtered"]
                elif os.path.exists(os.path.join(output_folder,"MAF_OncoKB")):
                    folders=["MAF_OncoKB"]
                elif os.path.exists(os.path.join(output_folder,"maf")):
                    folders=["maf"]
                else:
                    pass
            
                for folder in folders:
                    input_folder=os.path.join(output_folder,folder)
                    output_file=os.path.join(input_folder,"data_mutations_extended.txt")
                    concatenate_main(input_folder,"maf",output_file)

                
                if not args.Filter==None:
                    logger.info("Extracting data_mutations_extended from MAF_Filtered folder") 
                    os.system("cp "+os.path.join(output_folder,os.path.join("MAF_Filtered","data_mutations_extended.txt"))+" "+ output_folder )
                elif os.path.exists(os.path.join(output_folder,"MAF_OncoKB")):
                    logger.info("Extracting data_mutations_extended from MAF_OncoKB folder") 
                    os.system("cp "+os.path.join(output_folder,os.path.join("MAF_OncoKB","data_mutations_extended.txt"))+" "+ output_folder )
                else:
                    logger.info("Extracting data_mutations_extended from maf folder") 
                    os.system("cp "+os.path.join(output_folder,os.path.join("maf","data_mutations_extended.txt"))+" "+ output_folder )
                    
            
            
                
        ###########################################
        #      4. MAKE AND POPULATE TABLES        #
        ###########################################
        from Make_meta_and_cases import meta_case_main
        logger.info("It's time to create tables!")
        try:
            output_folder
        except NameError:
            output_folder=args.input
        meta_case_main(args.Cancer,output_folder)

            
        ############################
        #      5. VALIDATION       #
        ############################
        from ValidateFolder import validateFolderlog
        logger.info("Starting Validation Folder")
        validateFolderlog(output_folder)
        logger.success("The end! The study is ready to be uploaded on cBioportal")
        if os.path.exists(os.path.join(output_folder,"snv_filtered")):
            shutil.make_archive(os.path.join(output_folder,"snv_filtered"),"zip",os.path.join(output_folder,"snv_filtered"))
            shutil.rmtree(os.path.join(output_folder,"snv_filtered"))
        if os.path.exists(os.path.join(output_folder,"maf")):     
            shutil.make_archive(os.path.join(output_folder,"maf"),"zip",os.path.join(output_folder,"maf"))
            shutil.rmtree(os.path.join(output_folder,"maf"))
        
    ############################
    #         UPDATE           #
    ############################

    elif args.Update: 
        from Update_script import update_main 
        logger.info("Starting Update study")
        update_main(args.Path,args.NewPath,args.output_folder)

    ############################
    #         DELETE           #
    ############################

    elif args.Remove:
        from Delete_script import delete_main 
        logger.info("Starting Delete samples from study")  
        delete_main(args.Path,args.SampleList,args.output_folder) #args.overWrite

    ############################
    #         EXTRACT          #
    ############################

    elif args.Extract:
        from ExtractSamples_script import extract_main
        logger.info("Starting Extract samples from study")
        extract_main(args.Path,args.SampleList,args.output_folder) #args.overWrite

    
#################################################################################################################

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)


@logger.catch()
def main(): 
    logger.remove()
    logfile="Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
    logger.level("INFO", color="<green>")
    logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True, catch=True, backtrace=True, diagnose=True)
    logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w", backtrace=True, diagnose=True)
    logger.info("Welcome to VARAN")

    parser = MyArgumentParser(add_help=True, exit_on_error=True, usage=None, description='Argument of Varan script')
    

    # WALK BLOCK
    parser.add_argument('-c', '--Cancer', required=False,help='Cancer Name')
    parser.add_argument('-i', '--input', required=False, help='input folder tsv with data or tsv with path of data')
    parser.add_argument('-m', '--multiple', required=False,action='store_true',help='Multiple sample VCF?')
    parser.add_argument('-t', '--type', required=False, choices=['snv', 'cnv','fus','tab'],help='Select the type of file to parse')
    
    # ANNOTATION BLOCK
    parser.add_argument('-k', '--oncoKB', required=False,action='store_true',help='OncoKB annotation')

    # FILTER BLOCK
    parser.add_argument('-f', '--Filter', required=False, help='Select filter for SNV',default="")

    # p -> filter==PASS , b-> Benign , v-> vaf, o-> Oncokb , g -> gnomAD, q > Consequence, y-> polyphen
    # s -> clin_sig
    
    # UPDATE BLOCK
    parser.add_argument('-u', '--Update', required=False,action='store_true',help='Add this argument if you want to concatenate two studies')
    parser.add_argument('-n', '--NewPath', required=False,help='Path of new study folder to add')  
    
    # DELETE BLOCK
    parser.add_argument('-r', '--Remove', required=False,action='store_true',help='Add this argument if you want to remove samples from a study')
    # EXTRACT BLOCK
    parser.add_argument('-e', '--Extract', required=False,action='store_true', help='Add this argument if you want to extract samples from a study')
    
    # COMMON BLOCK 
    parser.add_argument('-o', '--output_folder', required=False,help='Output folder')
    parser.add_argument('-l', '--sampleList', required=False,help='Path of file with list of SampleIDs') # List
    parser.add_argument('-p', '--Path', required=False,help='Path of original study folder')
    
    parser.add_argument('-meta', '--Metafile', required=False,action='store_true',help='Create MetaFile')
    parser.add_argument('-fo', '--FilterOnly', required=False,action='store_true',help='Use Filter script')
    
    try:
        args = parser.parse_args()  
    except ValueError :
        logger.critical("Error Argument: Output is required")
        exit(1)

    
    
    if args.Update and not all([args.Path,args.NewPath]):
        logger.critical("To update a study, you need to specify both original and new folder paths")
        exit(1)
    if args.Extract and not all([args.Path,args.SampleList]):
        logger.critical("To extract samples from a study, you need to specify both original folder path and list samples")
        exit(1)
    if args.Remove and not all([args.Path,args.SampleList]):
        logger.critical("To remove samples from a study, you need to specify both original folder path and list samples")
        exit()
    
    if not any([args.Update ,args.Extract , args.Remove]) and args.Cancer==None:
            logger.critical("Error Argument: Cancer name is required")
            exit()  
    if not any([args.Update ,args.Extract , args.Remove]) and args.input==None:
            logger.critical("Error Argument: Input is required")
            exit()

    varan(args)

if __name__ == '__main__':
    main()