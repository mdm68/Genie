#!/usr/bin/env python

import pandas as pd
import synapseclient
import os
#import subprocess
import argparse
import getpass
import string
import numpy as np

def synapse_login():
    """
    This function logs into synapse for you if credentials are saved.  
    If not saved, then user is prompted username and password.
    """
    try:
        syn = synapseclient.login()
    except Exception as e:
        print("Please provide your synapse username/email and password (You will only be prompted once)")
        Username = raw_input("Username: ")
        Password = getpass.getpass()
        syn = synapseclient.login(email=Username, password=Password,rememberMe=True)
    return syn

#VALIDATING CLINICAL
def checkColExist(clinicalDF, key):
    """
    This function checks if the key exists as a header in the clinical dataframe
    
    :params clinicalDF:     Pandas clinical dataframe 
    :params key:            Expected header column name

    :returns:               An error message or an empty string
    """
    if clinicalDF.get(key) is None:
        error = "<p>clinical file must have %s column</p>" % key
    else:
        error = ""
    return(error)

#CHECKS IF THE MAPPING IS CORRECT
def checkMapping(clinicalDF, primaryName, secondaryName, mapping, required=False, fileType = "Patient"):
    """
    This function checks if the column exists then checks if the values in the column have the correct integer values
    
    :params primaryName:        Primary expected column name
    :params secondaryName:      Secondary expected column name
    :params mapping:            List of possible values

    :returns:                   A tuple warning, error
    """
    warning = ""
    error = ""
    if clinicalDF.get(primaryName) is not None:
        race = primaryName
    else:
        race = secondaryName
    checkCol = checkColExist(clinicalDF, race)
    if checkCol != "":
        if required:
            error = "%s: clinical file must have %s column.\n" % (fileType,primaryName)
        else:
            warning = "%s: clinical file doesn't have %s column. A blank column will be added\n" % (fileType,primaryName)
    else:
        if not all([i in mapping for i in clinicalDF[race]]):
            error = "%s: Please double check your %s column.  This column must be these values %s and blank.\n" % (fileType, primaryName,", ".join(map(str,mapping)))
    return(warning, error)

def validateClinical(clinicalFilePath,oncotree_mapping,clinicalSamplePath=None):
    """
    This function validates the clinical file to make sure it adhere to the clinical SOP.
    
    :params clinicalFilePath:              Flattened clinical file or patient clinical file
    :params clinicalSamplePath:            Sample clinical file if patient clinical file is provided

    :returns:                              Error message
    """

    message="Below are some of the issues with your files:\n"
    clinicalDF = pd.read_csv(clinicalFilePath,sep="\t",comment="#")
    clinicalDF.columns = [col.upper() for col in clinicalDF.columns]
    clinicalDF = clinicalDF.fillna("")
    total_error = ""
    warning = ""
    if clinicalSamplePath is None:
        clinicalSampleDF = clinicalDF.copy()
    else:
        clinicalSampleDF = pd.read_csv(clinicalSamplePath,sep="\t",comment="#")
        clinicalSampleDF = clinicalSampleDF.fillna("")
    
    #CHECK: SAMPLE_ID
    if clinicalSampleDF.get("SAMPLE_ID") is not None:
        sampleId = 'SAMPLE_ID'
    else:
        sampleId = 'GENIE_SAMPLE_ID'
    error = checkColExist(clinicalSampleDF, sampleId)
    if error != "":
        total_error = total_error + "Sample: clinical file must have SAMPLE_ID column.\n"

    #CHECK: AGE_AT_SEQ_REPORT
    if clinicalSampleDF.get("AGE_AT_SEQ_REPORT") is not None:
        age = "AGE_AT_SEQ_REPORT"
    else:
        age = "AGE_SEQ_REPORT_DAYS"
    error = checkColExist(clinicalSampleDF, age)
    if error == "":
        #Deal with HIPAA converted rows from DFCI
        #First for loop can't int(text) because there are instances that have <3435 
        clinicalSampleDF[age] = [text.replace(">","") if isinstance(text, str) else text for text in clinicalSampleDF[age]]
        clinicalSampleDF[age] = [int(text.replace("<","")) if isinstance(text, str) and text != "" else text for text in clinicalSampleDF[age]]
        if not all([isinstance(i, (int,float)) or i == "" for i in clinicalSampleDF[age]]) or np.median(clinicalSampleDF[age]) < 100:
            total_error = total_error + "Sample: Please double check your AGE_AT_SEQ_REPORT.  This is the interval in DAYS (integer) between the patient's date of birth and the date of the sequencing report that is associated with the sample.\n"
    else:
        total_error = total_error + "Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"

    #CHECK: ONCOTREE_CODE
    error = checkColExist(clinicalSampleDF, "ONCOTREE_CODE")
    if error == "":
        if not all(clinicalSampleDF['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])):
            unmapped_oncotrees = clinicalSampleDF['ONCOTREE_CODE'][~clinicalSampleDF['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
            total_error = total_error + "Sample: Please double check that all your ONCOTREE CODES exist in the mapping. You have %d samples that don't map. These are the codes that don't map: %s\n" % (len(unmapped_oncotrees),",".join(set(unmapped_oncotrees)))
    else:
        total_error = total_error + "Sample: clinical file must have ONCOTREE_CODE column.\n"
    
    #CHECK: SAMPLE_TYPE
    error = checkColExist(clinicalSampleDF, "SAMPLE_TYPE")
    if error == "":
        if not all([isinstance(i, (int,float)) or i in [1,2,3,4,5,6,7] for i in clinicalSampleDF['SAMPLE_TYPE']]):
            total_error = total_error+"Sample: Please double check your SAMPLE_TYPE column. This column must be integers 1-7.\n"
    else:
        total_error = total_error + "Sample: clinical file must have SAMPLE_TYPE column.\n"

    #CHECK: SEQ_ASSAY_ID
    error = checkColExist(clinicalSampleDF, "SEQ_ASSAY_ID")
    if error == "":
        if not all([i != "" for i in clinicalSampleDF['SEQ_ASSAY_ID']]):
            warning = warning + "Sample: Please double check your SEQ_ASSAY_ID columns, there are empty rows.\n"
    else:
        total_error = total_error + "Sample: clinical file must have SEQ_ASSAY_ID column.\n"

    #CHECK: BIRTH_YEAR
    if clinicalDF.get("BIRTH_YEAR") is not None:
        birth_year = "BIRTH_YEAR"
    else:
        birth_year = "PATIENT_BIRTH_YEAR"
    error = checkColExist(clinicalDF, birth_year)
    if error == "": 
        #Deal with HIPAA converted rows from DFCI
        #First for loop can't int(text) because there are instances that have <3435 
        clinicalDF[birth_year] = [text.replace(">","") if isinstance(text, str) else text for text in clinicalDF[birth_year]]
        clinicalDF[birth_year] = [int(text.replace("<","")) if isinstance(text, str) and text != "" else text for text in clinicalDF[birth_year]]
        if not all([isinstance(i, (int,float)) or i == "" for i in clinicalDF[birth_year]]):
            total_error = total_error + "Patient: Please double check your BIRTH_YEAR column.  This column must be integers or blank.\n"
    else:
        total_error = total_error + "Patient: clinical file must have BIRTH_YEAR column.\n"

    #CHECK: PATIENT_ID
    if clinicalDF.get("PATIENT_ID") is not None:
        patientId = "PATIENT_ID"
    else:
        patientId = "GENIE_PATIENT_ID"
    error = checkColExist(clinicalDF, patientId)
    if error != "":
        total_error = total_error + "Patient: clinical file must have PATIENT_ID column.\n"

    # Create patient Id in sample data
    if clinicalSampleDF.get(patientId) is None:
        clinicalSampleDF[patientId] = ["-".join(samp.split("-")[0:3]) for samp in clinicalSampleDF[sampleId]]
    if clinicalSamplePath is not None:
        clinicalDF = clinicalDF.merge(clinicalSampleDF, on=patientId,how="outer")
    #CHECK: All patients must have associated sample data 
    if not all(clinicalSampleDF[patientId].isin(clinicalDF[patientId])):
        total_error = total_error + "Sample: All patients must have associated sample information\n"
    #CHECK: All samples must have associated patient data 
    if sum(clinicalDF[patientId].isnull()) >0:
        total_error = total_error + "Patient: All samples must have associated patient information\n"

    #CHECK: PRIMARY_RACE
    warn, error = checkMapping(clinicalDF,"PRIMARY_RACE","NAACCR_RACE_CODE_PRIMARY",[1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,20,21,22,25,26,27,28,30,31,32,88,96,97,98,99,""])
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: SECONDARY_RACE
    warn, error = checkMapping(clinicalDF,"SECONDARY_RACE","NAACCR_RACE_CODE_SECONDARY",[1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,20,21,22,25,26,27,28,30,31,32,88,96,97,98,99,""])
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: TERTIARY_RACE
    warn, error = checkMapping(clinicalDF,"TERTIARY_RACE","NAACCR_RACE_CODE_TERTIARY",[1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,20,21,22,25,26,27,28,30,31,32,88,96,97,98,99,""])
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: SEX
    warn, error = checkMapping(clinicalDF,"SEX","NAACCR_SEX_CODE",[1,2,3,4,5,6,9,""], required=True)
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: ETHNICITY
    warn, error = checkMapping(clinicalDF,"ETHNICITY","NAACCR_ETHNICITY_CODE",[0,1,2,3,4,5,6,7,8,9,""])
    warning = warning + warn
    total_error  = total_error + error 

    #Complete error message
    if total_error == "":
        message = "There is nothing wrong with your file!\n"
    else:
        message = message + total_error
    if warning != "":
        message = message + "-------------WARNINGS-------------\n" + warning
    return(message)

#VALIDATING MAF
def validateMAF(filePath):
    """
    This function validates the clinical file to make sure it adhere to the clinical SOP.
    
    :params filePath:     Path to mutation file

    :returns:             Text with all the errors in the clinical file
    """


    first_header = ['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Tumor_Sample_Barcode']
    correct_column_headers = ['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Tumor_Sample_Barcode','t_alt_count','t_depth']
    optional_headers = ['t_ref_count','n_depth','n_ref_count','n_alt_count']
    
    mutationDF = pd.read_csv(filePath,sep="\t",comment="#")

    total_error = ""
    message = "Below are some of the issues with your file:\n"
    warning = ""

    if mutationDF.columns[0] not in first_header:
        total_error = total_error+"First column header must be one of these (Case sensitive): %s.\n" % ", ".join(first_header)

    if not all(mutationDF.columns.isin(correct_column_headers)):
        total_error = total_error + "Your mutation file must at least have these headers (Case sensitive): %s.\n" % ",".join([i for i in correct_column_headers if i not in mutationDF.columns.values])

    if not all(mutationDF.columns.isin(optional_headers)):
        warning = warning + "Your mutation file does not have the column headers that can give extra information to the processed mutation file (Case sensitive): %s.\n" % ", ".join([i for i in optional_headers if i not in mutationDF.columns.values ])      

    if total_error == "":
        message = "There is nothing wrong with your file!\n"
    else:
        message = message + total_error
    message = message + warning
    return(message)

### VALIDATING VCF
def contains_whitespace(x):
    """
    Helper function for validateVCF.  No whitespace is allowed in VCF files
    """
    return(sum([" " in i for i in x if isinstance(i, str)]))

# Resolve missing read counts

def validateVCF(filePath):
    """
    This function validates the VCF file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to VCF file

    :returns:             Text with all the errors in the VCF file
    """
    message="Below are some of the issues with your files:\n"
    
    REQUIRED_HEADERS = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    #FORMAT is optional
    total_error = ""
    with open(filePath,"r") as foo:
        import csv
        temp = csv.reader(foo, delimiter="\t")
        head = 1
        while head == 1:
            headers = temp.next()
            head = len(headers)

    vcf = pd.read_csv(filePath, sep="\t",comment="#",header=None,names=headers)
    
    if sum(vcf.columns.isin(REQUIRED_HEADERS)) != 8:
        total_error = total_error + "Your vcf file must have these headers: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n"

    if len(vcf.columns) > 8:
        if "FORMAT" not in vcf.columns:
            total_error = total_error + "Your vcf file must have FORMAT header if genotype columns exist.\n"
   
    #Require that they report variants mapped to either GRCh37 or hg19 without 
    #the chr-prefix. variants on chrM are not supported
    nochr = ["chr" in i for i in vcf['#CHROM'] if isinstance(i, str)]
    if sum(nochr) > 0:
        total_error = total_error + "Your vcf file must not have the chr prefix in front of chromosomes.\n"
    if sum(vcf['#CHROM'].isin(["chrM"])) > 0:
        total_error = total_error + "Your vcf file must not have variants on chrM.\n"

    #No white spaces
    temp = vcf.apply(lambda x: contains_whitespace(x), axis=1)
    if sum(temp) >0:
        total_error = total_error + "Your vcf file must not have any white spaces in any of the columns.\n"
    #I can also recommend a `bcftools query` command that will parse a VCF in a detailed way, 
    #and output with warnings or errors if the format is not adhered too
    if total_error == "":
        message = "There is nothing wrong with your file!\n"
    else:
        message = message + total_error

    return(message)

#VALIDATING CNV
def validateCNV(filePath):
    """
    This function validates the CNV (linear or discrete) file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to CNV file

    :returns:             Text with all the errors in the CNV file
    """
    total_error = ""
    message="Below are some of the issues with your files:\n"

    cnvDF = pd.read_csv(filePath,sep="\t",comment="#")
    if cnvDF.columns[0] != "Hugo_Symbol":
        total_error = total_error + "Your cnv file's first column must be Hugo_Symbol (Case sensitive)\n"
    
    cnvDF.drop(cnvDF.columns[[0]], axis=1, inplace=True)
    if cnvDF.get("Entrez_Gene_Id") is not None:
        del cnvDF['Entrez_Gene_Id']

    if not all(~cnvDF.isnull()):
        total_error = total_error + "Your cnv file must not have any empty values\n"
    
    if not all(cnvDF.applymap(lambda x: isinstance(x, float))):
        total_error = total_error + "All values must be numerical values\n"

    if total_error == "":
        message = "There is nothing wrong with your file!\n"
    else:
        message = message + total_error

    return(message)

def validateFusion(filePath):
    """
    This function validates the Fusion file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to Fusion file

    :returns:             Text with all the errors in the Fusion file
    """
    total_error = ""
    message="Below are some of the issues with your files:\n"

    fusionDF = pd.read_csv(filePath,sep="\t",comment="#")

    REQUIRED_HEADERS = ['Hugo_Symbol','Entrez_Gene_Id','Center','Tumor_Sample_Barcode','Fusion','DNA_support','RNA_support','Method','Frame']

    if not all(fusionDF.columns.isin(REQUIRED_HEADERS)):
        total_error = total_error + "Your fusion file must at least have these headers (Case sensitive): %s.\n" % ",".join([i for i in REQUIRED_HEADERS if i not in mutationDF.columns.values])
    
    if total_error == "":
        message = "There is nothing wrong with your file!\n"
    else:
        message = message + total_error

    return(message)

#Validate SEG files
def validateSEG(filePath):
    """
    This function validates the SEG file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to SEG file

    :returns:             Text with all the errors in the SEG file
    """
    total_error = ""
    message="Below are some of the issues with your files:\n"

    segDF = pd.read_csv(filePath,sep="\t",comment="#")

    REQUIRED_HEADERS = ['ID','chromosome','loc.start','loc.end','num.mark','seg.mean']
    
    if not all(segDF.columns.isin(REQUIRED_HEADERS)):
        total_error = total_error + "Your fusion file must at least have these headers: %s.\n" % ",".join([i for i in REQUIRED_HEADERS if i not in mutationDF.columns.values])
    
    if total_error == "":
        message = "There is nothing wrong with your file!\n"
    else:
        message = message + total_error

    return(message)


VALIDATE_MAPPING = {'maf':validateMAF,
                    'clinical':validateClinical,
                    'vcf':validateVCF,
                    'cnv':validateCNV,
                    'fusion':validateFusion,
                    'seg':validateSEG}

def perform_main(args):
    """
    This performs the validation of files

    :returns:   Text with the errors of the chosen file
    """
    syn = synapse_login()
    validate_func = VALIDATE_MAPPING[args.fileType]
    if args.fileType == "clinical":
        oncotree_mapping_ent = syn.tableQuery('SELECT * FROM syn7437073')
        oncotree_mapping = oncotree_mapping_ent.asDataFrame()
        
        if len(args.file) > 1:
            for filename in args.file:
                if not os.path.isfile(filename):
                    raise ValueError("File doesn't exist")
                if "patient" not in filename.lower() and "sample" not in filename.lower():
                    raise ValueError("You can only specify two files when validating clinical files. 'patient' and 'sample' must exist in either filenames")

            if "patient" in args.file[0].lower():
                message = validate_func(args.file[0],oncotree_mapping,args.file[1])
            else:
                message = validate_func(args.file[1],oncotree_mapping,args.file[0])
        else:
            if not os.path.isfile(args.file[0]):
                raise ValueError("File doesn't exist")
            message = validate_func(args.file[0],oncotree_mapping)
    else:
        message = validate_func(args.file[0])
    print(message)
    return(message)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Validate GENIE files')

    parser.add_argument("fileType", type=str, choices = ['maf','clinical','fusion','cnv','vcf','seg'],
                        help='File type that you are validating: maf, clinical, fusion, cnv, vcf, seg')
    parser.add_argument("file", type=str, nargs="+",
                        help='File(s) that you are validating.  If you validation your clinical files and you have both sample and patient files, you must provide both')

    args = parser.parse_args()
    perform_main(args)
