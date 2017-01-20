#!/usr/bin/env python

import pandas as pd
import synapseclient
import os
#import subprocess
import argparse
import getpass
import string
import httplib2 as http
import json

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

#DEPENDENCIES
#PANDAS
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
    return(syn)

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
    
    :params clinicalDF          Patient/sample/flattened clinical file
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
        if not all([i in mapping.tolist() for i in clinicalDF[race]]):
            error = "%s: Please double check your %s column.  This column must be these values %sor blank.\n" % (fileType, primaryName,", ".join(map(str,mapping)).replace(".0",""))
    return(warning, error)

#Getting the GENIE mapping synapse tables
def getGenieMapping(syn, synId):
    """
    This function gets the GENIE mapping tables
    
    :params synId:          Synapse Id of synapse table

    :returns:               Table dataframe
    """
    table_ent = syn.tableQuery('SELECT * FROM %s' %synId)
    table = table_ent.asDataFrame()
    table = table.fillna("")
    return(table)

#Validate genes
def hgncRestCall(path):
    headers = {'Accept': 'application/json',}

    uri = 'http://rest.genenames.org'

    target = urlparse(uri+path)
    method = 'GET'
    body = ''
    h = http.Http()
    response, content = h.request(target.geturl(),
                                  method,
                                  body,
                                  headers)
    if response['status'] == '200':
        data = json.loads(content)
        if len(data['response']['docs']) == 0:
            return(False, None)
        else:
            return(True, data['response']['docs'][0]['symbol'])
    else:
        #return(False, response['status'])
        return(False, None)

# Validation of gene names
def validateSymbol(gene):
    #if gene name is "MLL2" OR "MLL4" AND "chrom is 12", then the HUGO symbol is KMT2D
    #if gene name is "MLL2" OR "MLL4" AND "chrom is 19", then the HUGO symbol is KMT2B
    path = '/fetch/symbol/%s' %  gene
    verified, symbol = hgncRestCall(path)
    if not verified:
        path = '/fetch/prev_symbol/%s' %  gene
        verified, symbol = hgncRestCall(path)
    if not verified:
        path = '/fetch/alias_symbol/%s' %  gene
        verified, symbol = hgncRestCall(path)       
    if gene == symbol:
        return(True)
    else:
        if symbol is None:
            print("%s cannot be mapped" % gene)
        else:
            print("%s is remapped to %s" % (gene, symbol))
        return({gene:symbol})

def validateClinical(clinicalFilePath,oncotree_mapping,sampleType_mapping,ethnicity_mapping,race_mapping,sex_mapping,clinicalSamplePath=None):
    """
    This function validates the clinical file to make sure it adhere to the clinical SOP.
    
    :params clinicalFilePath:              Flattened clinical file or patient clinical file
    :params clinicalSamplePath:            Sample clinical file if patient clinical file is provided

    :returns:                              Error message
    """
    clinicalDF = pd.read_csv(clinicalFilePath,sep="\t",comment="#")
    clinicalDF.columns = [col.upper() for col in clinicalDF.columns]
    clinicalDF = clinicalDF.fillna("")
    total_error = ""
    warning = ""
    if clinicalSamplePath is None:
        clinicalSampleDF = clinicalDF.copy()
    else:
        clinicalSampleDF = pd.read_csv(clinicalSamplePath,sep="\t",comment="#")
        clinicalSampleDF.columns = [col.upper() for col in clinicalSampleDF.columns]
        clinicalSampleDF = clinicalSampleDF.fillna("")
    
    #CHECK: SAMPLE_ID
    if clinicalSampleDF.get("SAMPLE_ID") is not None:
        sampleId = 'SAMPLE_ID'
    else:
        sampleId = 'GENIE_SAMPLE_ID'
    error = checkColExist(clinicalSampleDF, sampleId)
    if error != "":
        total_error = total_error + "Sample: clinical file must have SAMPLE_ID column.\n"
    else:
        if sum(clinicalSampleDF[sampleId].isnull()) > 0:
            total_error = total_error + "Sample: There can't be any blank values for PATIENT_ID\n"

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
        if not all([isinstance(i, (int,float)) or i == "" for i in clinicalSampleDF[age]]) or pd.np.median(clinicalSampleDF[age]) < 100:
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

        if not all(clinicalSampleDF['SAMPLE_TYPE'].isin(sampleType_mapping['CODE'])):
        #if not all([isinstance(i, (int,float)) or i in [1,2,3,4,5,6,7] for i in clinicalSampleDF['SAMPLE_TYPE']]):
            total_error = total_error+"Sample: Please double check your SAMPLE_TYPE column. This column must be %s.\n" % ", ".join(map(str,sampleType_mapping['CODE']))
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
    else:
        if sum(clinicalDF[patientId].isnull()) > 0:
            total_error = total_error + "Patient: There can't be any blank values for PATIENT_ID\n"

    # Create patient Id in sample data
    if clinicalSampleDF.get(patientId) is None:
        clinicalSampleDF[patientId] = ["-".join(samp.split("-")[0:3]) for samp in clinicalSampleDF[sampleId]]
    #CHECK: All samples must have associated patient data 
    if not all(clinicalSampleDF[patientId].isin(clinicalDF[patientId])):
        total_error = total_error + "Sample: All samples must have associated patient information\n"
    #CHECK: All patients must have associated sample data 
    if not all(clinicalDF[patientId].isin(clinicalSampleDF[patientId])):
        total_error = total_error + "Sample: All patients must have associated sample information\n"

    #CHECK: PRIMARY_RACE
    warn, error = checkMapping(clinicalDF,"PRIMARY_RACE","NAACCR_RACE_CODE_PRIMARY",race_mapping['CODE'])
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: SECONDARY_RACE
    warn, error = checkMapping(clinicalDF,"SECONDARY_RACE","NAACCR_RACE_CODE_SECONDARY",race_mapping['CODE'])
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: TERTIARY_RACE
    warn, error = checkMapping(clinicalDF,"TERTIARY_RACE","NAACCR_RACE_CODE_TERTIARY",race_mapping['CODE'])
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: SEX
    warn, error = checkMapping(clinicalDF,"SEX","NAACCR_SEX_CODE",sex_mapping['CODE'], required=True)
    warning = warning + warn
    total_error  = total_error + error 

    #CHECK: ETHNICITY
    warn, error = checkMapping(clinicalDF,"ETHNICITY","NAACCR_ETHNICITY_CODE",ethnicity_mapping['CODE'])
    warning = warning + warn
    total_error  = total_error + error

    return(total_error, warning)

#VALIDATING MAF
def validateMAF(filePath):
    """
    This function validates the clinical file to make sure it adhere to the clinical SOP.
    
    :params filePath:     Path to mutation file

    :returns:             Text with all the errors in the clinical file
    """

    first_header = ['CHROMOSOME','HUGO_SYMBOL','TUMOR_SAMPLE_BARCODE']
    correct_column_headers = ['CHROMOSOME','START_POSITION','REFERENCE_ALLELE','TUMOR_SAMPLE_BARCODE','T_ALT_COUNT','T_DEPTH'] #T_REF_COUNT + T_ALT_COUNT = T_DEPTH
    optional_headers = ['T_REF_COUNT','N_DEPTH','N_REF_COUNT','N_ALT_COUNT']
    tumors = ['TUMOR_SEQ_ALLELE2','TUMOR_SEQ_ALLELE1']
    
    mutationDF = pd.read_csv(filePath,sep="\t",comment="#",na_values = ['-1.#IND', '1.#QNAN', '1.#IND', 
                             '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A', '#NA', 'NULL', 'NaN', 
                             '-NaN', 'nan','-nan',''],keep_default_na=False)
    mutationDF.columns = [col.upper() for col in mutationDF.columns]

    total_error = ""
    warning = ""
    
    #CHECK: First column must be in the first_header list
    if mutationDF.columns[0] not in first_header:
        total_error = total_error+"First column header must be one of these: %s.\n" % ", ".join(first_header)
    
    #CHECK: Everything in correct_column_headers must be in mutation file
    if not all([i in mutationDF.columns for i in correct_column_headers]):
        total_error = total_error + "Your mutation file must at least have these headers: %s. If you are missing T_DEPTH, you must have T_REF_COUNT!\n" % ",".join([i for i in correct_column_headers if i not in mutationDF.columns.values])
    
    #CHECK: Must have either Tumor_Seq_Allele1 or 2
    tumor_cols = [i for i in tumors if i in mutationDF.columns.values]
    if len(tumor_cols) == 0:
        total_error = total_error + "Your mutation file must at least have one of these headers: %s." % ",".join(tumors)

    for i in tumor_cols:
        if sum(mutationDF[i] == "NA") > 0:
            warning = warning + "Your %s column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n" % i

    if sum(mutationDF['REFERENCE_ALLELE'] == "NA") > 0:
        warning = warning + "Your REFERENCE_ALLELE column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n"

    #CHECK: Mutation file would benefit from columns in optional_headers
    if not all([i in mutationDF.columns for i in optional_headers]):
        warning = warning + "Your mutation file does not have the column headers that can give extra information to the processed mutation file: %s.\n" % ", ".join([i for i in optional_headers if i not in mutationDF.columns.values ])      
    
    #CHECK: mutation file must not have empty reference or variant alleles
    if sum(mutationDF['REFERENCE_ALLELE'].isnull()) > 0:
        total_error = total_error + "Your mutation file cannot have any empty REFERENCE_ALLELE values.\n"
    if sum(mutationDF['TUMOR_SEQ_ALLELE2'].isnull()) >0:
        total_error = total_error + "Your mutation file cannot have any empty TUMOR_SEQ_ALLELE2 values.\n"

    return(total_error, warning)

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
    REQUIRED_HEADERS = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    #FORMAT is optional
    total_error = ""
    warning = ""
    headers = None
    with open(filePath,"r") as foo:
        for i in foo:
            if i.startswith("#CHROM"):
                headers = i.replace("\n","").split("\t")
    if headers is not None:
        vcf = pd.read_csv(filePath, sep="\t",comment="#",header=None,names=headers)
    else:
        raise ValueError("Your vcf must start with the header #CHROM")
    
    if sum(vcf.columns.isin(REQUIRED_HEADERS)) != 8:
        total_error = total_error + "Your vcf file must have these headers: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n"

    if len(vcf.columns) > 8:
        if "FORMAT" not in vcf.columns:
            total_error = total_error + "Your vcf file must have FORMAT header if genotype columns exist.\n"
   
    #Require that they report variants mapped to either GRCh37 or hg19 without 
    #the chr-prefix. variants on chrM are not supported
    nochr = ["chr" in i for i in vcf['#CHROM'] if isinstance(i, str)]
    if sum(nochr) > 0:
        warning = warning + "Your vcf file must not have the chr prefix in front of chromosomes.\n"
    if sum(vcf['#CHROM'].isin(["chrM"])) > 0:
        total_error = total_error + "Your vcf file must not have variants on chrM.\n"

    #No white spaces
    temp = vcf.apply(lambda x: contains_whitespace(x), axis=1)
    if sum(temp) >0:
        warning = warning + "Your vcf file must not have any white spaces in any of the columns.\n"
    #I can also recommend a `bcftools query` command that will parse a VCF in a detailed way, 
    #and output with warnings or errors if the format is not adhered too
    return(total_error, warning)

#VALIDATING CNV
def validateCNV(filePath):
    """
    This function validates the CNV (linear or discrete) file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to CNV file

    :returns:             Text with all the errors in the CNV file
    """
    total_error = ""
    warning = ""
    cnvDF = pd.read_csv(filePath,sep="\t",comment="#")
    cnvDF.columns = [col.upper() for col in cnvDF.columns]

    if cnvDF.columns[0] != "HUGO_SYMBOL":
        total_error = total_error + "Your cnv file's first column must be Hugo_symbol\n"
    
    cnvDF.drop(cnvDF.columns[[0]], axis=1, inplace=True)
    if cnvDF.get("ENTREZ_GENE_ID") is not None:
        del cnvDF['ENTREZ_GENE_ID']

    if not all(~cnvDF.isnull()):
        total_error = total_error + "Your cnv file must not have any empty values\n"
    
    if not all(cnvDF.applymap(lambda x: isinstance(x, float))):
        total_error = total_error + "All values must be numerical values\n"

    print("VALIDATING GENE SYMBOLS")
    invalidated_genes = cnv["HUGO_SYMBOL"].drop_duplicates().apply(validateSymbol)

    return(total_error, warning)

def validateFusion(filePath):
    """
    This function validates the Fusion file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to Fusion file

    :returns:             Text with all the errors in the Fusion file
    """
    total_error = ""
    warning = ""

    fusionDF = pd.read_csv(filePath,sep="\t",comment="#")
    fusionDF.columns = [col.upper() for col in fusionDF.columns]

    REQUIRED_HEADERS = ['HUGO_SYMBOL','ENTREZ_GENE_ID','CENTER','TUMOR_SAMPLE_BARCODE','FUSION','DNA_SUPPORT','RNA_SUPPORT','METHOD','FRAME']

    if not all([i in fusionDF.columns for i in REQUIRED_HEADERS]):
        total_error = total_error + "Your fusion file must at least have these headers: %s.\n" % ",".join([i for i in REQUIRED_HEADERS if i not in mutationDF.columns.values])
    print("VALIDATING GENE SYMBOLS")
    invalidated_genes = fusionDF["HUGO_SYMBOL"].drop_duplicates().apply(validateSymbol)

    return(total_error, warning)

#Validate SEG/CBS files
def validateSEG(filePath):
    """
    This function validates the SEG file to make sure it adhere to the genomic SOP.
    
    :params filePath:     Path to SEG file

    :returns:             Text with all the errors in the SEG file
    """
    total_error = ""
    warning = ""

    segDF = pd.read_csv(filePath,sep="\t",comment="#")
    segDF.columns = [col.upper() for col in segDF.columns]

    REQUIRED_HEADERS = ['ID','CHROM','LOC.START','LOC.END','NUM.MARK','SEG.MEAN']
    
    if not all(segDF.columns.isin(REQUIRED_HEADERS)):
        total_error = total_error + "Your fusion file must at least have these headers: %s.\n" % ",".join([i for i in REQUIRED_HEADERS if i not in mutationDF.columns.values])
    print("VALIDATING GENE SYMBOLS")   
    invalidated_genes = segDF["HUGO_SYMBOL"].drop_duplicates().apply(validateSymbol)

    return(total_error, warning)

def validateBED(filePath):
    total_error = ""
    warning = ""

    bed = pd.read_csv(filePath, sep="\t",comment="#")
    

    print("VALIDATING GENE SYMBOLS")   
    invalidated_genes = bed["HUGO_SYMBOL"].drop_duplicates().apply(validateSymbol)


def validateFileName(args):
    VALIDATE_FILENAME = {'maf':"data_mutations_extended_%s.txt",
                         'clinical': ["data_clinical_supp_%s.txt", "data_clinical_supp_sample_%s.txt", "data_clinical_supp_patient_%s.txt"],
                         'vcf':"GENIE-%s-",
                         'cnv':"data_CNA_%s.txt",
                         'fusion':"data_fusions_%s.txt",
                         'seg':"genie_data_cna_hg19_%s.seg"}

    assert all([os.path.isfile(filename) for filename in args.file]), "Files must exist on the drive"
    if args.fileType == "clinical":
        formatting = [i % args.center for i in VALIDATE_FILENAME[args.fileType]]
        if len(args.file) > 1:
            assert len(set(args.file)) > 1, "Must submit two different filenames!"
            assert sum([os.path.basename(i) in formatting[1:3] for i in args.file]) == 2, "When submitting a patient and sample file, these must be named: %s!" % ", ".join(formatting[1:3]) 
        else:
            assert "patient" not in args.file[0] and "sample" not in args.file[0], "If you submit a patient or sample file, you must submit both at the same time: eg. python validateGENIE.py clinical data_clinical_supp_patient_SAGE.txt data_clinical_supp_sample_SAGE.txt SAGE"
            assert os.path.basename(args.file[0]) == formatting[0], "Clinical file must be named: %s!" % formatting[0]
    else:
        formatting = VALIDATE_FILENAME[args.fileType] % args.center
        if args.fileType == "vcf":
            assert os.path.basename(args.file[0]).startswith(formatting), "VCF filename must be in this format: GENIE-%s-patientId-sampleId!" % args.center 
        else:
            assert os.path.basename(args.file[0]) == formatting, "%s filename must be: %s!" % (args.fileType, formatting)

def perform_main(args):
    """
    This performs the validation of files

    :returns:   Text with the errors of the chosen file
    """

    VALIDATE_MAPPING = {'maf':validateMAF,
                'clinical':validateClinical,
                'vcf':validateVCF,
                'cnv':validateCNV,
                'fusion':validateFusion,
                'seg':validateSEG,
                'bed':validateBED}
    syn = synapse_login()
    #CHECK: Fail if filename is incorrect
    try:
        validateFileName(args)
    except AssertionError as e:
        raise ValueError("Your filename is incorrect!\n%s\nPlease change your filename before you run the validator again."  % e)
    
    validate_func = VALIDATE_MAPPING[args.fileType]
    if args.fileType == "clinical":
        oncotree_mapping = getGenieMapping(syn, "syn7437073")
        sampleType_mapping = getGenieMapping(syn, "syn7434273")
        ethnicity_mapping = getGenieMapping(syn, "syn7434242")
        race_mapping = getGenieMapping(syn, "syn7434236")
        sex_mapping = getGenieMapping(syn, "syn7434222")

        if len(args.file) > 1:
            if "patient" in args.file[0].lower():
                total_error, warning = validate_func(args.file[0],oncotree_mapping,sampleType_mapping,ethnicity_mapping,race_mapping,sex_mapping,clinicalSamplePath=args.file[1])
            else:
                total_error, warning = validate_func(args.file[1],oncotree_mapping,sampleType_mapping,ethnicity_mapping,race_mapping,sex_mapping,clinicalSamplePath=args.file[0])
        else:
            total_error, warning = validate_func(args.file[0],oncotree_mapping,sampleType_mapping,ethnicity_mapping,race_mapping,sex_mapping)
    else:
        total_error, warning = validate_func(args.file[0])
    
    #Complete error message
    message = "Below are some of the issues with your file:\n"
    if total_error == "":
        message = "There is nothing wrong with the contents of your file!\n"
    else:
        message = message + total_error
    if warning != "":
        message = message + "-------------WARNINGS-------------\n" + warning

    print(message)
    return(message)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Validate GENIE files')

    parser.add_argument("fileType", type=str, choices = ['maf','clinical','fusion','cnv','vcf','seg','bed'],
                        help='File type that you are validating: maf, clinical, fusion, cnv, vcf, seg, bed')
    parser.add_argument("file", type=str, nargs="+",
                        help='File(s) that you are validating.  If you validation your clinical files and you have both sample and patient files, you must provide both')
    parser.add_argument("center", type=str, choices = ['MSK','GRCC','DFCI','NKI','JHU','MDA','VICC','UHN'],
                        help='Contributing Center')
    args = parser.parse_args()

    perform_main(args)
