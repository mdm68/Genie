import synapseclient
import calendar
import time
import argparse
import os
import multiprocessing.dummy as mp
import pandas as pd
import re
import shutil
#from processing_functions import *
import synapseutils as synu
import validateGENIE as validate_genie
import processing_functions


CENTERS_INPUT_SYNIDS = {"NKI":"syn5017764",
						"DFCI":"syn5017765",
						"GRCC":"syn5017766",
						"JHU":"syn5017774",
						"MSK":"syn5018076",
						"UHN":"syn5017775",
						"VICC":"syn5017776",
						"MDA":"syn7264740"}

CENTER_STAGING_SYNIDS = dict(JHU = "syn5016900",
							 DFCI = "syn5016917",
							 GRCC = "syn5016904",
							 NKI = "syn5016919",
							 MSK = "syn5548760",
							 UHN = "syn5016894",
							 VICC = "syn5016921",
							 MDA = "syn7073135")

VALIDATE_MAPPING = {'maf':validate_genie.validateMAF,
		            'clinical':validate_genie.validateClinical,
		            'vcf':validate_genie.validateVCF,
		            'cnv':validate_genie.validateCNV,
		            'fusion':validate_genie.validateFusion,
		            'seg':validate_genie.validateSEG,
		            'bed':validate_genie.validateBED}

VALIDATE_FILENAME = {'maf':"data_mutations_extended_%s.txt",
                     'clinical': ["data_clinical_supp_%s.txt", "data_clinical_supp_sample_%s.txt", "data_clinical_supp_patient_%s.txt"],
                     'vcf':"GENIE-%s-",
                     'cnv':"data_CNA_%s.txt",
                     'fusions':"data_fusions_%s.txt",
                     'seg':"genie_data_cna_hg19_%s.seg",
                     'bed':"%s-"}

def validateFiles(synId, center):
	"""
	This function walks through each center's input directory and validates every file
	"""
	syn = synapseclient.login()
	walked = synu.walk(syn, synId)
	clinicalpair = []
	invalid = []
	validFiles = []
	for dirpath, dirname, filenames in walked:
		for name, synid in filenames:
			file = []
			if name.endswith(".vcf") and VALIDATE_FILENAME['vcf'] % center in name:
				fileType = "vcf"
			elif name.endswith(".bed") and VALIDATE_FILENAME['bed'] % center in name:
				fileType = "bed"
			elif  VALIDATE_FILENAME['maf'] % center == name:
				fileType = "maf"
			elif VALIDATE_FILENAME['cnv'] % center == name:
				fileType = "cnv"
			elif VALIDATE_FILENAME['fusions'] % center == name:
				fileType = "fusions"
			elif VALIDATE_FILENAME['seg'] % center == name:
				fileType = "seg"
			elif VALIDATE_FILENAME['clinical'][0] % center == name:
				fileType = "clinical"
			elif name in VALIDATE_FILENAME['clinical'][1:]:
				clinicalpair.append(synid)
				if len(clinicalpair) == 2:
					fileType = "clinical"
				else:
					fileType = None
			else:
				fileType = None
			if fileType is not None:
				if len(clinicalpair) == 2:
					#Need to pass in both filepath and synid for processing files
					file = [(syn.get(i).path,i) for i in clinicalpair]
					clinicalpair = []
				else:
					file = [(syn.get(synid).path,synid)]
				#Validation only takes in a list of filepaths
				paths = [i[0] for i in file]
				message, valid = validate_genie.main(fileType, paths, center)
			else:
				print("%s: Cannot be processed" % name)
				valid = False
				invalid.append(name)
			if not valid:
				invalid.append(name)
			else:
				validFiles.extend(file)
	print(", ".join(invalid) + " can't be processed!")
	return(validFiles)

########################################################################
#Processing all other files
########################################################################
def processFiles(validFiles, center, path_to_GENIE):
	centerStagingFolder = os.path.join(path_to_GENIE, center, "staging")
	for filePath, synid in validFiles:
		inputFile = syn.get(synid,downloadFile=False)
		filename = os.path.basename(filePath)
		newPath = os.path.join(centerStagingFolder, filename)
		store = True
		if filename.startswith('data_clinical_supp'):
			filename = processing_functions.formatClinical(filePath, center, newPath)
		elif filename.startswith("data_CNA"):
			filename = processing_functions.formatCNA(filePath, center, newPath)
		elif filename.endswith(".seg"):
			formattedSeg = processing_functions.formatSEG(filePath, center, newPath)
		elif filename.startswith("data_fusions"):
			filename = processing_functions.formatFusion(filePath, center, newPath)
		elif filename.startswith("data_mutations") or filename.endswith(".vcf"):
			print("Please run python processing_vcfmaf.py")
		elif filename.endswith(".bed"):
			continue
		#elif filename.lower().startswith("data_linear_cna"):
		#	newfilename = "data_linear_CNA_%s.txt" % center
		else:
			store = False
		if store:
			shutil.copyfile(filePath, newPath)
			stagingFile = storeFile(newPath, stagingID = CENTER_STAGING_SYNIDS[center], used = inputFile, center = center, annotations= inputFile.annotations)

def build_parser():
	"""Set up argument parser and returns"""
	parser = argparse.ArgumentParser(
		description='GENIE center inputs to database')
	parser.add_argument('--genie_path', metavar='/path/folder', default='/Users/ThomasY/sage_projects/Genie_processing',
			help='Path to GENIE folder')
	return parser

#p = mp.Pool(6)
args = build_parser().parse_args()
syn = synapseclient.login(silent=True)

# ----------------------------------------
# Start input to staging process
# ----------------------------------------
path_to_GENIE = args.genie_path

# for center in CENTERS_INPUT_SYNIDS:
# 	print("Center: " + center)
# 	#validation = validateAnnotations(temp)
# 	validFiles = getFiles(CENTERS_INPUT_SYNIDS[center], center)
# 	if len(validFiles) > 0:
# 		processFiles(validFiles,center)
# 	else:
# 		print("%s does not have any valid files to process" % center)

center = "NKI"
print("Center: " + center)
#validation = validateAnnotations(temp)
validFiles = validateFiles(CENTERS_INPUT_SYNIDS[center], center)
print(validFiles)
# if len(validFiles) > 0:
# 	processFiles(validFiles,center, path_to_GENIE)
# else:
# 	print("%s does not have any valid files to process" % center)


