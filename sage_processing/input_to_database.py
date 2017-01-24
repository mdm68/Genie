import synapseclient
import calendar
import time
import argparse
import os
import multiprocessing.dummy as mp
import pandas as pd
import re
#from processing_functions import *
import synapseutils as synu
import validateGENIE as validate_genie


CENTERS_INPUT_ID = {"NKI":"syn5017764",
					"DFCI":"syn5017765",
					"GRCC":"syn5017766",
					"JHU":"syn5017774",
					"MSK":"syn5018076",
					"UHN":"syn5017775",
					"VICC":"syn5017776",
					"MDA":"syn7264740"}

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


def getFiles(synId, center):
	syn = synapseclient.login()
	walked = synu.walk(syn, synId)
	clinicalpair = []
	invalid = []
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
					file = [syn.get(i).path for i in clinicalpair]
					clinicalpair = []
				else:
					file = [syn.get(synid).path]
				message, valid = validate_genie.main(fileType, file, center)
			else:
				print("%s: Cannot be processed" % name)
				valid = False
				invalid.append(name)
			if not valid:
				invalid.append(name)
	print(", ".join(invalid) + " can't be processed!")





def findNewFiles(args, id):
	"""Performs query query to find changed entities in id. """
	
	centersInputId = {"syn5017764":"NKI",
					  "syn5017765":"DFCI",
					  "syn5017766":"GRCC",
					  "syn5017774":"JHU",
					  "syn5018076":"MSK",
					  "syn5017775":"UHN",
					  "syn5017776":"VICC",
					  "syn7264740":"MDA"}
	
	walked = synu.walk(syn, id) 
	for dirpath, dirname, filename in walked:
		print dirpath
		print dirname
		print filename


	QUERY = "select id, name, versionNumber, modifiedOn, modifiedByPrincipalId, nodeType from entity where parentId=='%s' and modifiedOn>%i" 
	t = calendar.timegm(time.gmtime())*1000
	project = syn.get(id)
	#Determine the last audit time or overide with lastTime
	if args.days is None:  #No time specified
		args.days = project.get('lastAuditTimeStamp', None)
		if args.days is None:  #No time specified and no lastAuditTimeStamp set
			args.days = t - ONEDAY*1.1
		else: #args.days came from annotation strip out from list
			args.days = args.days[0]  
	print t, args.days, id, (t-args.days)/float(ONEDAY), 'days'
	results = list(syn.chunkedQuery(QUERY % (id, args.days)))
	#Add the project and other metadata
	for r in results:
		r['projectId'] = id
		r['projectName'] = centersInputId[id]
		r['date'] = synapseclient.utils.from_unix_epoch_time(r['entity.modifiedOn']).strftime("%b/%d/%Y %H:%M")
		r['user'] = syn.getUserProfile(r['entity.modifiedByPrincipalId'])['userName']
		r['type'] = r['entity.nodeType']
		
	#Set lastAuditTimeStamp
	if args.updateProject:
		project.lastAuditTimeStamp = t
		try:
			project = syn.store(project)
		except synapseclient.exceptions.SynapseHTTPError:
			pass
	return results


########################################################################
#Processing all other files
########################################################################
def processFiles(fileIDs,center, path_to_GENIE):
	centerInputFolder = os.path.join(path_to_GENIE,center,"input")
	centerStagingFolder = os.path.join(path_to_GENIE,center,"staging")
	for i in fileIDs:
		inputFile = syn.get(i)
		oldfilename = inputFile.name
		newPath = os.path.join(centerInputFolder, oldfilename)
		store = True
		needMeta = True
		if "meta" not in oldfilename:
			if oldfilename.lower().startswith(('genie_clinical_data','data_clinical',"jhu_genie_sample","jhu_genie_patient")):
				oldfilename = formatClinical(inputFile.path, center,newPath)
				if "patient" in oldfilename.lower():
					newfilename = "data_clinical_supp_patient_%s.txt" % center
				elif "sample" in oldfilename.lower():
					newfilename = "data_clinical_supp_sample_%s.txt" % center
				else:
					newfilename = "data_clinical_supp_%s.txt" % center

			elif oldfilename.lower().startswith("data_cna"):
				oldfilename = formatCNA(inputFile.path, center, newPath)
				newfilename = "data_CNA_%s.txt" % center

			elif oldfilename.lower().startswith("data_linear_cna"):
				newfilename = "data_linear_CNA_%s.txt" % center
			elif oldfilename.lower().endswith(".seg"):
				formattedSeg = formatSEG(inputFile.path, center, newPath)
				newfilename = "genie_data_cna_hg19_%s.seg" % center
			elif oldfilename.lower().startswith("data_fusions"):
				oldfilename = formatFusion(inputFile.path, center, newPath)
				newfilename = "data_fusions_%s.txt" % center
			elif oldfilename.lower().startswith("data_mutation"):
				print("NEW MUTATION FILE! Run processing_vcfmaf.py on: %s" % i)
				store = False
			else:
				store = False
		else:
			filename = oldfilename.split(".")[0]
			newfilename = "%s_%s.txt" %(filename,center)
			needMeta = False
			store=False
		if store:
			newInputPath = os.path.join(centerStagingFolder,newfilename)
			os.system('cp %s %s' %(os.path.join(centerInputFolder,oldfilename),newInputPath))
			stagingFile = storeFile(newInputPath, center_staging_synIDs[center], inputFile,center,inputFile.annotations)
			#if needMeta:
			#	writeMeta(stagingFile,center)

def build_parser():
	"""Set up argument parser and returns"""
	parser = argparse.ArgumentParser(
		description='Checks for new/modified entities in a project.')
	parser.add_argument('--userId', dest='userId',
						help='User Id of individual to send report, defaults to current user.')
	parser.add_argument('--projects', '-p', metavar='projects', type=str, nargs='*',
			help='Synapse IDs of projects to be monitored.')
	parser.add_argument('--days', '-d', metavar='days', type=float, default=None,
			help='Find modifications in the last days')
	parser.add_argument('--updateProject', dest='updateProject',  action='store_true',
			help='If set will modify the annotations by setting lastAuditTimeStamp to the current time on each project.')
	parser.add_argument('--config', metavar='file', dest='configPath',  type=str,
			help='Synapse config file with user credentials (overides default ~/.synapseConfig)')
	parser.add_argument('--genie_path', metavar='/path/folder', default='/Users/ThomasY/sage_projects/Genie_processing',
			help='Path to GENIE folder')
	return parser

p = mp.Pool(6)
args = build_parser().parse_args()
args.days = None if args.days is None else calendar.timegm(time.gmtime())*1000 - args.days*ONEDAY
if args.configPath is not None:
	syn=synapseclient.Synapse(skip_checks=True, configPath=args.configPath)
else:
	syn=synapseclient.Synapse(skip_checks=True)
syn.login(silent=True) 
args.userId = syn.getUserProfile()['ownerId'] if args.userId is None else args.userId

#query each project then combine into long list
entityList = p.map(lambda project: findNewFiles(args, project), args.projects)
entityList = [item for sublist in entityList for item in sublist]
#Filter out projects and folders
entityList = [e for e in entityList if e['entity.nodeType'] not in ['project', 'folder']]
print 'Total number of entities = ', len(entityList)

newEntities = pd.DataFrame(entityList)
print(newEntities)

if len(entityList) == 0:
	uniqueCenters = []
else:
	uniqueCenters = pd.unique(newEntities['projectName'])

# ----------------------------------------
# Start input to staging process
# ----------------------------------------
path_to_GENIE = args.genie_path
center_staging_synIDs = dict(JHU = "syn5016900",
							 DFCI = "syn5016917",
							 GRCC = "syn5016904",
							 NKI = "syn5016919",
							 MSK = "syn5548760",
							 UHN = "syn5016894",
							 VICC = "syn5016921")

for center in uniqueCenters:
	temp = newEntities[newEntities['projectName']==center]
	temp = temp.reset_index()
	print("Center: " + center)
	validation = validateAnnotations(temp)
	#if validation:
	# if center == "GRCC":
	# 	cbssegfiles = [temp['entity.id'][i] for i,name in enumerate(temp['entity.name']) if "cbs" in name]
	# 	sologrdfiles = [temp['entity.id'][i] for i,name in enumerate(temp['entity.name']) if "solo.grd" in name]
	# 	IDmapping = processCBSSEG(cbssegfiles, center, path_to_GENIE)
	# 	processLinearCNA(sologrdfiles,center, path_to_GENIE, IDmapping)
	processFiles(temp['entity.id'],center,path_to_GENIE)
	#else:
	#if not validation:
	#	print("%s files don't have correct annotations:" % center)
#		print(validation)


