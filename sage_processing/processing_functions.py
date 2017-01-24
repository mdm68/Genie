import synapseclient
from synapseclient import File, Table
import calendar
import time
import argparse
import os
import multiprocessing.dummy as mp
import pandas as pd
import numpy as np
import re
import subprocess

syn = synapseclient.login()

center_staging_synIDs = dict(JHU = "syn5016900",
							 DFCI = "syn5016917",
							 GRCC = "syn5016904",
							 NKI = "syn5016919",
							 MSK = "syn5548760",
							 UHN = "syn5016894",
							 VICC = "syn5016921",
							 MDA = "syn7073135")

clinical_column_mapping = {"GENIE_PATIENT_ID": "PATIENT_ID",
						   "NAACCR_SEX_CODE" : "SEX",
						   "NAACCR_RACE_CODE_PRIMARY" : "PRIMARY_RACE",
						   "NAACCR_RACE_CODE_SECONDARY" : "SECONDARY_RACE",
						   "NAACCR_RACE_CODE_TERTIARY" : "TERTIARY_RACE",
						   "NAACCR_ETHNICITY_CODE" : "ETHNICITY",
						   "GENIE_SAMPLE_ID" : "SAMPLE_ID",
						   "AGE_SEQ_REPORT_DAYS" : "AGE_AT_SEQ_REPORT",
						   "PATIENT_BIRTH_YEAR":"BIRTH_YEAR"}

########################################################################
# Check if GENIE ID is labelled correctly
########################################################################
def checkGenieId(ID,center):
	if str(ID).startswith("%s-" % center):
		return('GENIE-%s' % str(ID))
	elif not str(ID).startswith('GENIE-%s-' % center):
		return('GENIE-%s-%s' % (center, str(ID)))
	else:
		return(str(ID))

########################################################################
#Storing Files along with annotations
########################################################################
def storeFile(fileName, stagingID, used, center, annotations, meta=False):
	print("STORING FILES")
	fileEnt = File(fileName, parent = stagingID)
	fileEnt.center = center
	fileEnt.dataSubType = annotations.get("dataSubType",'')
	fileEnt.dataType = annotations.get("dataType",'')
	fileEnt.disease = 'cancer'
	fileEnt.fileType = annotations.get("fileType",'')
	fileEnt.organism = 'Homo Sapiens'
	fileEnt.platform = annotations.get("platform",'')
	fileEnt.tissueSource = annotations.get("tissueSource",'')
	fileEnt.consortium = 'GENIE'
	if meta:
		fileEnt.fileType = "txt"
		fileEnt.dataType = "meta"
	fileEnt.fileStage = "staging"
	ent = syn.store(fileEnt,annotations = used)
	return(ent)

####################################################################################
# UPDATING DATABASE
####################################################################################
def updateDatabase(database, new_dataset, databaseSynId, checkBy):
	"""
	Updates synapse tables by a row identifier with another dataset that has the same number and order of columns
	
	:param database:   	   The synapse table (pandas dataframe)
	:param new_dataset:    New dataset (pandas dataframe)
	:param databaseSynId   Synapse Id of the database table
	:param checkBy:        Column to compare both datasets by

	:returns:      		   Don't know yet	
	"""
	updatedSet = database.apply(lambda x: _updateRows(x, new_dataset, checkBy),axis=1)
	updatedSet = updatedSet[~updatedSet[checkBy].isnull()]
	#All new rows
	newSet =  new_dataset[~new_dataset[checkBy].isin(database[checkBy])]
	#All deleted rows (This assumes that all data that don't show up in the new uploaded data should be deleted...)
	deleteSets = database[~database[checkBy].isin(new_dataset[checkBy])]
	if not deleteSets.empty:
		deleteRows = syn.delete(Table(syn.get(databaseSynId), deleteSets))
	else:
		print("No deleted rows")
	#updatedSet = updatedSet.append(newSet)
	if not updatedSet.empty:
		table = syn.store(Table(syn.get(databaseSynId), updatedSet))	
	else:
		print("No updated rows")
	if not newSet.empty:
		table = syn.store(Table(syn.get(databaseSynId), newSet))	
	else:
		print("No new rows")

def _updateRows(database_record, new_dataset, checkBy):
	"""
	Checks if the rows are the same, if not then update
	
	:param database_record:   	   One row in the database (pandas series)
	:param new_dataset:    		   New dataset (pandas dataframe)
	:param checkBy:        		   Column to compare both datasets by

	:returns:      		   		   List of updated values
	"""
	#Make sure there is only one unique match
	new_record = new_dataset[new_dataset[checkBy] == database_record[checkBy]]
	#Must make -NAN because np.nan != np.nan
	new_record = new_record.fillna("-NAN")
	database_record = database_record.fillna("-NAN")
	if all([old == new for old, new in zip(database_record.values, new_record.values[0])]):
		#Make all values in the new record null so its removed in the upload
		for i in new_record:
			new_record[i] = np.nan
	else:
		#Have to use a for loop because there are integers in the series, can't do comparison
		for i in new_record:
			if str(new_record[i].values[0]) == "-NAN":
				new_record[i] = np.nan
	return(list(new_record.values[0]))

########################################################################
# Format SEG files
########################################################################
def formatSEG(filePath,center,newPath):
	seg = pd.read_csv(filePath, sep="\t")
	seg.columns = [col.upper() for col in seg.columns]
	newsamples = [checkGenieId(i,center) for i in seg['ID']]
	seg['ID'] = newsamples
	seg = seg.drop_duplicates()
	seg.rename({'LOC.START':'LOCSTART','LOC.END':'LOCEND','SEG.MEAN':'SEGMEAN'})
	seg['UNIQUE_KEY'] = seg['ID'] + seg['CHROM'] + seg['LOCSTART'].astype(str)  + seg['LOCEND'].astype(str)

	databaseSynId = "syn7893341"
	checkBy = "UNIQUE_KEY"
	seg_database = syn.tableQuery('SELECT * FROM %s where %s in (%s)' % (databaseSynId, checkBy, ",".join("'" + seg[checkBy]+"'")))
	seg_database = seg_database.asDataFrame()[cols]

	#newClinical[seqColumn + "_NUMERICAL"] = [int(year) if checkInt(year) else np.nan for year in newClinical[seqColumn]]
	Fusion = Fusion[cols]
	updateDatabase(seg_database, seg, databaseSynId, checkBy)
#def updateDatabase(database, new_dataset, databaseSynId, checkBy):

	seg.to_csv(newPath,sep="\t",index=False)
	return(os.path.basename(newPath))


########################################################################
# Format Clinical files
########################################################################
def formatClinical(filePath, center, newPath):
	print(filePath)
	clinical = pd.read_csv(filePath, sep="\t", comment="#")
	patient= False
	sample = False
	patientCols = ["PATIENT_ID","SEX","PRIMARY_RACE","SECONDARY_RACE",
				   "TERTIARY_RACE","ETHNICITY","BIRTH_YEAR","CENTER"]
	sampleCols = ["SAMPLE_ID","AGE_AT_SEQ_REPORT","ONCOTREE_CODE","SAMPLE_TYPE",
				  "SEQ_ASSAY_ID"]
	if "patient" in filePath.lower():
		newClinical = pd.DataFrame(columns=patientCols)
		patient = True
	elif "sample" in filePath.lower():
		newClinical = pd.DataFrame(columns=sampleCols)
		sample = True
	else:
		newClinical = pd.DataFrame(columns=patientCols + sampleCols)
		sample = True 
		patient = True

	for key in clinical_column_mapping:
		if clinical.get(key) is not None:
			clinical[clinical_column_mapping[key]] = clinical[key]

	new = clinical.merge(newClinical,how='outer')
	new = new.drop(new.columns[~new.columns.isin(newClinical.columns)],1)

	SEX_MAPPING_ENT = syn.tableQuery('SELECT * FROM syn7434222')
	SEX_MAPPING = SEX_MAPPING_ENT.asDataFrame()

	RACE_MAPPING_ENT = syn.tableQuery('SELECT * FROM syn7434236')
	RACE_MAPPING = RACE_MAPPING_ENT.asDataFrame()
	#Fill NA's for the ones that are uncoded
	RACE_MAPPING = RACE_MAPPING.fillna("")

	ETHNICITY_MAPPING_ENT = syn.tableQuery('SELECT * FROM syn7434242')
	ETHNICITY_MAPPING = ETHNICITY_MAPPING_ENT.asDataFrame()
	#Fill NA's for the ones that are uncoded
	ETHNICITY_MAPPING = ETHNICITY_MAPPING.fillna("")

	SAMPLE_TYPE_ENT = syn.tableQuery('SELECT * FROM syn7434273')
	SAMPLE_TYPE = SAMPLE_TYPE_ENT.asDataFrame()
	#Attach MSK to centers
	new = new.fillna("")
	newClinical = new.apply(lambda x: update_clinical(x, center, SEX_MAPPING, RACE_MAPPING, ETHNICITY_MAPPING, SAMPLE_TYPE,patient),1)	
	if newClinical.get("CENTER") is not None:
		newClinical['CENTER'] = center

	if patient:
		seqColumn = "BIRTH_YEAR"
		cols = patientCols
		cols.append(seqColumn + "_NUMERICAL")
		databaseSynId = "syn7517669"
		checkBy = "PATIENT_ID"
		#Can't do this checkBy query because we need to remove samples that don't exist in the dataset from the database
		clinical_database = syn.tableQuery('SELECT * FROM %s where %s in (%s)' % (databaseSynId, checkBy, ",".join("'" + newClinical[checkBy]+"'")))
		clinical_database = clinical_database.asDataFrame()[cols]
		newClinical[seqColumn + "_NUMERICAL"] = [int(year) if checkInt(year) else np.nan for year in newClinical[seqColumn]]
		newClinical = newClinical[patientCols].drop_duplicates()
		updateDatabase(clinical_database, newClinical, databaseSynId, "PATIENT_ID")
	if sample:
		seqColumn = "AGE_AT_SEQ_REPORT"
		cols = sampleCols
		cols.append(seqColumn + "_NUMERICAL")
		databaseSynId = "syn7517674"
		checkBy = "SAMPLE_ID"
		clinical_database = syn.tableQuery('SELECT * FROM %s where %s in (%s)' % (databaseSynId, checkBy, ",".join("'" + newClinical[checkBy]+"'")))
		clinical_database = clinical_database.asDataFrame()[cols]
		newClinical[seqColumn + "_NUMERICAL"] = [int(year) if checkInt(year) else np.nan for year in newClinical[seqColumn]]
		newClinical = newClinical[sampleCols]
		updateDatabase(clinical_database, newClinical, databaseSynId, "SAMPLE_ID")

	newClinical.to_csv(newPath, sep="\t", index=False)
	return(os.path.basename(newPath))

#Check if an item can become an integer
def checkInt(element):
	try:
	    int(element)
	    return(True)
	except ValueError:
	    return(False)

#Get mapping code
def getCODE(mapping, key):
	value = mapping['CBIO_LABEL'][mapping['CODE'] == key].values
	if len(value) >0:
		return(value[0])
	else:
		return("")

#Update clinical file with the correct mappings
def update_clinical(x, center,SEX_MAPPING,RACE_MAPPING,ETHNICITY_MAPPING,SAMPLE_TYPE, patient):	
	#Can only do this for patient data, because patient is tagged as True when its a flat clinical file or only patient file
	#No check for sample file
	if patient:
		#TRIM EVERY COLUMN MAKE ALL DASHES 
		#PATIENT ID
		x['PATIENT_ID'] = checkGenieId(x['PATIENT_ID'], center)
		# RACE
		if x.get('PRIMARY_RACE') is not None:
			x['PRIMARY_RACE'] = getCODE(RACE_MAPPING, x['PRIMARY_RACE'])
		if x.get('SECONDARY_RACE') is not None:
			x['SECONDARY_RACE'] = getCODE(RACE_MAPPING, x['SECONDARY_RACE'])
		if x.get('TERTIARY_RACE') is not None:
			x['TERTIARY_RACE'] = getCODE(RACE_MAPPING, x['TERTIARY_RACE'])
		# ETHNICITY
		if x.get('ETHNICITY') is not None:
			x['ETHNICITY'] = getCODE(ETHNICITY_MAPPING, x['ETHNICITY'])
		# BIRTH YEAR (Check if integer)
		if checkInt(x['BIRTH_YEAR']):
			x['BIRTH_YEAR'] = int(x['BIRTH_YEAR'])
		# SEX
		x['SEX'] = getCODE(SEX_MAPPING, x['SEX'])
	#TRIM EVERY COLUMN MAKE ALL DASHES 
	#SAMPLE ID
	if x.get('SAMPLE_ID') is not None:
		x['SAMPLE_ID'] = checkGenieId(x['SAMPLE_ID'], center)
	#AGE AT SEQ REPORT
	if x.get('AGE_AT_SEQ_REPORT') is not None:
		if checkInt(x['AGE_AT_SEQ_REPORT']):
			x['AGE_AT_SEQ_REPORT'] = int(x['AGE_AT_SEQ_REPORT'])
	#SEQ ASSAY ID
	if x.get('SEQ_ASSAY_ID') is not None:
		if not str(x['SEQ_ASSAY_ID']).startswith(center) and str(x['SEQ_ASSAY_ID']) != "":
			x['SEQ_ASSAY_ID'] = "%s-%s" % (center, str(x['SEQ_ASSAY_ID']))
		x['SEQ_ASSAY_ID'] = x['SEQ_ASSAY_ID'].replace('_','-')
	#SAMPLE_TYPE
	if x.get('SAMPLE_TYPE') is not None:
		x['SAMPLE_TYPE'] = getCODE(SAMPLE_TYPE, x['SAMPLE_TYPE'])
	#Trim spaces
	for i in x.keys():
		if isinstance(x[i],str):
			x[i] = x[i].strip(" ")
	return(x)

########################################################################
# Format Fusion files
########################################################################
def formatFusion(filePath, center, newPath):
	Fusion = pd.read_csv(filePath, sep="\t",comment="#")
	Fusion.columns = [col.upper() for col in Fusion.columns]
	#newCols = [' '.join(i.split("_")).title().replace(' ','_') for i in Fusion.columns.values]
	#Fusion.columns = newCols
	Fusion['CENTER'] = center
	newsamples = [checkGenieId(i,center) for i in Fusion['TUMOR_SAMPLE_BARCODE']]
	Fusion['TUMOR_SAMPLE_BARCODE'] = newsamples

	cols = ['HUGO_SYMBOL','ENTREZ_GENE_ID','CENTER','TUMOR_SAMPLE_BARCODE','FUSION','DNA_SUPPORT','RNA_SUPPORT','METHOD','FRAME','COMMENTS']
	if Fusion.get("COMMENTS") is None:
		Fusion['COMMENTS'] = ""

	Fusion['COMMENTS'] = Fusion['COMMENTS'].fillna("")
	Fusion['ENTREZ_GENE_ID'] = Fusion['ENTREZ_GENE_ID'].fillna(0)

	Fusion = Fusion.drop_duplicates()

	Fusion['UNIQUE_KEY'] = Fusion['TUMOR_SAMPLE_BARCODE'] + Fusion['CENTER'] + Fusion['HUGO_SYMBOL'] + Fusion['FUSION'] + Fusion['FRAME'] + Fusion['COMMENTS']
	Fusion['ENTREZ_GENE_ID'] = [int(float(i)) for i in Fusion['ENTREZ_GENE_ID']]

	#cols.append(seqColumn + "_NUMERICAL")
	databaseSynId = "syn7893268"
	checkBy = "UNIQUE_KEY"
	fusion_database = syn.tableQuery('SELECT * FROM %s where %s in (%s)' % (databaseSynId, checkBy, ",".join("'" + Fusion[checkBy]+"'")))
	fusion_database = clinical_database.asDataFrame()[cols]
	fusion_database['UNIQUE_KEY'] = 

	#newClinical[seqColumn + "_NUMERICAL"] = [int(year) if checkInt(year) else np.nan for year in newClinical[seqColumn]]
	Fusion = Fusion[cols]

	updateDatabase(clinical_database, newClinical, databaseSynId, "PATIENT_ID")

	Fusion.to_csv(newPath, sep="\t",index=False)
	return(os.path.basename(newPath))

########################################################################
# Format CNA files
########################################################################
def formatCNA(filePath, center, newPath):
	CNA = pd.read_csv(filePath, sep="\t",comment="#")
	CNA.columns = [col.upper() for col in CNA.columns]
	if CNA.get("ENTREZ_GENE_ID") is not None:
		del CNA['ENTREZ_GENE_ID']
	symbols = CNA['HUGO_SYMBOL']
	del CNA['HUGO_SYMBOL']
	#Update all to int
	CNA = CNA.fillna(0)
	#CNA = CNA.applymap(int)
	newsamples = [checkGenieId(i,center) for i in CNA.columns]
	CNA.columns = newsamples
	#Transpose matrix
	CNA = CNA.transpose()
	CNA.columns = symbols.tolist()
	CNA = CNA.reset_index()
	CNA = CNA.rename_axis({'index':'TUMOR_SAMPLE_BARCODE'},axis="columns")
	CNA['CENTER'] = center

	#remove the 0.0, 1.0 and 2.0
	# os.system("sed 's/[.]0//g' %s > %s" % (newPath + "temp", newPath))
	# os.remove(newPath + "temp")
	#Issue of MAX 152 columns in a table

	cols = patientCols		
	databaseSynId = "syn7517669"
	checkBy = "TUMOR_SAMPLE_BARCODE"
	clinical_database = syn.tableQuery('SELECT * FROM %s where CENTER in %s' % (databaseSynId, center))
	clinical_database = clinical_database.asDataFrame()[cols]
	newClinical[seqColumn + "_NUMERICAL"] = [int(year) if checkInt(year) else np.nan for year in newClinical[seqColumn]]
	newClinical = newClinical[patientCols].drop_duplicates()
	updateDatabase(clinical_database, newClinical, databaseSynId, "PATIENT_ID")


	CNA.to_csv(newPath, sep="\t",index=False)
	return(os.path.basename(newPath))

########################################################################
# GRCC FILES
########################################################################
def processCBSSEG(fileIDs,center, path_to_GENIE):
	temp = syn.chunkedQuery('select id, name from file where parentId == "syn5017766"')
	segFiles = [i['file.id'] for i in temp if i['file.name'].endswith(".cbs")]
	segDF = pd.DataFrame(columns = ['ID','CHROM','LOC.START','LOC.END','NUM.MARK','SEG.MEAN'])
	for entityId in segFiles:
		newsegEnt = syn.get(entityId)
		newseg = pd.read_csv(newsegEnt.path,sep="\t")
		newseg.columns = ['ID','CHROM','LOC.START','LOC.END','num.mark','seg.mean']
		newseg['ID'] = newsegEnt.name.replace(".cbs","")
		newseg['chromosome'] = newseg['chrom']
		del newseg['chrom']
		segDF = segDF.append(newseg)
	segDF.to_csv("genie_data_cna_hg19_GRCC.seg",sep="\t",index=False)
	return("genie_data_cna_hg19_GRCC.seg")

#########################################################################
##Process Linear CNA file
#########################################################################
def processLinear():
	temp = syn.chunkedQuery('select id, name from file where parentId == "syn5017766"')
	cbsFiles = [i['file.id'] for i in temp if "solo.grd" in i['file.name']]
	linearCNA = pd.DataFrame(columns = ["Hugo_symbol"])
	for entityId in cbsFiles:
		sologrdent = syn.get(entityId)
		newname = sologrdent.name.split(".")[0].replace("_solo","")
		genelogrratio = pd.DataFrame(columns=["Hugo_symbol",newname])
		sologrd = pd.read_csv(sologrdent.path,sep="\t")
		for i in sologrd.index:
			row = sologrd.iloc[i]
			if row['Genes.all']!='0':
				genes = re.findall(".+\((.+)\)",row['Genes.all'])[0]
				genes = set(genes.split(","))
				temp = pd.DataFrame(index= range(len(genes)),columns=["Hugo_symbol",newname])
				temp['Hugo_symbol'] = genes
				temp[newname] = row['log2(ratio)']
				genelogrratio = genelogrratio.append(temp)
		linearCNA = pd.merge(linearCNA, genelogrratio,on='Hugo_symbol',how='outer')
		linearCNA = linearCNA.drop_duplicates()
	linearCNA.to_csv("data_linear_CNA_GRCC.txt", sep="\t",index=False)