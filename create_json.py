import pandas as pd
import numpy as np
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow,Flow
from google.auth.transport.requests import Request
import os, glob
import pickle
import json
import re

SCOPES = ['https://www.googleapis.com/auth/spreadsheets']

#ID for Sequencing_Systems Epigenetics Group
spreadsheetId = '1gNHthk3lt6pYYSrFOPgjaQvy5AygMnjJakqyu77_FqY'

#this is the function to sync with google sheet
def main():
	global df, service
	creds = None
	if os.path.exists('token.pickle'):
		with open('token.pickle', 'rb') as token:
			creds = pickle.load(token)
	if not creds or not creds.valid:
		if creds and creds.expired and creds.refresh_token:
			creds.refresh(Request())
		else:
			flow = InstalledAppFlow.from_client_secrets_file(
				'/home/gdstantonlab/mxw010/Data/secrets/credentials.json', SCOPES) # here enter the name of your downloaded JSON file
			creds = flow.run_local_server(port=0)
		with open('token.pickle', 'wb') as token:
			pickle.dump(creds, token)
	service = build('sheets', 'v4', credentials=creds)
	# Call the Sheets API
	#Get the sheet ID
	sheet_metadata = service.spreadsheets().get(spreadsheetId=spreadsheetId).execute()
	sheets = sheet_metadata.get('sheets', '')
	for sheet in sheets:
		title = sheet.get("properties", {}).get("title", {})
		#pull information from the Experiments tab
		if re.search("Experiments", title):
			rangeName = title + '!A:Z'
			sheet_id = sheet.get("properties", {}).get("sheetId", {})
			#tab name is part of the range variable
			request = service.spreadsheets().values().get(spreadsheetId=spreadsheetId, range=rangeName)
			response = request.execute().get('values',[])
			if not response and not values_expansion:
				print('No data found.')
			df=pd.DataFrame(response[1:], columns=response[0])


main()

#backup: format is crap

#install packages:
#pip3 install --user google-auth-oauthlib
#pip3 install --user google-api-python-client
#pip3 install --user google-auth


#guide: https://www.twilio.com/blog/2017/02/an-easy-way-to-read-and-write-to-a-google-spreadsheet-in-python.html
#the scope needs to be changed to the following for this to work
# import gspread
# from oauth2client.service_account import ServiceAccountCredentials
 
# scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
# creds = ServiceAccountCredentials.from_json_keyfile_name('/home/gdstantonlab/mxw010/Data/secrets/stanton-bioinfor-1603567989841-a6075106e2cc.json', scope)
# client = gspread.authorize(creds)

# # Find a workbook by name and open the first sheet
# # Make sure you use the right name here.
# sheet = client.open("Sequencing_Systems Epigenetics Group").worksheet("Experiments")

# list_of_hashes = sheet.get_all_records()
# print(list_of_hashes)

#------start of the process

#for a particular study
study="GSL-LL-2096"

select = [ 'Sample ID (multiplex)', 'Library ID', 'Library Type', 'Reference Genome', 'Target', 'Treatment', 'Timepoint', 'Cell Type', 'Replicate','Sequencing Modality' ]

workbook = df[select][df['Sample ID (multiplex)'] == study]
exp_ID = np.unique(workbook['Sample ID (multiplex)'])[0]
#full path to fastq files
exp_ID_full = glob.glob("/home/gdstantonlab/lab/fastq/*" + study)[0]



#string manipulation:
#remove whitespace at the end of any cell
workbook.loc[:,workbook.columns] = workbook.apply(lambda x: x.str.strip())

#convert to all lower case
workbook.loc[:,'Treatment'] = workbook['Treatment'].str.upper()
workbook.loc[:,'Cell Type'] = workbook['Cell Type'].str.upper()
workbook.loc[:,'Reference Genome'] = workbook['Reference Genome'].str.upper()

#replace empty string in Treatment with "untreated"
workbook.loc[:,'Treatment'] = workbook['Treatment'].replace(r'^\s*$', 'untreated', regex=True)
#replace space in treatment with underscore
workbook.loc[:,'Treatment'] = workbook['Treatment'].replace(r' ', '_', regex=True)
workbook.loc[:,'Cell Type'] = workbook['Cell Type'].replace(r' ', '_', regex=True)
#replace empty string in Replicate with "rep1"
workbook.loc[:,'Replicate'] = workbook['Replicate'].replace(r'^\s*$', 'rep1', regex=True)
#replace empty string in Timepoint with 0hr
workbook.loc[:,'Timepoint'] = workbook['Timepoint'].replace(r'^\s*$', 'NA', regex=True)

#lower case replicate
workbook.loc[:,'Replicate'] = workbook['Replicate'].str.lower()

#remove all whitespace
workbook = workbook.replace(r' ',"",regex = True)
#convert n/a to NA
workbook = workbook.replace("n/a","NA")
#replace "/" with "_"
#workbook = workbook.replace(r"/", "_", regex = True)


if any(re.search("/", x) for x in workbook['Reference Genome']): 
	workbook['SpikeIn Genome'] = ''
	workbook['SpikeIn Cell Type'] = ''
	workbook['Primary Genome'] = ''
	workbook['Primary Cell Type'] = ''


#sequencing modality
#reference genomes+cell types
for i in range(0, len(workbook)):
#paired end of single end:
	if workbook.iloc[i]['Sequencing Modality'] is None:
		workbook.loc[workbook.index[i],'Sequencing Modality'] = 'pair_end'
	elif re.search("single",workbook.iloc[i]['Sequencing Modality'], re.IGNORECASE) or re.search("exo",workbook.iloc[i]['Library Type'], re.IGNORECASE):
		workbook.loc[workbook.index[i],'Sequencing Modality'] = 'single_end'
	#get primary genome, spike-in, and primary cell type
	else:
		workbook.iloc[i]['Sequencing Modality'] = 'paired_end'
	if re.search("/", workbook.iloc[i]['Reference Genome']):
		ref_gen = workbook.iloc[i]['Reference Genome']
		#python does not like it when you assign value over iloc
		workbook.loc[workbook.index[i],'Primary Genome'] = re.split('/', ref_gen)[0]
		workbook.loc[workbook.index[i],'SpikeIn Genome'] = re.split('/', ref_gen)[1]
		cell_type = workbook.iloc[i]['Cell Type']
		workbook.loc[workbook.index[i],'Primary Cell Type'] = re.split('/', cell_type)[0]
		workbook.loc[workbook.index[i],'SpikeIn Cell Type'] = re.split('/', cell_type)[1]


#common variables to be passed from BASH
#aligner, peak_caller, align_only, true_rep_only (depends on # of rep), always_use_pooled_ctl
#get reference genome:
# ref_gen = np.unique(workbook['Reference Genome'])
# if len(ref_gen) > 1:
# 	sys.exit("Aligning to more than 1 reference genome:" + ", ".join(ref_gen))
#work on how to align to multiple reference genome...

#makr analysis dir
#if not os.path.isfile("/Volumes/RESStanton/Analysis"):
#	os.mkdir("/Volumes/RESStanton/Analysis")

os.mkdir(study)
os.chdir(study)


#get workbook for chip-seq:
workbook_chip = pd.DataFrame()

for i in range(0,len(workbook)):
	if re.search("ChIP", workbook.iloc[i]['Library Type'], re.IGNORECASE):
		workbook_chip = workbook_chip.append(workbook.iloc[i])


#remove variables that doesn't vary
cond_var = ['Primary Cell Type', 'Treatment', 'Timepoint']
for comp in cond_var:
	if len(np.unique(workbook_chip[comp])) == 1:
		cond_var.remove(comp)
		

#conditions
conditions = workbook_chip[cond_var].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)


#figure out if this is TF or histone marks
#to get the target of ChIP and write the type of pipeline:
#convert all target to upper case
#can have multiple type of targets. result stored in a list
targets = set(j.upper() for j in np.unique(workbook_chip['Target'])).difference({"INPUT"})
targets = list(targets)
pipeline_types = [ "histone" if re.search('H2|H3|H4', repr(x)) else "tf" for x in targets]


for j in range(0,len(targets)):
	target = targets[j]
	os.mkdir(target)
	os.chdir(target)
	pipeline = pipeline_types[j]
	cond = conditions[workbook_chip['Target'].str.upper() == target]
	uniq_cond = np.unique(cond).tolist()
	#for combination of uniq setting + target, find all targets:
	for k in range(0,len(uniq_cond)):
		os.mkdir(uniq_cond[k])
		os.chdir(uniq_cond[k])
		data =pd.DataFrame()
		for i in range(0,len(workbook_chip)):
			y=workbook_chip.iloc[i][cond_var]
			if all(re.search(x, uniq_cond[k], re.IGNORECASE) for x in y) and re.search(target, workbook_chip.iloc[i]['Target'], re.IGNORECASE):
				data = data.append(workbook_chip.iloc[i][['Library ID', 'Target', 'Replicate']])
		#get input sample matching fastqs
		input = []
		for id in data['Library ID']:
			exact_cond = workbook_chip.loc[workbook_chip['Library ID'] == id,cond_var]
			input.append(workbook_chip.loc[workbook_chip[cond_var].eq(exact_cond).all(1), 'Library ID'].to_string(index=False).strip())
		#pair end
		if re.search("pair_end", workbook_chip.loc[workbook_chip['Library ID'] ==data.loc[data.index[0],'Library ID'],'Sequencing Modality'].to_string(index=False), re.IGNORECASE):
			pair_end = 'true'
		else:
			pair_end = 'false'
		primary_genome = workbook_chip.loc[workbook_chip['Library ID'] ==data.loc[data.index[0],'Library ID'],'Primary Genome'].to_string(index=False)
		#reference genome file for hg38 and mm10
		if re.search("HUMAN", primary_genome, re.IGNORECASE):
			genome_tsv = "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv"
		elif re.search("MOUSE", primary_genome, re.IGNORECASE):
			genome_tsv = "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/mm10.tsv"
		else:
			sys.exit("Not aligning to human/mouse genome; reference genome is: " + ref_gen)
		#make json file
		#os.getcwd() = getwd() 
		description = "ChIP-seq analysis for " + uniq_cond[k] + ". Target: " + target
		description = description + ". Samples: " + ", ".join(data['Library ID'].tolist()) + "."
		description = description + " Inputs: " + ", ".join(input) +"."
		seq_type = 'chip'
		item = { "chip.title": uniq_cond[k], "chip.description": description, "chip.pipeline_type": pipeline, "chip.genome_tsv": genome_tsv}
		#json_item = json.dumps(item, indent = 2, separators = (",",":"))
		other_par = {"chip.aligner" : "bowtie2", "chip.align_only" : 'false', "chip.paired_end" : 'true', "chip.ctl_paired_end" : 'true'}
		item.update(other_par)
		#what if no controls?
		#write fastqs for samples and controls:
		#number of replicates
		nrep=len(data)
		for i in range(1,nrep+1):
			rep = "rep" + str(i)
			sampleID = data.loc[data['Replicate'] == rep, 'Library ID'].item()
			inputID = input[i-1]
			if pair_end:
				read1_pattern = exp_ID_full + "/**/" + sampleID + "*R1*.fastq.gz"
				read1_pattern = read1_pattern.replace(' ',"")
				read2_pattern = exp_ID_full + "/**/" + sampleID + "*R2*.fastq.gz"
				read2_pattern = read2_pattern.replace(' ',"")
				fastqs_R1 = glob.glob(read1_pattern, recursive=True)[0].split()
				fastqs_R2 = glob.glob(read2_pattern, recursive=True)[0].split()
				read1 = "chip.fastqs_" + rep + "_R1"
				read2 =  "chip.fastqs_" + rep + "_R2"
				#controls
				control1_pattern = exp_ID_full + "/**/" + inputID + "*R1*.fastq.gz"
				control1_pattern = control1_pattern.replace(' ',"")
				control2_pattern = exp_ID_full + "/**/" + inputID + "*R2*.fastq.gz"
				control2_pattern = control2_pattern.replace(' ',"")
				control_R1 = glob.glob(control1_pattern, recursive=True)[0].split()
				control_R2 = glob.glob(control2_pattern, recursive=True)[0].split()
				control1 = "chip.ctl_fastqs_" + rep + "_R1"
				control2 =  "chip.ctl_fastqs_" + rep + "_R2"
				reads = {read1: fastqs_R1, read2: fastqs_R2, control1: control_R1, control2: control_R2}
			else:
				read1_pattern = exp_ID_full + "/**/" + sampleID + "*R1*.fastq.gz"
				read1_pattern = read1_pattern.replace(' ',"")
				fastqs_R1 = glob.glob(read1_pattern, recursive=True)[0].split()
				read1 = "chip.fastqs_" + rep + "_R1"
				#controls
				control1_pattern = exp_ID_full + "/**/" + inputID + "*R1*.fastq.gz"
				control1_pattern = control1_pattern.replace(' ',"")
				control_R1 = glob.glob(control1_pattern, recursive=True)[0].split()
				control1 = "chip.ctl_fastqs_" + rep + "_R1"
				reads = {read1: fastqs_R1, control1: control_R1}
			item.update(reads)
		#other parameters
		other_par = {"chip.peak_caller" : "macs2", "chip.true_rep_only" : 'true', "chip.always_use_pooled_ctl" : 'false'}
		item.update(other_par)
		#for some reason json.dump is not formatting the file correctly
		# json_item = json.dumps(item, indent=2, separators = (","," : "))
		# with open('atac.json', 'w') as outfile:
		# 	json.dump(json_item, outfile)
		file = open('working.json', 'w')
		file.write(json.dumps(item, indent=4, separators = (","," : "))) 
		file.close()  
		os.chdir("..")
	os.chdir("..")



#=====stopped here#

#if re.search(y['Cell Type'],uniq_cond[0], re.IGNORECASE) and y['Target'] == target and 
#if the combination of conditions and replicates does not match nrow of the data
#if len(np.unique(cond + "_" + workbook['Replicate'])) != len(workbook):
#	sys.exit("Condition configuration is wrong. Double check Treatment, Timepoint, Cell Type and Replicate!")





# #if the combination of conditions and replicates does not match nrow of the data
# if len(np.unique(cond + "_" + workbook['Replicate'])) != len(workbook):
# 	sys.exit("Condition configuration is wrong. Double check Treatment, Timepoint, Cell Type and Replicate!")


# for dir in uniq_cond: 
# 	os.mkdir(dir)
# 	os.chdir(dir)
# 	title = dir
# 	#os.getcwd() = getwd()
# 	description = "ATAC-seq analysis for " + title + ". Samples: " + ", ".join(workbook[cond.str.contains(dir)]['Sample ID (multiplex)'])
# 	pipeline_type = 'atac'
# 	item = { "title": title, "description": description, "pipeline_type": "atac", "genome_tsv": genome_tsv }
# 	json_item = json.dumps(item, indent = 2, separators = (",",":"))
# 	all_reps = workbook[cond.str.contains(dir)]
# 	for rep in all_reps['Replicate']:
# 		if not re.search("rep", rep):
# 			sys.exit("Replicate naming is wrong. Check the spreadsheet")


# 	for i in range(1,len(np.unique(all_reps['Replicate']))+1):
# 		rep = "rep" + str(i)
# 		sampleID = all_reps[all_reps['Replicate'].str.contains(rep)]['Library ID']
# 		read1_pattern = exp_ID_full + "/**/" + sampleID.to_string(index=False) + "*R1*.fastq.gz"
# 		read1_pattern = read1_pattern.replace(' ',"")
# 		read2_pattern = exp_ID_full + "/**/" + sampleID.to_string(index=False) + "*R2*.fastq.gz"
# 		read2_pattern = read2_pattern.replace(' ',"")
# 		fastqs_R1 = glob.glob(read1_pattern, recursive=True)[0]
# 		fastqs_R2 = glob.glob(read2_pattern, recursive=True)[0]
# 		read1 = "atac.fastqs_" + rep + "_R1"
# 		read2 =  "atac.fastqs_" + rep + "_R2"
# 		reads = {read1: fastqs_R1, read2: fastqs_R2}
# 		item.update(reads)

# 	#other parameters
# 	other_par = {"atac.paired_end" : paired_end, "atac.auto_detect_adapter" : 'true', "atac.align_cpu" : 8, 
# 				"atac.filter_cpu" : 4, "atac.bam2ta_cpu" : 4, "atac.jsd_cpu" : 4}
# 	item.update(other_par)
# 	#for some reason json.dump is not formatting the file correctly
# 	# json_item = json.dumps(item, indent=2, separators = (","," : "))
# 	# with open('atac.json', 'w') as outfile:
# 	# 	json.dump(json_item, outfile)
# 	file = open('atac.json', 'w')
# 	file.write(json.dumps(item, indent=4, separators = (","," : "))) 
# 	file.close()  
