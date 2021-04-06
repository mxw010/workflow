#install packages:
#pip3 install --user google-auth-oauthlib
#pip3 install --user google-api-python-client
#pip3 install --user google-auth


#guide: https://www.twilio.com/blog/2017/02/an-easy-way-to-read-and-write-to-a-google-spreadsheet-in-python.html
#the scope needs to be changed to the following for this to work


scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name('/home/gdlessnicklab/mxw010/Data/secrets/stanton-bioinfor-1603567989841-a6075106e2cc.json', scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet
# Make sure you use the right name here.
sheet = client.open("Sequencing_Systems Epigenetics Group").worksheet("Experiments")

list_of_hashes = sheet.get_all_records()
print(list_of_hashes)

#=============end here

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
				'/Users/mxw010/Documents/python_test/credentials.json', SCOPES) # here enter the name of your downloaded JSON file
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


#for a particular study
study="GSL-CJ-1833"

select = [ 'Sample ID (multiplex)', 'Library ID', 'Library Type', 'Reference Genome', 'Target', 'Treatment', 'Timepoint', 'Cell Type', 'Replicate','Sequencing Modality' ]

workbook = df[select][df['Sample ID (multiplex)'] == study]
exp_ID = np.unique(workbook['Sample ID (multiplex)'])[0]
exp_ID_full = glob.glob("/Volumes/RESStanton/fastq/*" + study)[0]

#string manipulation:
#convert to all lower case
workbook['Treatment'] = workbook['Treatment'].str.lower()
workbook['Cell Type'] = workbook['Cell Type'].str.upper()
workbook['Reference Genome'] = workbook['Reference Genome'].str.upper()
#remove all whitespace
workbook = workbook.replace(r' ',"",regex = True)
#convert n/a to NA
workbook = workbook.replace("n/a","NA")
#replace "/" with "_"
workbook = workbook.replace(r"/", "_", regex = True)
#replace empty string in Treatment with "untreated"
workbook['Treatment'] = workbook['Treatment'].replace(r'^\s*$', 'untreated', regex=True)
#replace space in treatment with underscore
workbook['Treatment'] = workbook['Treatment'].replace(r' ', '_', regex=True)
#replace empty string in Replicate with "rep1"
workbook['Replicate'] = workbook['Replicate'].replace(r'^\s*$', 'rep1', regex=True)
#replace empty string in Timepoint with 0hr
workbook['Timepoint'] = workbook['Timepoint'].replace(r'^\s*$', 'NA', regex=True)

#makr analysis dir
if not os.path.isfile("/Volumes/RESStanton/Analysis"):
	os.mkdir("/Volumes/RESStanton/Analysis")


os.mkdir(study)
os.chdir("/Volumes/RESStanton/Analysis/"+ study)

#paired end of single end:
if any("single" in x for x in workbook['Sequencing Modality']) | any("exo" in x for x in workbook['Library Type']):
	paired_end = 'false'
else:
	paired_end = 'true'


#common variables to be passed from BASH
#aligner, peak_caller, align_only, true_rep_only (depends on # of rep), always_use_pooled_ctl
#get reference genome:
ref_gen = np.unique(workbook['Reference Genome'])
if len(ref_gen) > 1:
	sys.exit("Aligning to more than 1 reference genome:" + ", ".join(ref_gen))
#work on how to align to multiple reference genome...

ref_gen = ref_gen[0]
#if there is spike in:
if re.search("/", ref_gen):
	ref_gen = re.split('/', ref_gen)[0]

#reference genome file for hg38 and mm10
if re.search("HUMAN", ref_gen):
	genome_tsv = "https://storage.googleapis.com/encode-pipeline-genome-data/hg38_caper.tsv"
elif re.search("MOUSE", ref_gen):
	genome_tsv = "https://storage.googleapis.com/encode-pipeline-genome-data/mm10_caper.tsv"
else:
	sys.exit("Not aligning to human/mouse genome; reference genome is: " + ref_gen)


#making json file for ATAC-seq:
if re.search("ATAC", np.unique(workbook['Library Type'])[0], re.IGNORECASE):

	#figure out what to name each directory
	#celltype-treatment-timepoint-replicate
	cond = ['Cell Type', 'Treatment', 'Timpoint', 'Replicate']
	for comp in cond:
		if len(np.unique(workbook[comp])) == 1:
			cond = cond - comp


####inside the ATAC loop
#remove variables that doesn't vary
cond_var = ['Cell Type', 'Treatment', 'Timepoint']
for comp in cond_var:
	if len(np.unique(workbook[comp])) == 1:
		cond_var.remove(comp)

#conditions
cond = workbook[cond_var].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

uniq_cond = np.unique(cond).tolist()

#if the combination of conditions and replicates does not match nrow of the data
if len(np.unique(cond + "_" + workbook['Replicate'])) != len(workbook):
	sys.exit("Condition configuration is wrong. Double check Treatment, Timepoint, Cell Type and Replicate!")


for dir in uniq_cond: 
	os.mkdir(dir)
	os.chdir(dir)
	title = dir
	#os.getcwd() = getwd()
	description = "ATAC-seq analysis for " + title + ". Samples: " + ", ".join(workbook[cond.str.contains(dir)]['Sample ID (multiplex)'])
	pipeline_type = 'atac'
	item = { "title": title, "description": description, "pipeline_type": "atac", "genome_tsv": genome_tsv }
	json_item = json.dumps(item, indent = 2, separators = (",",":"))
	all_reps = workbook[cond.str.contains(dir)]
	for rep in all_reps['Replicate']:
		if not re.search("rep", rep):
			sys.exit("Replicate naming is wrong. Check the spreadsheet")


	for i in range(1,len(np.unique(all_reps['Replicate']))+1):
		rep = "rep" + str(i)
		sampleID = all_reps[all_reps['Replicate'].str.contains(rep)]['Library ID']
		read1_pattern = exp_ID_full + "/**/" + sampleID.to_string(index=False) + "*R1*.fastq.gz"
		read1_pattern = read1_pattern.replace(' ',"")
		read2_pattern = exp_ID_full + "/**/" + sampleID.to_string(index=False) + "*R2*.fastq.gz"
		read2_pattern = read2_pattern.replace(' ',"")
		fastqs_R1 = glob.glob(read1_pattern, recursive=True)[0]
		fastqs_R2 = glob.glob(read2_pattern, recursive=True)[0]
		read1 = "atac.fastqs_" + rep + "_R1"
		read2 =  "atac.fastqs_" + rep + "_R2"
		reads = {read1: fastqs_R1, read2: fastqs_R2}
		item.update(reads)

	#other parameters
	other_par = {"atac.paired_end" : paired_end, "atac.auto_detect_adapter" : 'true', "atac.align_cpu" : 8, 
				"atac.filter_cpu" : 4, "atac.bam2ta_cpu" : 4, "atac.jsd_cpu" : 4}
	item.update(other_par)
	#for some reason json.dump is not formatting the file correctly
	# json_item = json.dumps(item, indent=2, separators = (","," : "))
	# with open('atac.json', 'w') as outfile:
	# 	json.dump(json_item, outfile)
	file = open('atac.json', 'w')
	file.write(json.dumps(item, indent=4, separators = (","," : "))) 
	file.close()  


# #For ChIP-seq
# if re.search("ChIP", np.unique(workbook['Library Type'])[0], re.IGNORECASE):
# 	#figure out if this is TF or histone marks
# 	#to get the target of ChIP and write the type of pipeline:
# 	#convert all target to upper case
# 	target = set(j.upper() for j in np.unique(workbook['Target'])).difference({"INPUT"})
# 	if re.search('H2|H3|H4', repr(target)):
# 		chip.pipeline_type = "histone"
# 	else:
# 		chip.pipeline_type = "tf"
