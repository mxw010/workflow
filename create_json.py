
import pandas as pd
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


def main():
	global response, service
	creds = None
	if os.path.exists('token.pickle'):
		with open('token.pickle', 'rb') as token:
			creds = pickle.load(token)
	if not creds or not creds.valid:
		if creds and creds.expired and creds.refresh_token:
			creds.refresh(Request())
		else:
			flow = InstalledAppFlow.from_client_secrets_file(
				'credentials.json', SCOPES) # here enter the name of your downloaded JSON file
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
			if not values_input and not values_expansion:
				print('No data found.')


main()

df=pd.DataFrame(response[1:], columns=response[0])
#those variables are usful in making JSON file for ENCODE
select = [ 'Sample ID (multiplex)', 'Library ID', 'Library Type', 'Reference Genome', 'Target', 'Treatment', 'Timepoint', 'Cell Type', 'Replicate','Sequencing Modality' ]


#for a particular study
study='GSL-BS-1616'
workbook = df[select][df['Sample ID (multiplex)'] == study]
exp_ID = np.unique(workbook['Sample ID (multiplex)'])[0]
exp_ID_full = glob.glob("/Volumes/RESStanton/fastq/*" + study)[0]


#seperate ChIP-seq from ATAC-seq
#work on case ignore match
#For ChIP-seq
if re.search("ChIP", np.unique(workbook['Library Type'])[0]):
	#figure out if this is TF or histone marks
	#to get the target of ChIP and write the type of pipeline:
	target = set(np.unique(workbook['Target'])).difference({"Input"})
	if re.search('H2|H3|H4', repr(target)):
		chip.pipeline_type = "histone"
	else:
		chip.pipeline_type = "tf"


# create treatment conditions:
# Treatment

for sampleID in workbook['Library ID']: 
	#for paired end sequencing. is single end sequencing really still helpful? Maybe for chip-exo
	read1_fuzzy = exp_ID_full + "/**/" + sampleID + "*R1*.fastq.gz"
	read2_fuzzy = exp_ID_full + "/**/" + sampleID + "*R1*.fastq.gz"
	#glob.glob = find in BASH
	read1_path = glob.glob(read1_fuzzy, recursive=True)[0]
	read2_path = glob.glob(read2_fuzzy, recursive=True)[0]
	print(read1_path, read2_path)

