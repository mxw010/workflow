import numpy as np
import pandas as pd
import os, glob
import pickle
import json
import re


design_sheet = pd.read_csv('Design_sheet.csv', header=0, index_col=0)
primary = pd.read_csv('PrimaryQC.csv', header=0,index_col=False)
spikeIn = pd.read_csv('SpikeInQC.csv', header=0,index_col=False)

primary.rename(columns={'Primary Genome':'Genome'}, inplace=True)
primary.rename(columns={'Primary Cell Type':'Cell'}, inplace=True)
spikeIn.rename(columns={'SpikeIn Genome':'Genome'}, inplace=True)
spikeIn.rename(columns={'SpikeIn Cell Type':'Cell'}, inplace=True)

basedir = '/home/gdstantonlab/mxw010/Data/Ben/MYOD/GSL-BS-2095'

uni_input = primary[primary['Target'] == "INPUT"]['Library ID']
matching_sample = [ design_sheet[design_sheet['Input'] == x].iloc[0]['Sample'] for x in uni_input ]
input_data = pd.DataFrame({'Sample':matching_sample, 'Input':uni_input})

primary['total_reads'] = spikeIn['total_reads'] =''
primary['mapped'] = spikeIn['mapped'] =''
primary['proper_mapped'] = spikeIn['proper_mapped'] = ''
primary['filtered'] = spikeIn['filtered'] = ''
primary['deduped'] = spikeIn['deduped'] = ''
primary['dup_rate'] = spikeIn['dup_rate'] = ''


for i in range(0,len(primary)):
	library = primary.iloc[i]['Library ID']
	target = primary.iloc[i]['Target']
	cell = primary.iloc[i]['Cell']
	if (target == 'INPUT'):
		continue
	primary_qc_path = glob.glob(basedir + "/" + target + "/" + cell + "*" + "/**/" + "align_primary" + "/**/" + 'call-qc_report/execution/qc.json', recursive=True)
	spikeIn_qc_path = glob.glob(basedir + "/" + target + "/" + cell + "*" + "/**/" + "spike_in" + "/**/" + 'call-qc_report/execution/qc.json', recursive=True)
	with open(primary_qc_path[0]) as f:
		primary_qc = json.load(f)
	with open(spikeIn_qc_path[0]) as f:
		spikeIn_qc = json.load(f)
	primary.loc[primary.index[i], 'total_reads'] = primary_qc['align']['samstat']['rep1']['total_reads']/2
	primary.loc[primary.index[i], 'mapped'] = primary_qc['align']['samstat']['rep1']['mapped_reads']/2
	primary.loc[primary.index[i], 'proper_mapped'] = primary_qc['align']['samstat']['rep1']['properly_paired_reads']/2
	primary.loc[primary.index[i], 'filtered'] = primary_qc['align']['dup']['rep1']['paired_reads']
	primary.loc[primary.index[i], 'deduped'] = primary_qc['align']['nodup_samstat']['rep1']['total_reads']/2
	primary.loc[primary.index[i], 'dup_rate'] = primary_qc['align']['dup']['rep1']['pct_duplicate_reads']
	spikeIn.loc[spikeIn.index[i], 'total_reads'] = spikeIn_qc['align']['samstat']['rep1']['total_reads']/2
	spikeIn.loc[spikeIn.index[i], 'mapped'] = spikeIn_qc['align']['samstat']['rep1']['mapped_reads']/2
	spikeIn.loc[spikeIn.index[i], 'proper_mapped'] = spikeIn_qc['align']['samstat']['rep1']['properly_paired_reads']/2
	spikeIn.loc[spikeIn.index[i], 'filtered'] = spikeIn_qc['align']['dup']['rep1']['paired_reads']
	spikeIn.loc[spikeIn.index[i], 'deduped'] = spikeIn_qc['align']['nodup_samstat']['rep1']['total_reads']/2
	spikeIn.loc[spikeIn.index[i], 'dup_rate'] = spikeIn_qc['align']['dup']['rep1']['pct_duplicate_reads']
	if any(re.search(library, x,re.IGNORECASE) for x in input_data['Sample']):
		input=input_data[input_data['Sample'] == library]['Input'].to_string(index=False).strip()
		orig_indx = [ key for key, val in enumerate(primary['Library ID']) if val == input ][0]
		primary.loc[primary.index[orig_indx], 'total_reads'] = primary_qc['align']['samstat']['ctl1']['total_reads']
		primary.loc[primary.index[orig_indx], 'mapped'] = primary_qc['align']['samstat']['ctl1']['mapped_reads']
		primary.loc[primary.index[orig_indx], 'proper_mapped'] = primary_qc['align']['samstat']['ctl1']['properly_paired_reads']
		primary.loc[primary.index[orig_indx], 'filtered'] = primary_qc['align']['dup']['ctl1']['paired_reads']*2
		primary.loc[primary.index[orig_indx], 'deduped'] = primary_qc['align']['nodup_samstat']['ctl1']['total_reads']
		primary.loc[primary.index[orig_indx], 'dup_rate'] = primary_qc['align']['dup']['ctl1']['pct_duplicate_reads']
		spikeIn.loc[spikeIn.index[orig_indx], 'total_reads'] = spikeIn_qc['align']['samstat']['ctl1']['total_reads']
		spikeIn.loc[spikeIn.index[orig_indx], 'mapped'] = spikeIn_qc['align']['samstat']['ctl1']['mapped_reads']
		spikeIn.loc[spikeIn.index[orig_indx], 'proper_mapped'] = spikeIn_qc['align']['samstat']['ctl1']['properly_paired_reads']
		spikeIn.loc[spikeIn.index[orig_indx], 'filtered'] = spikeIn_qc['align']['dup']['ctl1']['paired_reads']*2
		spikeIn.loc[spikeIn.index[orig_indx], 'deduped'] = spikeIn_qc['align']['nodup_samstat']['ctl1']['total_reads']
		spikeIn.loc[primary.index[orig_indx], 'dup_rate'] = spikeIn_qc['align']['dup']['ctl1']['pct_duplicate_reads']



primary.to_csv("primary_alignment_stats.csv", index=False)
spikeIn.to_csv("SpikeIn_alignment_stats.csv", index=False)

#get normalizing factor
factor = min(spikeIn['deduped'])/spikeIn['deduped']
downsample = pd.DataFrame({'sample':primary['Library ID'], 'factor':factor})
downsample.to_csv("downsample_percentage.txt", sep="\t", index=False, header=False)