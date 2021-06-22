import numpy as np
import pandas as pd
import os, glob
import pickle
import json
import re
import subprocess

design_sheet = pd.read_csv('../primary_alignment_stats.csv', header=0, index_col=False)
basedir='/home/gdstantonlab/mxw010/Data/Ben/MYOD/GSL-BS-2095/MYOD1'
design = design_sheet[design_sheet['Target'] != "INPUT"]

for i in range(0,len(design)): 
    sample = design.loc[design.index[i],]['Library ID']
    cell = design.loc[design.index[i],]['Cell']
    json_path = glob.glob(basedir + "/" + cell + "*" + "/**/" + 'align_primary/working.json', recursive=True)
    with open(json_path[0]) as f:
            primary_json = json.load(f)
    primary_json["chip.true_rep_only"] = 'false'
    #get the bam file:
    input_file = primary_json['chip.fastqs_rep1_R1'][0] 
    #bash command to get the SEG number from the fastq file
    bashCmd = "echo" + ' "' + input_file + '" | tr "/" "\n" | grep "SEG*"'
    bashCmd = bashCmd + ' | cut -d"_" -f1'
    p = subprocess.Popen(bashCmd, shell=True, stdout=subprocess.PIPE)
    output, error = p.communicate()
    #output is a byte, needs to get converted to string with decode()
    case_bam = glob.glob(basedir + "/bams/bam/" + output.decode().strip() + "*", recursive=True)
    #bam file for control
    input_file = primary_json['chip.ctl_fastqs_rep1_R1'][0] 
    #bash command to get the SEG number from the fastq file
    bashCmd = "echo" + ' "' + input_file + '" | tr "/" "\n" | grep "SEG*"'
    bashCmd = bashCmd + ' | cut -d"_" -f1'
    p = subprocess.Popen(bashCmd, shell=True, stdout=subprocess.PIPE)
    output, error = p.communicate()
    #output is a byte, needs to get converted to string with decode()
    control_bam = glob.glob(basedir + "/bams/bam/" + output.decode().strip() + "*", recursive=True)
    #write the bam file
    primary_json['chip.nodup_bams'] = case_bam
    primary_json['chip.ctl_nodup_bams'] = control_bam
    #get rid of fastq
    for element in primary_json.copy():
        if 'fastqs' in element or 'chip.align_only' in element:
            primary_json.pop(element)
    #write 
    file = open(json_path[0].replace("align_primary","analysis"), 'w')
    file.write(json.dumps(primary_json, indent=4, separators = (","," : "))) 
    file.close()