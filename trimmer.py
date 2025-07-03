#!/bin/python3

import subprocess, os, shutil

#This is script is for SE sequences

#Creating a variable for the path of the input accessions file (replace when necessary)
pe_files = "/scratch1/msc_2025/s2103976/hgps_pe_accessions.txt" 

#Creating a variable for the path of the directory that will contain all the trimmed files.
pe_output = "/scratch1/msc_2025/s2103976/hgps_pe_trimmed/"

paired_output = os.makedir(os.path.join(pe_input, "paired"))
unpaired_ouput = os.makedir(os.path.join(pe_output, "unpaired"))

#noting the directory containing the files that need to be trimmed
files_to_trim = "/scratch1/msc_2025/s2103976/hgps/"

#If the ourput directory already exists, delete it. Else, create it.
if os.path.isdir(pe_output):
    os.rmdir(pe_output)
    os.mkdir(pe_output)
else:
    os.mkdir(pe_output)


#Trimming the data

#Open the accessions file and note each accession
with open(pe_files) as to_trim:
    for accession in to_trim:
        accession = accession.strip()
        
        fwd = f"{accession}_1.fastq.gz"
        rev = f"{accession}_2.fastq.gz"

        #trimming the data
        for doc in files_to_trim:
            if doc.startswith(accession) and (doc.endswith("_1.fastq.gz") or doc.endswith("_2.fastq.gz")):
                fwd_paired = os.path.join(paired_output, f"{accession}_1_paired_trimmed.fastq.gz")
                fwd_unpaired = os.path.join(unpaired_output , f"{accession}_1_unpaired_trimmed.fastq.gz")
                rev_paired = os.path.join(paired_output, f"{accession}_2_paired_trimmed.fastq.gz")
                rev_unpaired = os.path.join(unpaired_output, f"{accession}_2_paired_trimmed.fastq.gz")
                
                #output_file = os.path.join(pe_output, trimmed_name)
                subprocess.run(["TrimmomaticPE", "-phred33", fwd, rev, fwd_paired, fwd_unpaired, rev_unpaired, rev_unpaired, "HEADCROP:15", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:30"], check = True)
        else: 
            print("No appropriate files found.")
