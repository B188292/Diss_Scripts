#!/bin/python3

import subprocess, os, shutil

#This is script is for SE sequences

#Creating a variable for the path of the input accessions file (replace when necessary)
pe_files = "/scratch1/msc_2025/s2103976/dcm/GSE120836.txt"

#Creating a variable for the path of the directory that will contain all the trimmed files.
pe_output = "/scratch1/msc_2025/s2103976/dcm/trimmed_GSE120836/"


#if os.path.isdir(pe_output):
#    os.rmdir(pe_output)
#    os.mkdir(pe_output)
#else:
#    os.mkdir(pe_output)

paired_output = os.path.join(pe_output, "paired2")
os.makedirs(paired_output, exist_ok=True)

unpaired_output = os.path.join(pe_output, "unpaired2")
os.makedirs(unpaired_output, exist_ok=True)

#noting the directory containing the files that need to be trimmed
files_to_trim = "/scratch1/msc_2025/s2103976/dcm/GSE120836/"

#Trimming the data

#Open the accessions file and note each accession
with open(pe_files) as to_trim:
    for accession in to_trim:
        accession = accession.strip()
        
        fwd = f"{accession}_1.fastq.gz"
        rev = f"{accession}_2.fastq.gz"

        #checking that the files exist
        if not os.path.exists(os.path.join(files_to_trim, fwd)) or not os.path.exists(os.path.join(files_to_trim, rev)):
            print(f"No appropriate files found for {accession}")
            break
        else:
            fwd_in = os.path.join(files_to_trim, fwd)
            rev_in = os.path.join(files_to_trim, rev)

            #creating the necessary files
            fwd_paired = os.path.join(paired_output, f"{accession}_1_paired_trimmed.fastq.gz")
            fwd_unpaired = os.path.join(unpaired_output , f"{accession}_1_unpaired_trimmed.fastq.gz")
            rev_paired = os.path.join(paired_output, f"{accession}_2_paired_trimmed.fastq.gz")
            rev_unpaired = os.path.join(unpaired_output, f"{accession}_2_unpaired_trimmed.fastq.gz")
                
            #output_file = os.path.join(pe_output, trimmed_name)
            subprocess.run(["TrimmomaticPE", "-phred33", fwd_in, rev_in, fwd_paired, fwd_unpaired, rev_paired, rev_unpaired, "HEADCROP:20", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:30"], check = True) 
            print(f"Trimmed data for {accession}")
