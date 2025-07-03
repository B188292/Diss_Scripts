#!/bin/python3

import subprocess, os, shutil

#This is script is for SE sequences

#Creating a variable for the path of the input accessions file (replace when necessary)
se_files = "/scratch1/msc_2025/s2103976/hgps_se_acc.txt" 

#Creating a variable for the path of the directory that will contain all the trimmed files.
se_output = "/scratch1/msc_2025/s2103976/hgps_se_trimmed3/"

#noting the directory containing the files that need to be trimmed
files_to_trim = "/scratch1/msc_2025/s2103976/hgps/"

#If the output directory already exists, delete it. Else, create it.
#if os.path.isdir(se_output):
#    os.rmdir(se_output)
#    os.mkdir(se_output)
#else:
#    os.mkdir(se_output)


#Trimming the data

#Open the accessions file and note each accession
with open(se_files) as se_to_trim:
    for accession in se_to_trim:
        accession = accession.strip()
        
        filename = f"{accession}.fastq.gz"

        if not os.path.exists(os.path.join(files_to_trim, filename)):
            print(f"No corresponding file found for {accession}")
            break
        else:
            input_file = os.path.join(files_to_trim, filename)
            output_file = os.path.join(se_output, f"{accession}_trimmed.fastq.gz")
            subprocess.run(["TrimmomaticSE", "-phred33", input_file, output_file, "HEADCROP:25", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:20"], check = True)
