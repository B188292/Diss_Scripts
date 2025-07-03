#!/bin/python3

import subprocess, os, shutil

test_folder = "/scratch1/msc_2025/s2103976/dcm/GSE269705/"
accession_file = "/scratch1/msc_2025/s2103976/dcm/GSE269705.txt"
working_directory = "/scratch1/msc_2025/s2103976/"

#if os.path.exists(test_folder):
#    os.rmdir(test_folder)
#    os.mkdir(test_folder)
#else:
#    os.mkdir(test_folder)

#This code has been adapted from https://rpubs.com/snijesh/sratool-kit-data-retrieval
with open(accession_file) as test_file:
    for accession in test_file:
        #remove any whitespace from the accession numbers just in case
        accession = accession.strip()
        #download the .sra file
        subprocess.run(["prefetch", accession], check= True)
        print(f"{accession} has been downloaded from GEO")
        sra_file = os.path.join(accession, f"{accession}.sra")
        #split and compress the data. --split-3 is used just in case there are single-end reads.
        subprocess.run(["fasterq-dump", "--split-3", sra_file], check=True)
        print(f"The data from {accession} has been extracted.")

        for i in os.listdir():
            if i.startswith(accession) and i.endswith(".fastq"):
                subprocess.run(["gzip", i], check = True)
                print(f"The data from {accession} has been extracted and compressed.")

        #move the files to the corresponding folder
        for files in os.listdir(working_directory):
            if files.startswith(accession):
                shutil.move(files, test_folder)
                print(f"Moved {files} to the {test_folder}")



