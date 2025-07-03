#!/bin/python3

import shutil, subprocess, os

#This script is for conducting the genome alignment of SINGLE-ended reads. 

#When you're making the alignment, you need to specify the reference genome. I have saved the locations of the reference genomes as variables to make this easier. the path at the moment contains my student name, which may raise issues for 
human_ref = "/scratch1/msc_2025/s2103976/reference_genomes/GCF_000001405.40_GRCh38.p14_genomic.fna"
human_annot = "/scratch1/msc_2025/s2103976/reference_genomes/GCF_000001405.40_GRCh38.p14_genomic.gtf"
#mouse_ref = "/localdisk/home/s2103976/dissertation/reference_genomes/GCF_000001635.27_GRCm39_genomic.fna"
#mouse_annot = "/localdisk/home/s2103976/dissertation/reference_genomes/GCF_000001635.27_GRCm39_genomic.gtf"

human_index = "/scratch1/msc_2025/s2103976/reference_genomes/h_index_59/"

#Provide the path for the file containing the accession numbers
accession_list = "/scratch1/msc_2025/s2103976/hgps_se_accessions_59.txt"
#Provide the path for the directory containing the files you wish to align
directory_of_files_to_align = "/scratch1/msc_2025/s2103976/hgps_se_trimmed/"

aligned_directory = "/scratch1/msc_2025/s2103976/hgps_se_aligned/"

with open(accession_list) as files_to_read:
    for accession in files_to_read:
        accession = accession.strip()

        for doc in os.listdir(directory_of_files_to_align):
            if doc.startswith(accession) and doc.endswith("_trimmed.fastq.gz"):
                file_to_align = os.path.join(directory_of_files_to_align, doc)
                output_prefix = os.path.join(aligned_directory, f"{accession}_{doc.split('.')[0]}_")
                subprocess.run(f"STAR --runThreadN 10 --sjdbOverhang 59 --sjdbGTFfile {human_annot} --genomeDir {human_index} --readFilesIn {file_to_align} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {output_prefix}", shell = True)
