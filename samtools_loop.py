#!/bin/python

import shutil, subprocess, os

accession_list = "/localdisk/home/s2103976/dissertation/hgp_test_1/control_subjects.txt"
bam_location = "/localdisk/home/s2103976/dissertation/hgp_test_1/aligned_genome_control/"
indexed_bam_loc = "/localdisk/home/s2103976/dissertation/hgp_test_1/indexed_controls/"

with open(accession_list) as accessions:
    for name in accessions:
        name = name.strip()

    for bam in os.listdir(bam_location):
            if bam.startswith(name) and bam.endswith(".bam"):
                # Define file paths for sorted and indexed BAM files
                sorted_bam = f"{bam_location}/sorted_{bam}"
                indexed_bam = f"{indexed_bam_loc}/{bam.replace('.bam', '_sorted.bam')}"
                
                # Run samtools sort to sort the BAM file
                subprocess.run(f"samtools sort {bam_location}/{bam} -o {sorted_bam}", shell=True)

                # Run samtools index to index the sorted BAM file
                subprocess.run(f"samtools index {sorted_bam}", shell=True)

                # Optionally, move the sorted BAM to the final directory
                subprocess.run(f"mv {sorted_bam} {indexed_bam}", shell=True)

                print(f"Processed {bam} -> sorted and indexed to {indexed_bam}")
