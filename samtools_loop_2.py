#!/bin/python3

import shutil, subprocess, os
accession_list = "/scratch1/msc_2025/s2103976/hgps_se_accessions.txt"
bam_location = "/scratch1/msc_2025/s2103976/hgps_se_aligned/"
indexed_bam_loc = "/scratch1/msc_2025/s2103976/hgps_se_indexed/"

with open(accession_list) as accessions:
    for name in accessions:
        name = name.strip()
        print(f"Processing {name}")

        for bam in os.listdir(bam_location):
                if bam.startswith(name) and bam.endswith(".bam"):
                    print(f"{name} has a corresponding BAM file.")
                    # Define file paths for sorted and indexed BAM files
                    input_bam = os.path.join(bam_location, bam)
                    sorted_bam = os.path.join(bam_location, f"{name}_sorted_.bam")
                    indexed_bam = os.path.join(indexed_bam_loc, f"{name}_sorted.bam")
                
                    # Run samtools sort to sort the BAM file
                    print("Sorting the BAM file")
                    subprocess.run(f"samtools sort {input_bam} -o {sorted_bam}", shell=True)

                    # Run samtools index to index the sorted BAM file
                    print("Indexing the BAM file.")
                    subprocess.run(f"samtools index {sorted_bam}", shell=True)
    
                    #Move the sorted BAM to the final directory
                    subprocess.run(f"mv {sorted_bam} {indexed_bam}", shell=True)

                    print(f"Processed {bam} -> sorted and indexed to {indexed_bam}")
