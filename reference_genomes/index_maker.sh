#!/bin/bash

#Creating the human index
human_index="/scratch1/msc_2025/s2103976/reference_genomes/h_index_59/"
human_fasta="/scratch1/msc_2025/s2103976/reference_genomes/GCF_000001405.40_GRCh38.p14_genomic.fna"
human_annot="/scratch1/msc_2025/s2103976/reference_genomes/GCF_000001405.40_GRCh38.p14_genomic.gtf"

echo "Creating the human reference genome..."
STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir ${human_index} \
	--genomeFastaFiles ${human_fasta} \
	--sjdbGTFfile ${human_annot} \
	--sjdbOverhang 59
echo "Created the human reference genome."

#mouse_index="/scratch1/msc_2025/s2103976/reference_genomes/mouse_index/"
#mouse_fasta="/localdisk/home/s2103976/dissertation/reference_genomes/GCF_000001635.27_GRCm39_genomic.fna"
#mouse_annot="/localdisk/home/s2103976/dissertation/reference_genomes/GCF_000001635.27_GRCm39_genomic.gtf"

#echo "Creating the mouse reference genome..."
#STAR --runThreadN 10 \
#        --runMode genomeGenerate \
#        --genomeDir ${mouse_index} \
#        --genomeFastaFiles ${mouse_fasta} \
#        --sjdbGTFfile ${mouse_annot} \
#        --sjdbOverhang 100
#echo "Created the mouse reference genome."
