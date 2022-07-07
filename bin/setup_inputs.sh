#!/bin/bash -x

set -euo pipefail

GENOME_ID=${1}
GENOME_DIR=${GENOME_ID}

GFF_FILE=${GENOME_DIR}/${GENOME_ID}.gff
PFAM_FILE=${GENOME_DIR}/${GENOME_ID}.pfam.tab.txt

SPLIT_GFF_FILE=${GENOME_DIR}/${GENOME_ID}.gff.split
SORTED_GFF_FILE=${GENOME_DIR}/${GENOME_ID}.gff.sorted

if [[ -s ${GFF_FILE} ]]; then
    if [[ $(grep "##FASTA" ${GFF_FILE}) ]]; then
        echo "This GFF file contained fasta as well. Splitting ..."
        LINE_BEFORE_FASTA=$(awk '/##FASTA/ {print FNR-1}' ${GFF_FILE})
        head -n ${LINE_BEFORE_FASTA} ${GFF_FILE} > ${SPLIT_GFF_FILE}
    else
        cp ${GFF_FILE} ${SPLIT_GFF_FILE}
    fi
    echo "Sorting the GFF file by contig, start and end coordinates"
    sort -k1,1 -k4,5n ${SPLIT_GFF_FILE} > ${SORTED_GFF_FILE}
else
    echo "Expected GFF file not found!"
    exit 1
fi

NUM_PFAM_FILES=$(find ${GENOME_DIR}/ -name "*pfam*" | wc -l | awk '{print $1}')

if [[ -s ${PFAM_FILE} ]]; then
    echo "Found expected PFam file."
elif [[ ${NUM_PFAM_FILES} -gt 1 ]]; then
    echo "This Genome seems to have ${NUM_PFAM_FILES} PFam files. None of which match the expected naming convention. Exiting."
    exit 1
elif [[ ${NUM_PFAM_FILES} -eq 0 ]]; then
    echo "This Genome does not seem to have a PFam file. Exiting."
    exit 1
elif [[ ${NUM_PFAM_FILES} -eq 1 ]]; then
    echo "Found 1 PFam file that does not match the expected naming convention. Will use this PFam file."
    FOUND_PFAM=$(find ${GENOME_DIR}/ -name "*pfam*")
    mv ${FOUND_PFAM} ${PFAM_FILE}
fi
