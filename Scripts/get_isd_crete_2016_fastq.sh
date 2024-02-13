#!/bin/bash

# Name: get_isd_crete_2016_fastq.sh
# Purpuse: get the ISD Crete 2016 project PRJEB21776 sequences
# from ENA ftp service
# Author: Savvas Paragkamian

# load the metadata file which was downloaded from the ENA website
# and keep the ftp_fasta field, $7 that has the url of each fasta file.
# In total there are 280 .fastq.gz files.

runs_accession=`gawk -F"\t" '(NR>1){split($7,fastq, ";") ; for (i in fastq){print fastq[i]}}' \
    ~/isd-crete/ena_metadata/filereport_read_run_PRJEB21776_tsv.txt`

# go to directory
cd ~/isd-crete/ena_data

# for each ftp link wget and append to a log file
for i in $runs_accession
do
    echo "proceeding to download " $i 
    wget -c --append-output=../ena_metadata/get_fastq.log $i
    sleep 3
done

echo "finish the ISD 2016 retrieval"


