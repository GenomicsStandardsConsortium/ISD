#!/usr/bin/env python3

###############################################################################
# script name: get_isd_crete_2016_attributes.py
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to retrieve all available attributes from 
# ENA assigned to each sample of an ENA project.
###############################################################################
# usage:./get_isd_crete_2016_attributes.py path/to/ena_project_file.tsv ena_project_directory
###############################################################################

import requests, sys, time
from random import uniform
from datetime import date


def attr_request(accession):

    ena_url = "https://www.ebi.ac.uk/ena/browser/api/xml/"
    url = ena_url + accession + "?download=true&includeLinks=true"

    start_time = time.time()
    connection_timeout = 60 # seconds

    while True:
        try:
            attr = requests.get(url = url)
            break
        except requests.ConnectionError:
            if time.time() > start_time + connection_timeout:
                raise Exception('Unable to get updates after {} seconds of ConnectionErrors'.format(connection_timeout))
            else:
                time.sleep(1)

    code = attr.status_code

    if code == 204:
        print(str(accession) + "\t" + str(code))

    if code == 400:
        print(str(accession) + "\t" + str(code))

    if code == 200:
        print(str(accession) + "\t" + str(code))
        attributes = attr.text
        return(attributes)

    time.sleep(uniform(1,2))

def write_xml(attributes,accession, ena_project_directory):
    
    # load to xml
    newxml = str(ena_project_directory) + "/" + "attributes_" + str(accession) + ".xml"
    
    with open(newxml, "w") as xmlFile:
        xmlFile.write(attributes)

############################# main code ######################################

# user input
ena_project_file = sys.argv[1]
ena_project_directory = sys.argv[2]

samples = []

with open(ena_project_file) as file:
    next(file) #skip the first line
    for line in file:
        
        l=line.split('\t')
        samples.append(l)

print("there are " + str(len(samples)) + " samples accession ids")

## Iterate and make the GET requests to server

for a in samples:
    accession = a[0]
    
    attributes = attr_request(accession)
    if attributes != None:
        write_xml(attributes,accession,ena_project_directory)

    else:
        continue

