#!/usr/bin/env python3

###############################################################################
# script name: ena_xml_to_csv.py
# developed by: Savvas Paragkamian
# framework: ISD Crete 2016
###############################################################################
# GOAL:
# Aim of this script is to transform the xml available attributes from 
# ENA assigned to each sample of ISD Crete 2016 to csv.
###############################################################################
# usage:./scripts/ena_xml_to_csv.py
###############################################################################
import xml.etree.ElementTree as ET
import csv
import os,sys

# Directory containing the XML files
xml_dir = 'ena_samples_attr/'

# Output TSV file
output_file = 'ena_metadata/ena_isd_2016_attributes.tsv'

# List to store the extracted data
data = []

# Function to extract xml elements paired with a leading tag.
# for example the following children are in one list
# (geographic location (depth), 5 cm, m)
#          <SAMPLE_ATTRIBUTE>
#               <TAG>geographic location (depth)</TAG>
#               <VALUE>5 cm</VALUE>
#               <UNITS>m</UNITS>
#          </SAMPLE_ATTRIBUTE>
def groub_xml_elements(xml_list,tag):
   attribs = []
   temp_list = []
   for i in xml_list:
       
       if i.tag == tag:
           temp_list = []
           temp_list.append(i.text)
           attribs.append(temp_list)

       else :
           temp_list.append(i.text)

   return(attribs)

# Iterate over each XML file in the directory
# ENA has this structure of XML 
# ['IDENTIFIERS', 'TITLE', 'SAMPLE_NAME', 'DESCRIPTION', 'SAMPLE_LINKS', 'SAMPLE_ATTRIBUTES']
for filename in os.listdir(xml_dir):
    if filename.endswith('.xml'):
        # Parse the XML file
        
        tree = ET.parse(os.path.join(xml_dir, filename))
        root = tree.getroot()

        # here we store all the categories of the XML file
        xml_elem = []
        for sample in root.findall("./SAMPLE/"):
            xml_elem.append(sample.tag)
        
        title = []
        for titl in root.findall("./SAMPLE/TITLE"):
            title.append(titl.tag)
            title.append(titl.text)
        
        link_attributes = root.findall("./SAMPLE/SAMPLE_LINKS/SAMPLE_LINK/XREF_LINK/")
        links = groub_xml_elements(link_attributes, "DB")
        
        sample_attr = root.findall("./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/")
        attribs = groub_xml_elements(sample_attr, "TAG")

        # Extract data from XML and store in a dictionary
        xml_data = {}
        xml_data['filename'] = filename  # store the filename
        xml_data['title'] = title # store the title list
        xml_data['links'] = links # store the links list
        xml_data['attributes'] = attribs # store the sample attributes
        
        # Append the dictionary to the data list
        data.append(xml_data)

# Get the list of keys from the first dictionary in the data list
fieldnames = ['file','TAG', 'VALUE', 'UNITS']

# Write the data to TSV file
with open(output_file, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    
    # Write the header
    writer.writerow(fieldnames)
    
    # Write the data rows per dictionary/file
    for dictionary in data:
        
        #name the lists
        filename = dictionary['filename']
        title = dictionary['title']
        links = dictionary['links']
        attribs = dictionary['attributes']

        # iterate the links of the xml and append the filename to the end
        for link in links:

            writer.writerow([str(filename)] + link)

        # iterate the sample attributes of the xml and 
        # append the filename to the end
        for attrib in attribs:

            writer.writerow([str(filename)] + attrib)
            
print(f"Conversion complete. The TSV file '{output_file}' has been created.")
