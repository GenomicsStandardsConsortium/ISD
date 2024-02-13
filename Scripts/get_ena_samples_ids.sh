#!/bin/bash -l
#
# Name: get_ena_sample_ids.sh
# Purpuse: get sample ids from project ids from ENA
# Author: Savvas Paragkamian

## example of running the script
## ./scripts/get_ena_samples_ids.sh -f example_projects.txt -d ena_samples

## Usage of the script
usage="Use the parameter -f for the path of txt file (inside the repository, \
and first field with the project ids) and -d for the name of the new \
directory that the results will be saved in.\n\
Example: ./scripts/get_ena_sample_ids.sh -f example_projects.txt -d ena_samples \n"

## User input PDF file
while getopts "f:d:" option
do
   case "$option" in
      f)   projects_file="${OPTARG}";;
      d)   directory="${OPTARG}";;
      ?|:)   echo -e "$usage" ; exit 1;;
      *)   echo -e "option ${OPTARG} unknown. \n$usage" ; exit 1 ;;
   esac
done

## Detect if no options were passed
if [ $OPTIND -eq 1 ]; 
    then echo -e "No options were passed. \n$usage "; exit 1; 
fi

## Detect if a single option is missing
if [ -z "${projects_file}" ]; then
    echo -e "Option -f empty. $usage"; exit 1;
fi

if [ -z "${directory}" ]; then
    echo -e "Option -d empty. $usage"; exit 1;
fi


echo $projects_file
echo $directory

##################################### Start of script ##########################
# keep the first column of the file with the project ids

runs_accession=`cut -f 1 $projects_file`

# ENA API url prefix and suffix
ena_pre='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='
ena_suf='&result=read_run&fields=all&format=tsv'

# for each id build the url and make the api call and save to a file
cd $directory

for i in $runs_accession
do
    echo "proceeding to download " $i 
    mkdir $i
    mkdir $i/sequences
    mkdir $i/ena_samples_attr
    metadata="${i}/ena_samples-${i}.txt"
    ena_call="$ena_pre""$i""$ena_suf"
    
    curl -o $metadata $ena_call
    echo "project sample ids retrieved"
    sleep 5

    # get the xml files
    echo "project sample attributes retrieving"
    ../scripts/get_ena_samples_attributes.py $metadata $i/ena_samples_attr

    # transform the xml files
    ../scripts/ena_xml_to_csv.py $i/ena_samples_attr $i/ena_samples_attributes-$i.tsv

    
    sleep 2

done

echo -e "text from $projects_file is in $directory \n"
