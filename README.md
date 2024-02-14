# Island Sampling Day (ISD)
## Chapter one: Crete

ISD is a soil microbiome study to explore environmental drivers of community
composition across ecosystems in the island of Crete, Greece.

A snapshot of the prokaryotic diversity across the island of Crete. Aim to:

1. Put metadata into action through a citizen science project
2. describe the first full island microbial community assessment
3. shed light on phylotypes found not to be all the same across Crete’s distinct terrestrial habitats (<i>intra comparison</i>)
4. present a “local” (island) vs “global” (world, e.g. EMP) soil phylotype comparative analysis (<i>inter comparison</i>)

This repository holds ISD molecular ecology analysis pertinent scripts.
ISD sequence data are available as a <a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21776" target="blank">European Nucleotide Archive repository project</a>.

## Contents

* [Methods](#methods)
* [Scripts](#scripts)
* [Data retrieval](#data-retrieval)
* [Data integration](#data-integration)
* [Inference and Taxonomy](#inference-and-taxonomy)
* [Analysis](#analysis)
* [Software](#software)
* [Hardware](#hardware)
* [Citation](#citation)
* [Licence](#licence)

## Methods

## Scripts
The scripts of the analysis are in the `Scripts` folder and cover the following tasks:

* Search ENA for samples in the island of Crete
* Get ISD metadata and sequences
* HPC jobs and parameter files
* Filtering, clustering/denoising and taxonomic assigninments
* Biodiversity analysis
* Figures

## Data retrieval

### Metadata

Retrieve the filereport of ENA for the ISD Crete project id PRJEB21776

```
curl -o Data/Metadata/filereport_read_run_PRJEB21776_tsv.txt \
    'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB21776&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp,bam_ftp&format=tsv'

```

Then for each sample retrieve all attributes with the `get_isd_crete_2016_attributes.py` 
script. These are all in `xml` format. Using the `ena_xml_to_csv.py` these files 
are aggregated to an long format tsv file.

### Sequences

To retrieve the sequences use the `get_isd_crete_2016_fastq.sh` script after
downloading the metadata.

Total reads (forward and reverse) = 121232490 in 140 samples

Primers used are FWD: 5'-ACTCCTACGGGAGGCAGCAG-3' REV: 5'-GGACTACHVGGGTWTCTAAT-3'

Numbers of Ns in reads : there are many reads with Ns and the `dada` function
doesn't accept so after the filtering they are all removed.

## Data integration

Crete has been sampled multiple times for its' environmental microbiome.
The ISD project is the first one on this scale.

To investigate which samples are in Crete and are related to microbiome we 
did an analysis. 
* get all samples metadata from ENA
* filter locations to be in Crete, island
* filter the metagenomic samples that are also terrestrial and visualise

## Hardware

Most computations were performed in the Zorbas HPC facility of [IMBBC-HCMR](https://hpc.hcmr.gr),
see here for more [info](https://doi.org/10.1093/gigascience/giab053).

## Licence

GNU GPLv3 license (for 3rd party scripts separate licenses apply).
