# Island Sampling Day (ISD)
## Chapter one: Crete

ISD is a soil microbiome study to explore environmental drivers of community
composition across ecosystems in the island of Crete, Greece.

A snapshot of the prokaryotic diversity across the island of Crete. Aim to:
<ol type="a">
  <li>describe the first full island microbial community assessment</li>
  <li>shed light on phylotypes found not to be all the same across Crete’s distinct terrestrial habitats (<i>intra comparison</i>)</li>
  <li>present a “local” (island) vs “global” (world, e.g. EMP) soil phylotype comparative analysis (<i>inter comparison</i>)</li>
</ol>
This repository holds ISD molecular ecology analysis pertinent scripts.
ISD sequence data are available as a <a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21776" target="blank">European Nucleotide Archive repository project</a>.

## Contents

* [Methods](#methods)
* [Scripts](#scripts)
* [Data integration](#data-integration)
* [Data retrieval](#data-retrieval)
* [Sequences](#sequences)
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

## Data integration

Crete has been sampled multiple times for its' environmental microbiome.

The ISD project is the first one on this scale.


## Sequences
Total reads (forward and reverse) = 121232490 in 140 samples

Primers used are FWD: 5'-ACTCCTACGGGAGGCAGCAG-3' REV: 5'-GGACTACHVGGGTWTCTAAT-3'

Numbers of Ns in reads : there are many reads with Ns and the `dada` function
doesn't accept so after the filtering they are all removed.

## Hardware

Most computations were performed in the Zorbas HPC facility of [IMBBC-HCMR](https://hpc.hcmr.gr),
see here for more [info](https://doi.org/10.1093/gigascience/giab053).

## Licence

GNU GPLv3 license (for 3rd party scripts separate licenses apply).
