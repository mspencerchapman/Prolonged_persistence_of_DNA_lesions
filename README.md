# Prolonged_persistence_of_DNA_lesions
Code accompanying the manuscript 'Prolonged persistence of mutagenic DNA lesions in stem cells'

## Table of Contents

- [General notes](#general-notes)
- [System requirements](#system-requirements)
  - [Objects that need to be downloaded from Mendeley Data](#objects-that-need-to-be-downloaded-from-mendeley-data)
- [Notes on specific stages of data generation](#notes-on-specific-stages-of-data-generation)
  - [01_Running the core algorithm](#01_running-the-core-algorithm)
  - [02_Phasing and LOH analysis](#02_phasing-and-loh-analysis)
  - [03_Compiling fitered mutation list](#03_compiling-fitered-mutation-list)
  - [04_Simulation_scripts](#04_simulation_scripts)
  - [05_ABC (Approximate Bayesian Computation)](#05_abc-(approximate-bayesian-computation))
    - [Example code for generating the simulated populations](#example-code-for-generating-the-simulated-populations)
  - [06_Generating Figures](#06_generating-figures)
- [OTHER FOLDERS](#other-folders)
  - [Original analysis scripts](#original-analysis-scripts)
  - [Miscellaneous scripts](#miscellaneous-scripts)
  - [Plots](#plots)
  - [Data](#data)

# General notes
The code is provided for all analyses starting from a mutation list with read count matrices, and a phylogeny.  The mutation calling and tree inference is not covered here and is included within the original manuscripts. \
The numbered folders go through the different steps in the analysis, roughly in the order in which they were performed. However, (most) intermediate data is available such that each stage of the analysis can be readily performed without having to re-run earlier stages. \
In most scripts you will need to adjust some file paths to reflect those on your own system \
e.g. path to the cloned github repository, or the path to an available genome file (in the analyses this is GRCh37, ensembl version i.e. no 'chr' in the chromosome names) \
Smaller data objects are provided within the github repository. Larger data objects need to be downloaded from Mendeley Data at  https://data.mendeley.com/datasets/9tw3kbj2cw/1 \
If there are any queries, feel free to contact me via email: ms56@sanger.ac.uk

# System requirements
Much of the code can be run on a standard personal computer using R v4.4.0. There are several packages used that are listed in each individual script, all of which are available either via CRAN, BioConductor, or github. Most are installed up within the scripts if not already installed.
Some of the steps need a high performance computing (HPC) cluster as outlined below. Some of the code is designed to use an HPC cluster using the LSF system. For alternative job submission systems these will need to be altered.
The code has been tested on linux Ubuntu 22, and Mac OS 14.5 Sonoma. No non-standard hardware is required.

## Objects that need to be downloaded from Mendeley Data
Download and place the objects within the cloned repository using the data structure below.

Data/input_data \
Data/simulation_results/MAV_sim_results.tsv

# Notes on specific stages of data generation

## 01_Running the core algorithm
This includes the core analysis script 'Detect_persistent_lesions.R' which is designed to be set off via the command line.  There are several options to point to the necessary input files, and flags to force re-analysis, remove duplicates, and to include a dummy ancestral branch (see below).

| Option (short) | Option (long) | Description |
|----|-----|-----------------------|
|-w|--working_dir|Working directory for analysis - should be set to the cloned github directory. If not set will default to the current directory|
|-s|--sample_id|Sample ID to include in output file names|
|-t|--tree_path|path for the tree file - mandatory|
|-f|--filtering_muts_path|path for the filtered_muts file - mandatory, and needs to be in a very specific structure|
|-f|--genome_file|path to your local genome file (needs to be GRCh37, ensemble version)|
|-a|--ancestral|Set this flag if input trees contain an ancestral tip|
|-d|--duplicate_remove|Remove any duplicate samples that are present (i.e. samples that come from the same clone)|
|-r|--re_run|Force rerun of TREEMUT and the PVV assessment, even if conditions are met for re-importing existing intermediate files|
|-c|--cores|Number of cores to use for parallel steps|
|-o|--output_dir|output directory for files|

This bash code snippet can be used to set off the script locally to analyse the 3 clonal haematopoiesis individuals, as this can be done in a short space of time using the resources available on most personal laptop/ desktop computers. Note that the data/input_data directory must be downloaded from Mendeley Data for this to run.

```bash
#Set file paths - the STUDY_DIR should be the cloned github directory
STUDY_DIR='PATH/TO/CLONED/REPOSITORY/'
GENOME_FILE='GRCh37/GENOME/FILE/PATH/genome.fa'

#Define the secondary file paths (these are as follows if the github directory has been cloned)
#Note that the 'input_data' directory needs to be downloaded from Mendeley Data as per the README file
LESION_SEG_INPUT_DIR="${STUDY_DIR}Data/input_data/"
MAIN_ANALYSIS_SCRIPT_PATH=${STUDY_DIR}"01_Running_the_core_algorithm/Detect_persistent_lesions.R"
OUTPUT_DIR=${STUDY_DIR}output/

setopt shwordsplit #Only relevant if using zsh on mac
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
DATASET=MF
for EXP_ID in $SAMPLES; do
FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_noMixed.tree
Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done

```

The file Setting_off_persistent_lesion_analysis.sh contains the code to set off the 'Detect_persistent_lesions.R'script on an LSF compute farm for the different datasets in the study, including information for approximately the amount of memory required for each dataset. <strong>Copy and paste the lines individually rather than run the whole script.</strong> \

If trying to run on your own data you may need to edit the script to account for the format of your data.
Smaller datasets (e.g. trees of <100-150 samples) can be run locally on a personal computer.

Note that the analysis of the liver MAVs is done outside of this core script as the original samples are non-clonal and therefore much of the approach is bespoke and simplified. Therefore this is in the separate script 'Liver_MAV_analysis.R'.

## 02_Phasing and LOH analysis
The scripts provided here show how the phasing and LOH analyses were performed.
There is a core phasing script underlying this - this is run in julia and will require considerable tweaking to get functioning locally as there are multiple fixed file paths to large genome reference files.

## 03_Compiling fitered mutation list
This is a single script that takes the data from the previous two stages, as well as a comprehensive list of sample metadata to generate the final set of data objects for all downstream analyses.

## 04_Simulation_scripts
A selection of scripts that run various simulations using the phylogenies from the data as their starting point. Most are calculating the expected distributions of MAVs/ PVVs from alternative mechanisms (i.e. not from a persistent DNA lesion)

|Script name|Description|
|-----|-------------------------------------|
|MAV_simulations_2.R| NULL model of multi-allelic variants occurring from independent mutation acquisition|
|PVV_simulation_2.R| NULL model of phylogeny-violating variants occurring from independent mutation acquisition|
|PVV_simulation_reversion_2.R|Models somatic reversion events as a possible cause of phylogeny-violating variants|
|Tree_structure_PVV_sensitivity_2.R| Models the proportion of introduced lesions that would be anticipated to be detected given the phylogeny structure|

Note that the output from these scripts is provided in the github repo in Data/simulation_results

## 05_ABC (Approximate Bayesian Computation)
All the scripts required to run the ABC - these need to be run in the correct order (as below).
CONTROLLING_THE_ABC.sh contains the commands to set off the jobs on the command line (using LSF) in the correct order. <strong>Copy and paste the lines individually rather than run the whole script.</strong> \

|Script name|Description|
|-----|-------------------------------------|
|generate_simpops.R| Each run of this script will generate a simulated complete haematopoietic population for the ABC using the R package 'rsimpop'. 40 are necessary for the ABC. Run with the index of the population as a trailing argument i.e. Rscript generate_simpops.R 13 |
|ABC_simulation_new_INTRODUCE_LESIONS.R|Introduce lesions into each simulated population. This takes two trailing arguments (1) the index of the simpop to introduce lesions into and (2) the mean duration of the introduced lesions (in years)|
|ABC_simulation_new_INTRODUCE_LESIONS_BOOST.R|A script to see if there are sufficient PVVs generated from the first run of lesion generation, and introduce additional lesions if not.
|ABC_simulation_new_parallel.R|Script for the main simulation of the ABC. Designed to run with 4 cores with analysis for the four phylogenies running in parallelel. Each run of the simulation typically takes ~24hours.
|ABC_new_extract_info.R|Code to extract the pertinent lesion info from each simulation, generate the summary statistics and then run the ABC. Easiest to run interactively. This is mostly duplicated in the script: 'Generate_Fig5a.R' and ' Generate_ExtDatFig11.R'|
|ABC_new_GROUND_TRUTH.R|Extracts the data from the 'ground truth' simulations to assess the accuracy of the approach|
|ABC_new_PPCs_extract_info.R|Extracts the data from the 'posterior predictive checks' simulations to assess the degree to which the posterior accurately emulates the summary statistics of the data.|
|ABC_simulation_new_PPCs.R|Script for running the PPC simulations (identical to the main simulation run but with different output directory!)|

In addition there are some bash wrapper scripts (lesion_boost_wrapper.sh AND PPC_wrapper.sh) to allow various steps to be run as an LSF array.

### Example code for generating the simulated populations
```bash
#Go into the cloned repository, as script contains relative file paths
cd PATH/TO/CLONED/REPOSITORY/

#Set off script to generate simulated population number 1 (it will be saved with index 1)
#This will typically take 30-60 minutes, therefore best to do in parallel on a compute farm
Rscript 05_ABC/generate_simpops.R 1

```

## 06_Generating Figures
Individual scripts to generate all the figures/ extended data figures from the manuscript, divided by into folders by figure. All the data for the figures is available from the github/ Mendeley data without having the re-run all the previous steps.

# OTHER FOLDERS
## Original analysis scripts
These are the original scripts used to analyse the data. They are not quite as 'tidy' as those within the '06_Generating_figures/' folder, and much of this is duplications of the same analyses. However, some scripts contain code for some additional plots and statistical analyses that are not included in the manuscript and are therefore included for reference.

## Miscellaneous scripts
Some additional miscellaneous scripts:
Running_HDP/ - scripts for running the HDP signature extraction
Manual_fix_for_KX003_tree.R - short script for altering the phylogeny configuration
Bootstrap_tree_analysis_readcounts.R - script for analysing the bootstrap readcounts
Running_IQtree_and_SCITE.R - script for setting off IQtree and SCITE phylogeny inference algorithms and comparing them to the MPBoot phylogenies.

## Plots
This is where all plots from the figure generation scripts are saved.

## Data
Contains various forms of data required for the analysis.
Some data is from reference files required for the analysis. Some are intermediate files generated during the analysis but which are time-consuming to run, therefore results are provided.

|Directory| Description|
|------|-------------------------------|
|Data/reference_files/ | Reference files used in various analyses (does not include genome files) |
|Data/metadata_files/ | Individual-level metadata, and sample-level metadata |
|Data/VCFs/ | Mutation VCFs generated in/ required for various mutational signature analyses |
|Data/ASCAT_LOH_analysis/ | Results of the ASCAT and LOH scripts in the folder 02_Phasing and LOH analysis/ |
|Data/phasing_results/ | Results of the phasing scripts in the folder 02_Phasing and LOH analysis/ |
|Data/lesion_segregation_analysis/ | Includes intermediate data from lesion segregation analysis as this is very time-consuming to run from scratch. |
|Data/ABC_results/ |  Intermediate data from the ABC simulations to allow the ABC to be re-run and the plots generated|
