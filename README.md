# Prolonged_persistence_of_DNA_lesions
Code accompanying the manuscript 'Prolonged persistence of mutagenic DNA lesions in stem cells'

# General notes
The code is provided for all analyses starting from a mutation list with read count matrices, and a phylogeny.  The mutation calling and tree inference is not covered here and is included within the original manuscripts. \
The numbered folders go through the different steps in the analysis, roughly in the order in which they were performed. However, (most) intermediate data is available such that each stage of the analysis can be readily performed without having to re-run earlier stages. \
Smaller data objects are provided within the github repository. Larger data objects need to be downloaded from Mendeley Data at  doi: 10.17632/9tw3kbj2cw.1 \
If there are any queries, feel free to contact me via email: ms56@sanger.ac.uk

## Objects that need to be downloaded from Mendeley Data
data/input_data \
data/simulation_results/MAV_sim_results.tsv

# Notes on specific stages of data generation

## 01_Running the core algorithm
This includes the core analysis script 'Detect_persistent_lesions.R' which is designed to be set off via the command line.  There are several options to point to the necessary input files, and flags to force re-analysis, remove duplicates, and to include a dummy ancestral branch (should be used if there is such a branch in the provided tree).
The file Setting_off_persistent_lesion_analysis.sh contains the code to set off the analysis script for the different datasets in the study.
If trying to run on your own data you may need to edit the script to account for the format of your data.

## 02_Phasing and LOH analysis
The scripts provided here show how the phasing and LOH analyses were performed.
There is a core phasing script underlying this - this is run in julia and will require considerable tweaking to get functioning locally as there are multiple fixed file paths.

## 03 Comipiling fitered mutation list
This is a single script that takes the data from the previous two stages, as well as a comprehensive list of sample metadata to generate the final set of data objects for all downstream analyses.

## 04 Simulation_scripts
A selection of scripts that run various simulations:
- MAV_simulations_2.R; script to model the NULL model of multi-allelic variants occurring from independent mutation acquisition
- PVV_simulation_2.R; script to model the NULL model of phylogeny-violating variants occurring from independent mutation acquisition
- PVV_simulation_reversion_2.R; script to model somatic reversion events as a possible cause of phylogeny-violating variants
- Tree_structure_PVV_sensitivity_2.R; script to model the proportion of introduced lesions that would be anticipated to be detected given the phylogeny structure

Note that the output from these scripts is provided in the github repo in Data/simulation_results

## 05 ABC (Approximate Bayesian Computation)
All the scripts required to run the ABC - these need to be run in the correct order (as below). Lesion_duration_ABC_new.sh contains the commands to set off the jobs on the command line.

1. generate_simpops.R - each run of this script will generate a simulated population like the 40 used in the ABC. Run with the index of the simulation as a trailing argument i.e. Rscript generate_simpops.R 13
2. ABC_simulation_new_INTRODUCE_LESIONS.R - introduce lesions into each simulated population. This takes two trailing arguments (1) the index of the simpop to introduce lesions into and (2) the mean duration of the introduced lesions (in years).
3. ABC_simulation_new_INTRODUCE_LESIONS_BOOST.R - a script to see if there are sufficient PVVs generated from the first run of lesion generation, and introduce additional lesions if not.
4. ABC_simulation_new_parallel.R - script for the main simulation of the ABC.

## 06 Generating Figures
Individual scripts to generate all the figures/ extended data figures from the manuscript. All the data for the figures is available from the github/ Mendeley data without having the re-run all the previous steps.


# OTHER FOLDERS
## Original analysis scripts
These are the original scripts used to analyse the data. They are not quite as 'tidy' as those within the '06_Generating_figures/' folder, and much of this is duplications of the same analyses. However, some scripts contain code for some additional plots and analyses that are not included in the manuscript and are therefore included for reference.

## Miscellaneous scripts

## Plots
This is where all plots from the figure generation scripts are saved.

## Data
Data/reference_files/ - includes reference files used in various analyses \
Data/metadata_files/ - includes individual-level metadata, and sample-level metadata \
Data/VCFs/ - includes mutation VCFs generated in/ required for various mutational signature analyses \
Data/ASCAT_LOH_analysis/ - results of the ASCAT and LOH scripts in the folder 02_Phasing and LOH analysis/ \
Data/phasing_results/ - results of the phasing scripts in the folder 02_Phasing and LOH analysis/ \
Data/lesion_segregation_analysis/ - includes intermediate data from lesion segregation analysis as this is very time-consuming \
Data/ABC_results/ - includes intermediate data from the ABC simulations to allow the ABC to be re-run and the plots generated\
