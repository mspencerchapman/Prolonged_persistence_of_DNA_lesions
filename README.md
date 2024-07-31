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

## 01 Running the core algorithm
sdf

## 02 Phasing and LOH analysis
The scripts provided here allow this to be re-run locally, or better, on a compute farm with parallel computing capability. All input data & scripts are provided.

## 03 Comipiling fitered mutation list
This is a single script that takes the data from the previous two stages, as well as a comprehensive list of sample metadata to generate the final set of data objects for all downstream analyses.

## 04 Simulation_scripts
A selection of scripts that run various simulations:
- MAV_simulations_2.R; script to model the NULL model of multi-allelic variants occurring from independent mutation acquisition
- PVV_simulation_2.R; script to model the NULL model of phylogeny-violating variants occurring from independent mutation acquisition
- PVV_simulation_reversion_2.R; script to model somatic reversion events as a possible cause of phylogeny-violating variants
- Tree_structure_PVV_sensitivity_2.R; script to model the proportion of introduced lesions that would be anticipated to be detected given the phylogeny structure

## 05 ABC (Approximate Bayesian Computation)
These scripts are written in julia. They need to be run on data using the three separate baitsets individually. This is because for each baitset, the samples coming from different individuals to the one in which the mutation was called are used to estimate the locus-specific sequencing error rates.

## 06 Generating Figures
These scripts must be run on a compute farm, and use nextflow for job submission.
While they are set up for submission on LSF, they could readily be adapted to other systems.


# OTHER FOLDERS
## Original analysis scripts
These are the original scripts used to analyse the data. They are not quite as 'tidy' as those within the '05 Generating_figures/' folder, and much of this is duplications of the same analyses. However, some scripts contain code for some additional plots and analyses that are not included in the manuscript and are therefore included for reference.

## Other simulation scripts
Simulation scripts used for analyses other than the main 'engrafting cell number' ABC.
This includes:
1. estimating the phylogenetic age \
2. estimating the effect that increased T-cell clone longevity may have on clonal composition relative to the myeloid fraction

## data
data/reference_files/ - includes reference files used in various analyses \
data/metadata_files/ - includes individual-level metadata, and sample-level metadata \
data/tree_and_mutation_files/ - includes all saved objects relating to tree structures or mutation information from the WGS \
data/SV_and_CNA_data/ - includes summaries of structural variants and copy number alterrations from GRIDSS and ASCAT. Loss-of-Y information is derived from the mean coverage data. \
data/HDP/ - data files relating mutational signature extraction with HDP \
data/APOBEC_VCFs - vcf files containing only the likely APOBEC mutations from branches affected by APOBEC/ \
data/targeted_sequencing_data - raw and inferred data related to the targeted sequencing. The 'data_by_baitset' folder includes \
data/ABC_simulation_results - the posterior results from the models with different ABCs
