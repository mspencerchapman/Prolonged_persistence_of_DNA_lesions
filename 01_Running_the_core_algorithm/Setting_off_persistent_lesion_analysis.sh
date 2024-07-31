#Set file paths - the STUDY_DIR should be the cloned github directory
STUDY_DIR='/Users/ms56/R_work/Prolonged_persistence_of_DNA_lesions/'
GENOME_FILE='/Users/ms56/R_work/reference_files/genome.fa'

#Define the secondary file paths (these are as follows if the github directory has been cloned)
#Note that the 'input_data' directory needs to be downloaded from Mendeley Data as per the README file
LESION_SEG_INPUT_DIR="${STUDY_DIR}Data/input_data/"
SCRIPT_DIR=${STUDY_DIR}'02 Phasing_and_LOH_analysis/'
MAIN_ANALYSIS_SCRIPT_PATH=${STUDY_DIR}"01_Running_the_core_algorithm/Detect_persistent_lesions.R"
OUTPUT_DIR=${STUDY_DIR}output/


##Most of the datasets need to be run on a compute farm, but some smaller trees can be run locally
#Here is an example of running the algorithm for the 3 clonal haematopoiesis trees (which are fairly small)
#Note that these trees do include the ancestral tip

setopt shwordsplit
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
DATASET=MF
for EXP_ID in $SAMPLES; do
FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_noMixed.tree
Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done

#-----------------------------------------------------------------------------
#### SET OFF JOBS THE MAIN SCRIPT BY STUDY (as file name structures vary) ####
#-----------------------------------------------------------------------------
# Some jobs require a compute farm as they take a very long time/ require large amounts of memory
# However, for smaller trees could easily be run locally, removing the 'bsub' command
#i.e. just run as: Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s ${EXP_ID} -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -d -o ${OUTPUT_DIR}${DATASET}

####TRANSPLANT SAMPLES
#These trees do include the ancestral tip
DATASET=MSC_BMT
MEM=8000
QUEUE=normal
for EXP_NO in 11 13 21 24 25 28 31 38 40 41; do
	EXP_ID=Pair${EXP_NO}
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_m40_postMS_reduced_a_j_pval_post_mix.tree
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_m40_postMS_reduced_a_j_pval_post_mix
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s ${EXP_ID} -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -d -o ${OUTPUT_DIR}${DATASET} -a
done

####MY FOETAL SAMPLES
#These trees do not include the ancestral tip
DATASET=MSC_fetal
MEM=4000
QUEUE=normal
for EXP_ID in 8pcw 18pcw; do
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/Tree_${EXP_ID}.tree
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/Filtered_mut_set_annotated_${EXP_ID}
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -d -o ${OUTPUT_DIR}${DATASET}
done

####BRONCHIAL ORGANOID DATASET
#These trees do not include the ancestral tip
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}KY_samples.txt)
DATASET=KY
MEM=16000
QUEUE=basement
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/Filtered_muts_${EXP_ID}
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/${EXP_ID}_rmix_consense_tree_no_branch_lengths_1811.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET}
done

####NORMAL HAEMATOPOIESIS DATASET
#These trees do not include the ancestral tip
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}EM_samples.txt)
DATASET=EM
QUEUE=long
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_standard_rho01.tree
	MEM=24000
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET}
done

####CLONAL HAEMATOPOIESIS DATASET
#These trees do include the ancestral tip
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
DATASET=MF
MEM=4000
QUEUE=normal
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_noMixed.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done

####MYELOPROLIFERATIVE NEOPLASM DATASET
#These trees do include the ancestral tip

SAMPLES=$(ls input_data/NW|grep '.RDS'|cut -c1-6)
DATASET=NW
MEM=4000
QUEUE=normal
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/filtered_muts_${EXP_ID}
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done

#--------------------SETTING OFF PHASING/ SIMULATION JOBS--------------------
# Note that the phasing scripts will need considerable work to set up locally, but the scripts give info as to the logic
#MAV Phasing
MEM=16000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/PhaseMAV.log.%J -e ${STUDY_DIR}err_files/PhaseMAV.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J phaseMAV Rscript ${SCRIPT_DIR}Phasing_analysis_MAVs_2.R

#PVV Phasing
MEM=16000
QUEUE=normal
bsub -o ${STUDY_DIR}log_files/PhasePVV.log.%J -e ${STUDY_DIR}err_files/PhasePVV.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J phasePVV Rscript ${SCRIPT_DIR}Phasing_analysis_PVVs_2.R

#PVV Simulation - independent events
MEM=48000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/PVV_simulation.log.%J -e ${STUDY_DIR}err_files/PVV_simulation.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J PVVsim Rscript ${SCRIPT_DIR}PVV_simulation_2.R

#PVV Simulation - somatic reversion
MEM=16000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/PVV_revsim.log.%J -e ${STUDY_DIR}err_files/PVV_revsim.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J PVV_revsim Rscript ${SCRIPT_DIR}PVV_simulation_reversion_2.R

#MAV Simulation
MEM=16000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/MAV_simulation.log.%J -e ${STUDY_DIR}err_files/MAV_simulation.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J MAVsim Rscript ${SCRIPT_DIR}MAV_simulation_2.R

#LOH analysis - PVVs
MEM=16000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/LOH_analysis.log.%J -e ${STUDY_DIR}err_files/LOH_analysis.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J LOH_analysis Rscript ${SCRIPT_DIR}LOH_analysis_PVVs_2.R

#ASCAT LOH analysis - PVVs
MEM=32000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/LOH_ASCAT.log.%J -e ${STUDY_DIR}err_files/LOH_ASCAT.err.%J -q $QUEUE -R "select[mem>=$MEM] span[hosts=1] rusage[mem=$MEM]" -M${MEM} -n1 -J LOH_ASCAT Rscript ${SCRIPT_DIR}ASCAT_LOH_analysis_PVVs_2.R

#ASCAT LOH analysis - MAVs
MEM=16000
QUEUE=yesterday
bsub -o ${STUDY_DIR}log_files/LOH_MAV_ASCAT.log.%J -e ${STUDY_DIR}err_files/LOH_ASCAT.err.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J LOH_ASCAT Rscript ${SCRIPT_DIR}ASCAT_LOH_analysis_MAVs_2.R


#------------------------------
##### IQ TREE ANALYSIS ######
#------------------------------

OUTPUT_DIR=${STUDY_DIR}output_iqtree/

####EMILY's HSC SAMPLES iqtrees
SAMPLES=$(cat input_data/EM_samples.txt|tail -n6)
DATASET=EMiq
MEM=32000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET}
done

DATASET=EMiq
EXP_ID=KX004_5_01
MEM=128000
FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET}


####MSC_BMT HSC SAMPLES iqtrees
SAMPLES=$(cat input_data/MSC_BMT_samples.txt)
DATASET=MSC_BMTiq
MEM=32000

for EXP_NO in 11 13 21 24 25 28 31 38 40 41; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_Pair${EXP_NO}_m40_postMS_reduced_a_j_pval_post_mix
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_fasta_Pair${EXP_NO}.fa.allocated_muts.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/Pair${EXP_NO}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/Pair${EXP_NO}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s Pair${EXP_NO} -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done

####MARGA's SAMPLES iqtrees - these trees do include the ancestral tip
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
DATASET=MFiq
MEM=8000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done


#------------------------------
#### SCITE TREE ANALYSIS ######
#------------------------------

OUTPUT_DIR=${STUDY_DIR}output_scite/

####EMILY's HSC SAMPLES scite
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}EM_samples.txt)
DATASET=EMsc
MEM=64000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/scite_tree_allocated_muts_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET}
done

DATASET=EMsc
EXP_ID=KX004_5_01
MEM=128000
FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/scite_tree_allocated_muts_${EXP_ID}.tree
bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET}


####MSC_BMT HSC SAMPLES scite
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MSC_BMT_samples.txt)
DATASET=MSC_BMTsc
MEM=32000

for EXP_NO in 11 13 21 24 25 28 31 38 40 41; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_Pair${EXP_NO}_m40_postMS_reduced_a_j_pval_post_mix
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/scite_tree_allocated_muts_Pair${EXP_NO}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/Pair${EXP_NO}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/Pair${EXP_NO}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s Pair${EXP_NO} -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done

####MARGA's SAMPLES scite
#These trees do include the ancestral tip
DATASET=MFiq
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
MEM=8000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -g $GENOME_FILE -o ${OUTPUT_DIR}${DATASET} -a
done
