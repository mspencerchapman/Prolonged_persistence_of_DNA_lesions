#Running the lesion segregation script
STUDY_DIR='/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/'
LESION_SEG_INPUT_DIR='/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/'
SCRIPT_DIR=$STUDY_DIR
MAIN_ANALYSIS_SCRIPT_PATH='/lustre/scratch126/casm/team154pc/ms56/my_programs/Lesion_segregation_mutation_summaries_MNVs_4.R'
OUTPUT_DIR=${STUDY_DIR}output2/

####MY TRANSPLANT SAMPLES

DATASET=MSC_BMT
MEM=8000

for EXP_NO in 11 13 28 31 38 40 41; do
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_Pair${EXP_NO}_m40_postMS_reduced_a_j_pval_post_mix.tree
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_Pair${EXP_NO}_m40_postMS_reduced_a_j_pval_post_mix
	bsub -o ${STUDY_DIR}log_files/${DATASET}/Pair${EXP_NO}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/Pair${EXP_NO}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s Pair${EXP_NO} -t $TREE_PATH -f $FILTERED_MUTS_PATH -d -o ${OUTPUT_DIR}${DATASET} -a
done

####MY FOETAL SAMPLES

DATASET=MSC_fetal
MEM=16000
for EXP_ID in 8pcw 18pcw; do
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/Tree_${EXP_ID}.tree
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/Filtered_mut_set_annotated_${EXP_ID}
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -d -o ${OUTPUT_DIR}${DATASET}
done

####DAVID KENT SAMPLES

DATASET=DK
MEM=8000
for EXP_ID in BCL002 BCL009 CGD1 CGD2; do
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_m40_postMS_reduced_pval_post_mix.tree
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_m40_postMS_reduced_pval_post_mix
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -d -o ${OUTPUT_DIR}${DATASET}
done

####KENICHI's BRONCHIAL SAMPLES

SAMPLES=$(cat input_data/KY_samples.txt)
DATASET=KY
MEM=16000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/Filtered_muts_${EXP_ID}
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/${EXP_ID}_rmix_consense_tree_no_branch_lengths_1811.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}
done

####EMILY's HSC SAMPLES

SAMPLES=$(cat input_data/EM_samples.txt)
DATASET=EM

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_standard_rho01.tree
	MEM=24000
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}
done

####PHIL'S COLON SAMPLES

SAMPLES=$(cat input_data/PR_samples.txt)
DATASET=PR
MEM=16000
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/${EXP_ID}/Filtered_muts_${EXP_ID}
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/${EXP_ID}/snp_tree_with_branch_length_polytomised.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}
done

####MARGA's SAMPLES
#These trees do include the ancestral tip

SAMPLES=$(cat input_data/MF_samples.txt)
DATASET=MF
MEM=4000
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/filtered_muts_${EXP_ID}_noMixed
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_noMixed.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET} -a
done

####NICK's SAMPLES
#These trees do include the ancestral tip

SAMPLES=$(ls input_data/NW|grep '.RDS'|cut -c1-6)
DATASET=NW
MEM=4000
for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/filtered_muts_${EXP_ID}
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET} -a
done

#--------------------SETTING OFF PHASING/ SIMULATION JOBS--------------------

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
##### ANALYSE ALL MUTATIONS FOR OVERDISPERSION ######
#------------------------------

ANALYSIS_SCRIPT_COMPLETE='/lustre/scratch126/casm/team154pc/ms56/my_programs/Lesion_segregation_mutation_summaries_MNVs_all.R'
OUTPUT_COMPLETE_DIR=${STUDY_DIR}output_complete/

####EMILY's HSC SAMPLES
DATASET=EM
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}${DATASET}_samples.txt)
MEM=64000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/tree_${EXP_ID}_standard_rho01.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript ${ANALYSIS_SCRIPT_COMPLETE} -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_COMPLETE_DIR}${DATASET}
done

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
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}
done

DATASET=EMiq
EXP_ID=KX004_5_01
MEM=128000
FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}


####MSC_BMT HSC SAMPLES iqtrees
SAMPLES=$(cat input_data/MSC_BMT_samples.txt)
DATASET=MSC_BMTiq
MEM=32000

for EXP_NO in 11 13 21 24 25 28 31 38 40 41; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_Pair${EXP_NO}_m40_postMS_reduced_a_j_pval_post_mix
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_fasta_Pair${EXP_NO}.fa.allocated_muts.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/Pair${EXP_NO}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/Pair${EXP_NO}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s Pair${EXP_NO} -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET} -a
done

####MARGA's SAMPLES iqtrees - these trees do include the ancestral tip
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
DATASET=MFiq
MEM=8000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET} -a
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
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}
done

DATASET=EMsc
EXP_ID=KX004_5_01
MEM=128000
FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_standard_rho01
TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/scite_tree_allocated_muts_${EXP_ID}.tree
bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET}


####MSC_BMT HSC SAMPLES scite
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MSC_BMT_samples.txt)
DATASET=MSC_BMTsc
MEM=32000

for EXP_NO in 11 13 21 24 25 28 31 38 40 41; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_Pair${EXP_NO}_m40_postMS_reduced_a_j_pval_post_mix
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/scite_tree_allocated_muts_Pair${EXP_NO}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/Pair${EXP_NO}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/Pair${EXP_NO}.err.%J -q basement -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s Pair${EXP_NO} -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET} -a
done

####MARGA's SAMPLES scite
#These trees do include the ancestral tip
DATASET=MFiq
SAMPLES=$(cat ${LESION_SEG_INPUT_DIR}MF_samples.txt)
MEM=8000

for EXP_ID in $SAMPLES; do
	FILTERED_MUTS_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/annotated_mut_set_${EXP_ID}_pval
	TREE_PATH=${LESION_SEG_INPUT_DIR}${DATASET}/iqtree_tree_allocated_muts_${EXP_ID}.tree
	bsub -o ${STUDY_DIR}log_files/${DATASET}/${EXP_ID}.log.%J -e ${STUDY_DIR}err_files/${DATASET}/${EXP_ID}.err.%J -q normal -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J lesion_seg Rscript $MAIN_ANALYSIS_SCRIPT_PATH -w ${STUDY_DIR} -s $EXP_ID -t $TREE_PATH -f $FILTERED_MUTS_PATH -o ${OUTPUT_DIR}${DATASET} -a
done
