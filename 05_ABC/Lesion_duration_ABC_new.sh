ABC_DIR='/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_new/'

##ABC framework
cd $ABC_DIR

#----------------------------------------------------
###1. GENERATE THE POPULATIONS
#----------------------------------------------------

#Create 40 aged HSC simulated populations for the ABC (age = 75)
#AND 4 additional populations for assessing the performance of the ABC
MEM=32000
for i in {1..44}; do
	echo $i
	bsub -o $PWD/simpop_logs/log.%J -e $PWD/simpop_logs/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J simpop${i} Rscript generate_simpops.R $i
done

#----------------------------------------------------
###2. INTRODUCE SIMULATED LESIONS INTO EACH OF THESE
#----------------------------------------------------

#Introduce lesions into each of these
#For each one, vary the lesion duration between 0.5 and 5 years
MEM=24000
for i in `seq 0.5 0.1 5`; do 
	for k in {1..44}; do 
		echo "called with mean lesion duration=$i and population number=$k"
		bsub -o $PWD/introduce_lesion_logs/log.%J -e $PWD/introduce_lesion_logs/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J les_sim Rscript ABC_simulation_new_INTRODUCE_LESIONS.R $k $i
	done
done

#Now 'BOOST' the number of lesions to get at least 60,000 detectable PVVs
#This works out how many PVVs were generated with the initial lesions, and therefore infers how many additional lesions are needed to get to 60,000
#This is done with a simple wrapper script to allow each 'mean lesion duration' to be run as an array
MEM=24000
for i in `seq 0.7 0.1 5`; do 
	bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n5 -J "lsim[1-44]%5" bash lesion_boost_wrapper.sh $i
done

#The lower lesion durations need additional memory as so many lesions need to be introduced to get 60,000 PVVs
MEM=48000
i=0.6
bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n5 -J "lsim[1-44]%5" bash lesion_boost_wrapper.sh $i

MEM=100000
i=0.5
bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n5 -J "lsim[1-44]%5" bash lesion_boost_wrapper.sh $i


##If have an issue with some individual populations, can use these loops to run the script without the wrapper
MEM=48000
for i in `seq 0.5 0.1 0.8`; do
	for k in 41 42 43 44; do 
		bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n5 -J "lsim_ind" Rscript ABC_simulation_new_INTRODUCE_LESIONS_BOOST.R $k $i
	done
done

#----------------------------------------------------
###3. USE THE PREVIOUS INFO TO RUN THE ACTUAL DATA SIMULATION
#----------------------------------------------------

#Assess lesion capture for various lesion durations
MEM=16000
for i in `seq 0.5 0.1 5`; do 
	bsub -o $PWD/capture_logs/log.%J -e $PWD/capture_logs/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n4 -J "lescap${i}[1-55]%10" Rscript ABC_simulation_new_parallel.R $i
done


#----------------------------------------------------
###4. ASSESS PERFORMANCE OF THE ABC ON NEW SIMULATED DATA (treating populations 41 - 44 as the data)
#----------------------------------------------------

MEM=16000
for i in `seq 2.1 0.1 5`; do 
	bsub -o $PWD/capture_logs/GT.log.%J -e $PWD/capture_logs/GT.log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n4 -J "GT_${i}" Rscript ABC_simulation_new_GROUND_TRUTH.R $i
done


#----------------------------------------------------
###5. EXTRACT THE KEY DATA FROM THE 'LESION CAPTURE' SIMULATIONS FOR RUNNING THE ABC
#----------------------------------------------------
MEM=100000
bsub -o $PWD/extractlog.%J -e $PWD/extractlog.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "extractinf" Rscript ABC_new_extract_info.R

#----------------------------------------------------
###6. RUN THE POSTERIOR PREDICTIVE CHECKS
#----------------------------------------------------
#Run the posterior predictive checks
#Run simulations for each parameter of the 227 posterior values
MEM=16000
bsub -o $PWD/capture_logs/PPClog.%J -e $PWD/capture_logs/PPClog.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n4 -J "PPC[1-227]%50" bash PPC_wrapper.sh