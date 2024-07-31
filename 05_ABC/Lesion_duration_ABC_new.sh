##ABC framework
cd /lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_new

#Create 100 aged HSC simulated populations (age = 75)
for i in {1..100}; do
	echo $i
	bsub -o $PWD/simpop_logs/log.%J -e $PWD/simpop_logs/log.%J -q long -R 'select[mem>=32000] span[hosts=1] rusage[mem=32000]' -M32000 -n1 -J simpop${i} Rscript generate_simpops.R $i
done

#Introduce lesions into each of these
#For each one, vary the lesion duration between 0.5 and 5 years
for i in `seq 0.5 0.1 5`; do 
	for k in {41..44}; do 
		echo "called with mean lesion duration=$i and population number=$k"
		bsub -o $PWD/introduce_lesion_logs/log.%J -e $PWD/introduce_lesion_logs/log.%J -q long -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' -M24000 -n1 -J les_sim Rscript ABC_simulation_new_INTRODUCE_LESIONS.R $k $i
	done
done

#Introduce lesions into each of these
#For each one, vary the lesion duration between 0.5 and 5 years
for i in `seq 0.5 0.1 0.6`; do 
	bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R 'select[mem>=48000] span[hosts=1] rusage[mem=48000]' -M48000 -n5 -J "lsim[1-40]%5" bash lesion_boost_wrapper.sh $i
done

#OR boost lesion numbers for individual
i=0.5
for k in 41 42 43 44; do 
	bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R 'select[mem>=100000] span[hosts=1] rusage[mem=100000]' -M100000 -n5 -J "lsim_ind" Rscript ABC_simulation_new_INTRODUCE_LESIONS_BOOST.R $k $i
done

MEM=48000
for i in `seq 0.5 0.1 0.8`; do
	for k in 41 42 43 44; do 
		bsub -o $PWD/introduce_lesion_logs/boosts/log.%J -e $PWD/introduce_lesion_logs/boosts/log.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n5 -J "lsim_ind" Rscript ABC_simulation_new_INTRODUCE_LESIONS_BOOST.R $k $i
	done
done


#Assess lesion capture for various lesion durations
for i in `seq 2.1 0.1 5`; do 
	for k in {1..3}; do 
		bsub -o $PWD/capture_logs/log.%J -e $PWD/capture_logs/log.%J -q long -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J les_cap Rscript ABC_simulation_new.R $i
	done
done

#Assess lesion capture for various lesion durations
for i in `seq 0.5 0.1 5`; do 
	bsub -o $PWD/capture_logs/log.%J -e $PWD/capture_logs/log.%J -q long -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n4 -J "lescap${i}[1-45]%10" Rscript ABC_simulation_new_parallel.R $i
done


for i in `seq 2.1 0.1 5`; do 
	bsub -o $PWD/capture_logs/GT.log.%J -e $PWD/capture_logs/GT.log.%J -q long -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n4 -J "GT_${i}" Rscript ABC_simulation_new_GROUND_TRUTH.R $i
done


MEM=100000
bsub -o $PWD/extractlog.%J -e $PWD/extractlog.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "extractinf" Rscript ABC_new_extract_info.R

#Run the posterior predictive checks
MEM=16000
bsub -o $PWD/capture_logs/PPClog.%J -e $PWD/capture_logs/PPClog.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n4 -J "PPC[1-227]%50" bash PPC_wrapper.sh
bsub -o $PWD/capture_logs/PPClog.%J -e $PWD/capture_logs/PPClog.%J -q long -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n4 -J "PPC[1-227]%50" bash PPC_wrapper2.sh