n<-commandArgs(trailingOnly=T)
print(n)
library(ape)
library(rsimpop)

##Import driver mutation parameters (from posterior of E. Mitchell et al, 2022)
params_path=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/Clonal_dynamics_of_HSCT/data/posterior_sample.txt","/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/posterior_sample.txt")
param_posterior<-read.delim(params_path,stringsAsFactors = F)

simpop_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_new/simulated_populations/"
system(paste0("mkdir -p ",simpop_dir))

age=75 #Set age to run to
HSC_pop_size=1e5
HSC_symmetric_division_rate=1/(2*365)
param_idx=sample(1:nrow(param_posterior),1)
number_drivers_per_year = param_posterior$number_drivers_per_year[param_idx]
gamma_shape = param_posterior$gamma_shape[param_idx]
gamma_rate = param_posterior$gamma_rate[param_idx]
fitness_threshold = 0.05

##Function to generate gamma distribution based fitness
genGammaFitness=function(fitness_threshold,shape,rate){
  function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
}
fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)

##Put parameters into a single list (to save later)
params=list(age=age,
            HSC_symmetric_division_rate=HSC_symmetric_division_rate,
            HSC_pop_size=HSC_pop_size,
            number_drivers_per_year=number_drivers_per_year,
            gamma_shape=gamma_shape,
            gamma_rate=gamma_rate,
            fitness_threshold=fitness_threshold)

pop_final=run_driver_process_sim(initial_division_rate=0.1,
                                 final_division_rate = HSC_symmetric_division_rate,
                                 target_pop_size = HSC_pop_size,
                                 nyears = age,
                                 fitnessGen=fitnessGammaFn,
                                 drivers_per_year = number_drivers_per_year)

saveRDS(object = pop_final,paste0(simpop_dir,"simulated_populations_",n,".Rds"))
print("Job completed, object saved")
