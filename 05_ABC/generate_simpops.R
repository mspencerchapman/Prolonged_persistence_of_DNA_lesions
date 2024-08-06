n<-commandArgs(trailingOnly=T)
print(n)

#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ape","remotes")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("rsimpop", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/rsimpop")
  library("rsimpop",character.only=T,quietly = T, warn.conflicts = F)
}
options(stringsAsFactors = FALSE)

#========================================#
# Read parameters and ensure directories are set ####
#========================================#

#Set directories
root_dir="./" #Or else set to the repository path
simpop_dir=paste0(root_dir,"05_ABC/simpops/")
system(paste0("mkdir -p ",simpop_dir))

##Import driver mutation parameters (from posterior of E. Mitchell et al, 2022)
params_path=paste0(root_dir,"Data/reference_files/ABC_parameters/driver_parameter_posterior_sample.txt")
param_posterior<-read.delim(params_path,stringsAsFactors = F)

#========================================#
# Set all the parameters for the simulation ####
#========================================#

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

##Put parameters into a single list (to save later if desired)
params=list(age=age,
            HSC_symmetric_division_rate=HSC_symmetric_division_rate,
            HSC_pop_size=HSC_pop_size,
            number_drivers_per_year=number_drivers_per_year,
            gamma_shape=gamma_shape,
            gamma_rate=gamma_rate,
            fitness_threshold=fitness_threshold)

#========================================#
# Run the simulation ####
#========================================#

pop_final=run_driver_process_sim(initial_division_rate=0.1,
                                 final_division_rate = HSC_symmetric_division_rate,
                                 target_pop_size = HSC_pop_size,
                                 nyears = age,
                                 fitnessGen=fitnessGammaFn,
                                 drivers_per_year = number_drivers_per_year)

saveRDS(object = pop_final,paste0(simpop_dir,"simulated_populations_",n,".Rds"))
print("Job completed, object saved")
