# HDP Flow II: run single chains
# Tim Coorens, Feb 2020

#Would normally run this script in parallel ~20 - 30 times using the bash script
# Can then combine the single chains using the script: hdp_combine_results.R

#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("devtools")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

options(stringsAsFactors = F)

#========================================#
# Set paths for running ####
#========================================#

root_dir<-"~/R_work/Prolonged_persistence_of_DNA_lesions/"
HDP_folder=paste0(root_dir,"/Data/HDP/")

#HDP_folder="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/mutsig_extract/"

for(data_set in c("KY","EM")) {
  
  #========================================#
  # Set paths for running ####
  #========================================#
  setwd(paste0(HDP_folder,data_set))
  
  for(n in 1:20) {
    cat(n,sep="\n")
    lower_threshold=50 #Only samples (or tree branches) with at least 50 mutations will be included
    #n=as.numeric(commandArgs(T)[1]) #This is the 'chain number' for saving the output
    
    #Read in the input files ---
    mutations=read.table("trinuc_mut_mat.txt")
    key_table=read.table("key_table.txt")
    
    #If requiring a minimum number of mutations:
    sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
    mutations=mutations[!rownames(mutations)%in%sample_remove,]
    key_table=key_table[!key_table$Sample%in%sample_remove,]
    
    #Hierarchy is set per patient, can change if wanted ----
    freq=nrow(mutations)
    
    hdp_mut <- hdp_init(ppindex = c(0, rep(1,length(freq)),rep(2:(length(freq)+1), times=freq)), # index of parental node
                        cpindex = c(1, rep(2,length(freq)),rep(3:(length(freq)+2), times=freq)), # index of the CP to use
                        hh = rep(1, 96), # prior is uniform over 96 categories
                        alphaa = rep(1,length(freq)+2), # shape hyperparameters for 2 CPs
                        alphab = rep(1,length(freq)+2))  # rate hyperparameters for 2 CPs
    
    hdp_mut <- hdp_setdata(hdp_mut, 
                           dpindex = (length(freq)+2):numdp(hdp_mut), # index of nodes to add data to
                           mutations)
    
    hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10,seed=n*300)
    
    chain=hdp_posterior(hdp_activated,
                        burnin=20000,
                        n=100,
                        seed=n*1000,
                        space=200,
                        cpiter=3)
    
    #Save output chain ----
    saveRDS(chain,paste0("hdp_chain_",n,".Rdata"))
  }
}


