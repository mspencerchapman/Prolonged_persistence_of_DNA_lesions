#THIS SCRIPT ATTEMPTS TO GET A MEASURE OF HOW MUCH THE TREE STRUCTURE IS LIKELY TO MAKE IT POSSIBLE TO PICK UP PVVs/ MAVs
#It works out (for a particular lesion duration) how likely that lesion is to be present at a subsequent node if it has
#equal probability of being picked up all the way along a branch

#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr")
bioconductor_packages=c("GenomicRanges","IRanges","Rsamtools")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}
options(stringsAsFactors = FALSE)

#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir="~/R_work/Prolonged_persistence_of_DNA_lesions/"
data_dir=paste0(root_dir,"Data/")
plots_dir=paste0(root_dir,"plots/")
output_dir=paste0(root_dir,"/output/")
source(paste0(root_dir,"/Data/Prolonged_persistence_functions.R"))
ref_table=read.csv(paste0(root_dir,"/Data/metadata/Individual_ref.csv"))
genome_file=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/reference_files/genome.fa","/nfs/cancer_ref02/human/GRCh37d5/genome.fa") ##Set to the location of GRCh37 genome file

Phasing_MAV_file_path=paste0(data_dir,"Phasing_results_MAVs_all")
Phasing_PVV_file_path=paste0(data_dir,"Phasing_results_PVVs_all")
ASCAT_PVV_file_path=paste0(data_dir,"ASCAT_LOH_analysis_PVVs_all")
ASCAT_MAV_file_path=paste0(data_dir,"ASCAT_LOH_analysis_MAVs_all")
SN_phasing_file_path=paste0(data_dir,"Phasing_results_MAVs_SN.csv")

lesion_seg_input_dir=paste0(data_dir,"input_data/")
lesion_seg_output_dir=output_dir

## Import the mutations & summary_df files ----
mutations_filt_file=paste0(data_dir,"mutations_filtered.tsv")
summary_df_file=paste0(data_dir,"summary_df.tsv")

mutations<-readr::read_delim(mutations_filt_file,col_types = "cccccccciicccciiiiiiccccccnccccccciccccci")
summary_table_df<-read.delim(summary_df_file)
sample_ref=read.csv(paste0(data_dir,"metadata/Individual_ref.csv"))

#========================================#
# Define CUSTOM FUNCTION for script ####
#========================================#

squash_tree=function(tree,cut_off=50,from_root=F) {
  if(from_root){
    idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
    new_edge_lengths=nodeHeights(tree)[idxs_to_squash,2]-cut_off #work-out the edge lengths that these should be such that they finish at the cut-off
    tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
    
    tree$edge.length[nodeHeights(tree)[,2]<=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
    return(tree)
  } else {
    tree$edge.length[nodeHeights(tree)[,1]>=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
    idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
    new_edge_lengths=cut_off - nodeHeights(tree)[idxs_to_squash,1] #work-out the edge lengths that these should be such that they finish at the cut-off
    tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
    return(tree)
  }
}

#========================================#
# TEST THE CAPTURE RATE of SIMULATED PVVs ####
#========================================#

#This first script tests what proportion of persistent DNA lesions would create a potentially detectable PVV
# It goes through a number of steps:
# 1. It "squashes" the early embryonic tree, as the hypothesis is that the causative lesion is primarily in the quiescent HSC cell state/ bone marrow environment
# 2. It models a "persistent DNA lesion" being introduced randomly along branches of the tree (the branch is chosen randomly in proportion to the branch length, and the position along the branch is chosen randomly [runif(0,1)])
# 3. The persistent lesion is also assigned a "lesion duration" - random draw from a gamma distribution with a mean selected from the prior
# 4. The lesion is tested to see if it lasts for ≥2 nodes, with the path down the tree chosen at random (i.e. daughter branches of each node selected at random)
# 5. If the lesion does cross ≥2 nodes, the resulting bases (Alt or Ref) in each subclade are chosen at random according to the "alt_base_probability". The maximum chance of a PVV is if this is set to 0.66.
# 6.The resulting bases are interrogated to see if a PVV would be detectable & what is the implied "minimum lesion duration" (which may be much less than the ACTUAL lesion duration)

#This data is summarised in terms of:
#1. The proportion of introduced persistent lesions that result in detectable PVVs
#2. The mean MLD of captured PVVs
#3. Characteristics of the distribution of MLDs of PVVs, when approximated as a gamma distribution with rate & shape parameters

mean_lesion_duration=20
alt_base_prob=0.66
distribution="poisson"

Capture_rate_list=lapply(data_sets, function(dataset) {
  Sample_IDs<-mutations%>%filter(data_set==dataset)%>%pull(Sample_ID)%>%unique()
  data_set_out=lapply(Sample_IDs,function(sample) {
    print(paste("Starting analysis for sample",sample))
    sample_info=get_file_paths_and_project(dataset,Sample_ID=sample,input_data_dir = lesion_seg_input_dir)
    tree=read.tree(sample_info$tree_file_path)
    
    if(!sample%in%c("CB001_3_01","CB002_2_01","8pcw","18pcw")){
      tree<-squash_tree(tree,cut_off=60,from_root = T)
    }
    #n=round(0.15*sum(tree$edge.length))
    n=1e6
    
    #Get vectors of length n for each variable
    if(distribution=="poisson"){
      lesion_duration=rpois(n=n,lambda=mean_lesion_duration)
    } else if(distribution=="gamma") {
      lesion_duration=rgamma(n=n,shape=1,rate=(1/mean_lesion_duration))
    }
    
    nodes=base::sample(x = tree$edge[,2],size = n,replace = T,prob = tree$edge.length)
    where_in_branch=runif(n=n)
    
    #Use these vectors to work out if lesion is present at time of dichotomy
    out=mapply(function(node,where,lesion_duration) {
      if(node%in%1:length(tree$tip.label)|all(tree$edge[tree$edge[,1]==node,2] %in% 1:length(tree$tip.label))){
        return("FAIL")
      } else {
        edge_length=tree$edge.length[tree$edge[,2]==node]
        time_to_next_node=(1-where)*edge_length
        duration_passed=time_to_next_node
        nodes_traversed=vector()
        n_nodes_traversed=0
        edge_length_vec=vector()
        next_node=node
        while(duration_passed<lesion_duration&!next_node%in%1:length(tree$tip.label)) {
          nodes_traversed=c(nodes_traversed,next_node)
          n_nodes_traversed=n_nodes_traversed+1
          next_node=sample(tree$edge[tree$edge[,1]==next_node,2],size=1)
          edge_length=tree$edge.length[tree$edge[,2]==next_node]
          duration_passed=duration_passed+edge_length
          edge_length_vec=c(edge_length_vec,edge_length)
        }
        
        #Assign paired bases to each of the "pure subclades"
        bases=sample(c("Alt","Ref"),replace=T,size=length(edge_length_vec)+1,prob = c(alt_base_prob,1-alt_base_prob))
        names(bases)<-c(nodes_traversed,nodes_traversed[length(nodes_traversed)])
        #print(bases)
        
        #Need at least 2 Alt bases to be detectable
        if(sum(bases=="Alt")>1){
          #Trim any initial reference bases which will be invisible on the phylogeny
          bases<-bases[min(which(bases=="Alt")):length(bases)]
          #print(bases)
          #Trim any final bases that are identical to previous bases - also invisible on phylogeny
          while(length(bases)>2&bases[length(bases)]==bases[length(bases)-1]){
            bases<-head(bases,n=-1)
          }
          #print(bases)
          if(length(bases)>2) {
            proven_traversed_nodes=head(tail(names(bases),n=-1),n=-1)
            MLD=sum(sapply(proven_traversed_nodes,function(node) tree$edge.length[tree$edge[,2]==node]))
            return(MLD)
          } else {
            stop(return("FAIL"))
          }
        } else {
          stop(return("FAIL"))
        }
      }
    }, node=nodes,where=where_in_branch,lesion_duration=lesion_duration)
    
    capture_rate=1-(table(out)["FAIL"]/length(out))
    captured_MLDs=as.numeric(out[out!="FAIL"])
    gamma_params_mld=estimate_gamma_params(value_vec = captured_MLDs)

    return(data.frame(Sample_ID=sample,
                      nsim=n,
                      mean_lesion_duration=mean_lesion_duration,
                      capture_rate=capture_rate,
                      tot_score=capture_rate*sum(tree$edge.length),
                      mean_MLD=mean(captured_MLDs),
                      gamma_shape=gamma_params_mld$shape_vec,
                      
                      gamma_rate=gamma_params_mld$rate_vec))
  })
  data_set_out=dplyr::bind_rows(data_set_out)
  return(data_set_out)
})

Capture_rate_df<-dplyr::bind_rows(Capture_rate_list)

write.table(Capture_rate_df,file=paste0(root_dir,"Data/simulation_results/capture_rate_sim_poisson.tsv"),row.names=F,quote=F,sep="\t")