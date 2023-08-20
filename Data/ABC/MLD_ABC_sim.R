#!/software/R-3.6.1/bin/Rscript
#THIS SCRIPT ATTEMPTS TO GET A MEASURE OF HOW MUCH THE TREE STRUCTURE IS LIKELY TO MAKE IT POSSIBLE TO PICK UP PVVs/ MAVs
#It works out (for a particular lesion duration) how likely that lesion is to be present at a subsequent node if it has
#equal probability of being picked up all the way along a branch

library(GenomicRanges)
library(IRanges)
library("Rsamtools")
library("MASS")
library(stringr)
library(dplyr)
library(tidyr)
library(ape)
options(stringsAsFactors = FALSE)

my_working_directory<-getwd()

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(rev(R_function_files[-2]),source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Set data file paths
data_sets=c("NW","MF","EM","KY","PR","MSC_BMT","MSC_fetal")
out_list=lapply(data_sets,function(data_set) {
  files=list.files(paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/output2/",data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
  mut_tables_list=lapply(files,read.delim)
  mut_tables_df=Reduce(rbind,mut_tables_list)
  mut_tables_df$data_set=data_set
  return(mut_tables_df)
})
mutations=Reduce(rbind,out_list)
mutations$Ref[mutations$Ref=="TRUE"]<-"T"
mutations$Alt1[mutations$Alt1=="TRUE"]<-"T"
mutations$Alt2[mutations$Alt2=="TRUE"]<-"T"

genome_file="/lustre/scratch119/casm/team154pc/ms56/genome.fa"

#This first script tests what proportion of persistent DNA lesions would create a potentially detectable PVV
# It goes through a number of steps:
# 1. It "squashes" the early embryonic tree, as the hypothesis is that the causative lesion is primarily in the quiescent HSC cell state/ bone marrow environment (not during embryonic development)
# 2. It models a "persistent DNA lesion" being introduced randomly along branches of the tree (the branch is chosen randomly in proportion to the branch length, and the position along the branch is chosen randomly [runif(0,1)])
# 3. The persistent lesion is also assigned a "lesion duration" - random draw from a gamma distribution with a mean selected from the prior
# 4. The lesion is tested to see if it lasts for ≥2 nodes, with the path down the tree chosen at random (i.e. daughter branches of each node selected at random)
# 5. If the lesion does cross ≥2 nodes, the resulting bases (Alt or Ref) in each subclade are chosen at random according to the "alt_base_probability". The maximum chance of a PVV is if this is set to 0.66.
# 6.The resulting bases are interrogated to see if a PVV would be detectable & what is the implied "minimum lesion duration" (which may be much less than the ACTUAL lesion duration)

#This data is summarised in terms of:
#1. The proportion of introduced persistent lesions that result in detectable PVVs
#2. The mean MLD of captured PVVs
#3. Characteristics of the distribution of MLDs of PVVs, when approximated as a gamma distribution with rate & shape parameters

alt_base_prob=0.54 #taken from the data
distribution="gamma"

res_list=lapply(seq(2,100,by=2),function(mean_lesion_duration) {
  print(paste("Running simulation with mean lesion duration of",mean_lesion_duration))
  dataset="EM"
  Sample_IDs<-c("KX004_5_01","KX007_2_01","KX008_2_01")
  data_set_out=lapply(Sample_IDs,function(sample) {
    print(paste("Running simulation for",sample))
    print(paste("Starting analysis for sample",sample))
    #Find the relevant project and tree file
    sample_info=get_file_paths_and_project(dataset,Sample_ID=sample)
    tree=read.tree(sample_info$tree_file_path)
    
    if(!sample%in%c("CB001_3_01","CB002_2_01","8pcw","18pcw")){
      tree<-squash_tree(tree,cut_off=60,from_root = T)
    }
    #n=round(0.15*sum(tree$edge.length))
    n=5e6
    
    #Get vectors of length n for each of the 3 variables (1. lesion duration, 2. which node, 3. where in branch)
    if(distribution=="poisson"){
      lesion_duration=rpois(n=n,lambda=mean_lesion_duration)
    } else if(distribution=="gamma") {
      lesion_duration=rgamma(n=n,shape=1,rate=(1/mean_lesion_duration))
    }
    nodes=base::sample(x = tree$edge[,2],size = n,replace = T,prob = tree$edge.length)
    where_in_branch=runif(n=n)
    
    #Use these vectors to work out if lesion results in a detectable PVV
    out=mapply(function(node,where,lesion_duration) {
      #If the lesion occurs on a private branch, it will never result in a detectable PVV
      if(node%in%1:length(tree$tip.label)|all(tree$edge[tree$edge[,1]==node,2] %in% 1:length(tree$tip.label))){
        return("FAIL")
      } else {
        #If lesion is on a shared branch, work out how far along the branch it has occurred and
        #how many nodes it will "cross" i.e. for how many captured cell divisions was the lesion present
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
    
    return(list(MLDs=captured_MLDs,
                 df=data.frame(Sample_ID=sample,
                               nsim=n,
                               mean_lesion_duration=mean_lesion_duration,
                               capture_rate=capture_rate,
                               tot_score=capture_rate*sum(tree$edge.length),
                               mean_MLD=mean(captured_MLDs),
                               gamma_shape=gamma_params_mld$shape_vec,
                               gamma_rate=gamma_params_mld$rate_vec)))
  })
  return(data_set_out)
})

saveRDS(object = res_list,file = "MLD_ABC_output.RDS")

