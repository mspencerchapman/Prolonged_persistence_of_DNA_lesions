##Script to introduce 100,000 lesions into the different aged simpop populations
#The population and the mean lesion duration are set outside the script

args<-commandArgs(trailingOnly = T)
print(args)
pop_number=as.numeric(args[1])
mean_lesion_duration_years=as.numeric(args[2])
print(pop_number)
print(mean_lesion_duration_years)

simpop_dir="simulated_populations/"
lesion_sim_dir="lesion_simulation_files/"
output_file=paste0(lesion_sim_dir,"lesion_simulation_",pop_number,"_",mean_lesion_duration_years,".Rds")

#### Assessing accuracy of ABC framework for inferring parameters

##------Outline for approach--------
# (1) Simulate phylogenies of ageing blood
# (2) Introduce persistent lesions of set durations (according to a gamma distribution, with varying mean)
# (3) Visualize PVVs resulting from these and their durations
# (4) Work back using the same ABC framework to see if we recover the correct lesion durations


library(ape)
library(rsimpop)
my_working_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/Lesion_Seg/ABC_new","/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_new")
setwd(my_working_dir)

##Import driver mutation parameters (from posterior of E. Mitchell et al, 2022)
params_path=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/Clonal_dynamics_of_HSCT/data/posterior_sample.txt","/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/posterior_sample.txt")
param_posterior<-read.delim(params_path,stringsAsFactors = F)


# Define custom functions -------------------------------------------------

#Do updated 'get_subsampled_tree' function that maintains original tip labels (needed to extract bulk clonal fractions)
get_subsampled_tree2=function (tree, N, tips = tree$edge[c(which(tree$state == 0 & 
                                                                   tree$edge[, 2] <= length(tree$tip.label)), sample(which(tree$state != 
                                                                                                                             0 & tree$edge[, 2] <= length(tree$tip.label)), N)), 2]) 
{
  N = length(tips)
  tip.labels<-tree$tip.label[sort(tips)]
  
  #Create the tree with the tips kept in the same order
  tree_same_order<-keep.tip(tree,tip=tip.labels)
  tmp = rsimpop:::C_subsample_pop(tree, sort(tips))
  tmp$tip.label = sprintf("s%d", 1:N)
  class(tmp) = c("simpop", "phylo")
  
  #Now reorder the tips
  node_trans<-all.equal(tree_same_order,tmp,use.tip.label=F,index.return = T)
  node_trans<-node_trans[node_trans[,2]<=length(tree_same_order$tip.label),]
  new_tips<-tip.labels[node_trans[,2]]
  tmp$tip.label<-new_tips
  
  #tmp$tip.label = tip.labels
  
  checkValidPhylo(tmp)
  tmp$is_combined = tree$is_combined
  tmp
}

find_latest_acquisition_node=function(tree,pos_samples){
  #Get list of ancestral nodes for all samples
  ancestral_nodes_list=lapply(pos_samples,function(Sample) {
    get_ancestral_nodes(node = which(tree$tip.label==Sample),edge=tree$edge)
  })
  #Find nodes that are ancestral to all the samples
  common_nodes=Reduce(intersect,ancestral_nodes_list)
  #Which of these is the most recent (i.e. has the maximum node height)
  nodeheights<-sapply(common_nodes,function(node) nodeheight(tree = tree,node = node))
  MRCA_node<-common_nodes[which.max(nodeheights)]
  return(MRCA_node)
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

get_ancestral_nodes=function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}

assess_base_order_for_PVV<-function(bases) {
  first_alt_idx=min(which(bases=="Alt"))
  bases_trimmed<-bases[first_alt_idx:length(bases)]
  
  while(length(bases_trimmed)>2&bases_trimmed[length(bases_trimmed)]==bases_trimmed[length(bases_trimmed)-1]){
    bases_trimmed<-head(bases_trimmed,n=-1)
  }
  return(list(detected=ifelse(length(bases_trimmed)<=2,FALSE,TRUE),
              bases_trimmed=bases_trimmed))
}


##---STEP 1: Read in the specified simpop------------------
simpop<-readRDS(paste0(simpop_dir,"simulated_populations_",pop_number,".Rds"))

##---STEP 1a: Read in the previously generated lesions------------------
previous_lesions<-readRDS(output_file)
n_existing_PVV<-length(previous_lesions$lesion_info_PVVs)
target_PVV=60000
boost_required=target_PVV/n_existing_PVV
cat(paste("There are",n_existing_PVV,"PVVs so far.", target_PVV,"PVVs are desired."),sep="\n")

if(boost_required>1) {
  ##---STEP 2: Introduce the persistent lesions---------------
  #Here we have complete trees of 100,000 cells
  #There are time stamps (in days) of when each node in the tree existed
  
  time_tree<-get_elapsed_time_tree(simpop)
  
  #Decide the lesion duration distribution and base pairing probabilities
  id<-ids::random_id(1,bytes = 8)
  #mean_lesion_duration_years=runif(1,min=0.1,max=3)
  mean_lesion_duration_days=mean_lesion_duration_years*365
  gamma_shape=1
  alt_base_prob=0.513 #taken from the data
  
  
  #How many lesions will we introduce into the phylogeny
  # - need to introduce more if the lesions are short, otherwise won't capture enough to get a distribution
  if(mean_lesion_duration_years<=2) {
    n_lesions=5e6
  } else {
    n_lesions=2e6
  }
  
  n_lesions=round(n_lesions*(boost_required-1))
  cat(paste("Introducing an additional",n_lesions,"lesions to try and reach desired number of PVVs."),sep="\n")
  
  #Introduce the lesions one by one
  cat(paste("Running with 5 cores"),sep="\n")
  lesion_info=parallel::mclapply(1:n_lesions,mc.cores=5,function(k) {
    if(k%%1000==0) {cat(k,sep="\n")}
    #1. select a branch where lesion is generated; equally likely to be generated throughout life
    selected_branch<-sample(time_tree$edge[,2],size=1,replace=T,prob = time_tree$edge.length)
    selected_branch_length=time_tree$edge.length[time_tree$edge[,2]==selected_branch]
    position_on_branch<-runif(1)
    
    #2. How long will the lesion last?
    lesion_duration=rgamma(n=1,shape=gamma_shape,rate=(1/mean_lesion_duration_days))
    
    #3. Now 'chart the course' of the lesion through the phylogeny
    # Need to record (1) which nodes are passed, (2) the nodes defining the clades with 'fixed' alt/ ref alleles
    # The lesion will 'stop' one the lesion duration is reached
    time_to_next_node=(1-position_on_branch)*selected_branch_length
    nodes_traversed=vector(); edge_length_vec=vector(); pure_subclades=vector()
    next_node=selected_branch
    duration_passed=time_to_next_node
    while(duration_passed<lesion_duration&!next_node%in%1:length(time_tree$tip.label)) {
      nodes_traversed=c(nodes_traversed,next_node)
      daughters<-time_tree$edge[time_tree$edge[,1]==next_node,2]
      next_node=sample(daughters,size=1,replace = T)
      pure_subclades=c(pure_subclades,daughters[!daughters==next_node])
      edge_length=time_tree$edge.length[time_tree$edge[,2]==next_node]
      duration_passed=duration_passed+edge_length
      edge_length_vec=c(edge_length_vec,edge_length)
    }
    
    #The node path and the 'pure subclades' vectors need to include the final node in the path (even though the lesion never actually reaches that node)
    node_path=c(nodes_traversed,next_node)
    pure_subclades=c(pure_subclades,next_node)
    
    #Return the list of this info
    res<-list(selected_branch=selected_branch,
              selected_branch_length=selected_branch_length,
              position_on_branch=position_on_branch,
              lesion_duration=lesion_duration,
              nodes_traversed=nodes_traversed,
              node_path=node_path,
              edge_length_vec=edge_length_vec,
              pure_subclades=pure_subclades) #The assumption is that the lesion may be 'lost' by translesion resynthesis
    return(res)
  })
  
  #Now add the base pairing info on top of this - i.e. whether the ref or alt base is paired during cell divisions
  lesion_info2=lapply(lesion_info,function(list) {
    if(length(list$pure_subclades)<3) {
      list$PVV_outcome<-"None"
    } else {
      bases=sample(c("Alt","Ref"),replace=T,size=length(list$pure_subclades),prob = c(alt_base_prob,1-alt_base_prob))
      names(bases)<-list$pure_subclades
      if(sum(bases=="Alt")>1){
        #Trim any initial reference bases which will be invisible on the phylogeny
        first_alt_idx=min(which(bases=="Alt"))
        bases_trimmed<-bases[first_alt_idx:length(bases)]
        
        #Accordingly trim the number of nodes that are 'proven' to have been traversed
        proven_traversed_nodes=list$nodes_traversed[first_alt_idx:length(list$nodes_traversed)]
        
        #print(bases)
        #Trim any final bases that are identical to previous bases - also invisible on phylogeny
        while(length(bases_trimmed)>2&bases_trimmed[length(bases_trimmed)]==bases_trimmed[length(bases_trimmed)-1]){
          bases_trimmed<-head(bases_trimmed,n=-1)
          proven_traversed_nodes<-head(proven_traversed_nodes,n=-1)
        }
        #print(bases)
        if(length(bases_trimmed)>2) {
          list$PVV_outcome<-"PVV"
          list$pure_subclade_bases<-bases
          list$proven_traversed_nodes<-proven_traversed_nodes
          list$MMLD<-sum(sapply(proven_traversed_nodes[-1],function(node) time_tree$edge.length[time_tree$edge[,2]==node]))
        } else {
          list$PVV_outcome<-"None"
          list$pure_subclade_bases<-bases
        }
      } else {
        list$PVV_outcome<-"None"
        list$pure_subclade_bases<-bases
      }
    }
    return(list)
  })
  
  ##---STEP 3: Include only lesions that result in a potentially detectable PVV & save---------------
  lesion_info_PVVs=lesion_info2[which(sapply(lesion_info2,function(list) list$PVV_outcome)=="PVV")]
  cat(paste(length(lesion_info_PVVs),"new PVVs detectable."),sep="\n")
  
  lesion_set=list(pop_number=pop_number,
                  lesion_simulation_id=id,
                  mean_lesion_duration_years=mean_lesion_duration_years,
                  mean_lesion_duration_days=mean_lesion_duration_days,
                  gamma_shape=gamma_shape,
                  alt_base_prob=alt_base_prob,
                  n_lesions=n_lesions,
                  lesion_info_PVVs=lesion_info_PVVs)
  
  #Read in the previous lesions sets & update with the newly generated info
  lesion_set$lesion_simulation_id=c(previous_lesions$lesion_simulation_id,lesion_set$lesion_simulation_id)
  lesion_set$n_lesions=sum(previous_lesions$n_lesions,lesion_set$n_lesions)
  lesion_set$lesion_info_PVVs=c(previous_lesions$lesion_info_PVVs,lesion_set$lesion_info_PVVs)
  
  saveRDS(lesion_set,file=output_file)
  
} else {
  cat("There are already adequate PVVs. Ending script.")
}

