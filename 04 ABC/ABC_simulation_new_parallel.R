#### Assessing accuracy of ABC framework for inferring parameters

#The trailing argument is the filename for the lesion file
args<-commandArgs(trailingOnly = T)
print(args)
mean_LD_years=as.numeric(args[1])
test=F


# Import packages and define the custom functions ---------------------------------------------
library(ape)
library(rsimpop)
library(parallel)

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

get_ancestor_node=function(node,tree,degree=1){ #to get the 1st degree ancestor (i.e. the direct parent) use degree=1.  Use higher degrees to go back several generations.
  curr<-node
  for(i in 1:degree){
    curr=tree$edge[which(tree$edge[,2]==curr),1]
    if(curr==(1+length(tree$tip.label))) {stop(return(curr))}
  }
  return(curr)
}

#'Trim' the TLS base order, to remove bases that would not be detected in the phylogeny
# Removes (1) Reference bases at the beginning of the vector & (2) Matching bases (ref-ref OR alt-alt) at the end of the vector
assess_base_order_for_PVV<-function(bases) {
  
  if(!any(bases=="Alt")) {stop(return(list(detected=FALSE,bases_trimmed=NA)))}
  
  #Do the first trimming of any reference bases at the beginning
  first_alt_idx=min(which(bases=="Alt"))
  bases_trimmed<-bases[first_alt_idx:length(bases)]
  
  #Then trim the final base if it matches the preceding base - this will be not be detectable as two separate events
  if(length(bases_trimmed)>1) {
    while(length(bases_trimmed)>2&bases_trimmed[length(bases_trimmed)]==bases_trimmed[length(bases_trimmed)-1]){
      bases_trimmed<-head(bases_trimmed,n=-1)
    }
  }
  
  #Any detectable PVV will still have a vector of length ≥3 (most common is Alt-Ref-Alt, or Alt-Alt-Ref [equivalent] but longer vectors are possible)
  return(list(detected=ifelse(length(bases_trimmed)<=2,FALSE,TRUE),
              bases_trimmed=bases_trimmed))
}

estimate_gamma_params=function(value_vec,log_rate_range=c(-2,1),shape_range=c(1,5)) {
  require(tidyr)
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rate_vec = 10^(seq(log_rate_range[1],log_rate_range[2],by=0.1)) # rho will be bounded within 1e-6 and 0.89
  shape_vec=seq(shape_range[1],shape_range[2],0.05)
  params_grid=expand_grid(rate_vec,shape_vec)
  ll = sapply(1:nrow(params_grid), function(i) {shape=params_grid$shape_vec[i]; rate=params_grid$rate_vec[i];sum(dgamma(x=value_vec, shape=shape,rate=rate,log = T))})
  return(params_grid[which.max(ll),])
}


# Step 1: Define directories & import the 3 phylogenies ----------------------------------------
# Import 3 phylogenies - older phylogenies correspond to the 4 used to estimate the parameters of the PVVs
my_working_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/Lesion_Seg/ABC_new","/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_new")
setwd(my_working_dir)

simpop_dir="simulated_populations/"
lesion_sim_dir="lesion_simulation_files/"
captured_lesion_dir="captured_lesion_files/"

trees_to_use=sample(1:40,size = 4)
cat("Importing simpop population files",sep="\n")
simulated_populations<-lapply(trees_to_use,function(n) readRDS(paste0(simpop_dir,"simulated_populations_",n,".Rds")))  


# Step 2: Import the lesion info ------------------------------------------
lesion_info_list=lapply(trees_to_use,function(pop_number) readRDS(paste0(lesion_sim_dir,"lesion_simulation_",pop_number,"_",mean_LD_years,".Rds")))


# Step 3: Downsample trees & assess for PVV capture --------------------
data_ids=c("KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01")
tree_sizes=c(328,922,315,367) #These are the number of samples in KX003, KX004, KX007 and KX008

if(test) {
  nPVVs_in_data=c(2,2,2,2)
} else {
  nPVVs_in_data=c(33,80,9,22)
}

captured_PVV_id<-ids::random_id(1,bytes=8)
cat(paste("UID:",captured_PVV_id),sep="\n")

captured_lesion_set_list<-mclapply(1:4,mc.cores=4,FUN=function(j) {
  
  #Assign core variable
  nsamp=tree_sizes[j]
  nPVVs=nPVVs_in_data[j]
  pop_number=trees_to_use[j]
  simpop=simulated_populations[[j]]
  lesion_set=lesion_info_list[[j]]
  
  
  cat(paste("Assessing population",which(trees_to_use==pop_number),"with",nsamp,"samples in the tree."))
  lesion_info_PVVs=lesion_set$lesion_info_PVVs
  time_tree=get_elapsed_time_tree(simpop)
  tree_ss<-get_subsampled_tree2(simpop,N=nsamp)
  
  #Set the target number of PVVs as the number in the data + 20% - this gives a buffer incase the captured durations of some are >100 and are removed (as per filtering on the data)
  n_target_PVVs=round(nPVVs*1.2)
  
  #Now cycle through the theoretically detectable PVVs and see if any of the descendants from the lineage (or "pure subclades") are detectable in the downsampled tree
  cat(paste("There are",length(lesion_info_PVVs),"PVVs to assess to see if detected. However, detection will stop once/if the desired",n_target_PVVs,"PVVs are found."),sep="\n")
  lesion_info_PVVs2=vector(mode="list")
  n_detected_PVV_counter<-0
  reorder_vec=sample(1:length(lesion_info_PVVs),replace = F)
  for(k in 1:length(lesion_info_PVVs)) {
    if(k%%100==0) {cat(paste("Population",j,":",k, "lesions assessed.",n_detected_PVV_counter,"/",n_target_PVVs,"PVVs captured so far."),sep="\n")}
    idx<-reorder_vec[k] #so that don't always cycle through the lesions from the beginning
    list<-lesion_info_PVVs[[idx]]
    pure_subclades<-list$pure_subclades
    subclades_detected=sapply(pure_subclades,function(node) {
      subclade_detected=any(tree_ss$tip.label%in%getTips(time_tree,node))
      return(subclade_detected)
    })
    
    if(any(subclades_detected)) {
      list$captured_pure_subclades<-sapply(pure_subclades[subclades_detected], function(node) {
        captured_pure_subclade_samples=tree_ss$tip.label[tree_ss$tip.label%in%getTips(time_tree,node)]
        assigned_node<-find_latest_acquisition_node(tree = tree_ss,pos_samples = captured_pure_subclade_samples)
        return(assigned_node)
      })
      
      list$captured_pure_subclade_bases<-list$pure_subclade_bases[subclades_detected]
      names(list$captured_pure_subclade_bases)<-list$captured_pure_subclades
      
      if(length(list$captured_pure_subclades)<=2) {
        list$detected_as_PVV<-F
      } else {
        list$detected_as_PVV<-assess_base_order_for_PVV(list$captured_pure_subclade_bases)$detected
        
        #If the PVV is detected, increase the counter by 1
        if(list$detected_as_PVV) {n_detected_PVV_counter<-n_detected_PVV_counter+1; cat(paste("Population",j,":",n_detected_PVV_counter,"of",n_target_PVVs,"detected"),sep="\n")}
      }
      
    } else {
      list$captured_pure_subclades<-"None"
      list$detected_as_PVV<-F
    }
    lesion_info_PVVs2[[k]]<-list
    if(n_detected_PVV_counter==n_target_PVVs) {cat(paste("Population",j,":","Reached the target number of captured PVVs. Stopping after assessment of",k,"PVVs."),sep="\n");break}
  }
  
  if(n_detected_PVV_counter<n_target_PVVs) {cat(paste("Population",j,":","All lesions assessed and target number of PVVs not reached."),sep="\n")}
  
  #Keep only the captured lesions
  captured_lesion_info_PVVs<-lesion_info_PVVs2[sapply(lesion_info_PVVs2,function(list) list$detected)]
  
  cat(paste("Population",j,":","Total recorded number of captured PVVs:",length(captured_lesion_info_PVVs)),sep="\n")
  
  captured_lesion_set=list(pop_number=lesion_set$pop_number,
                           lesion_simulation_id=lesion_set$lesion_simulation_id,
                           mean_lesion_duration_years=lesion_set$mean_lesion_duration_years,
                           mean_lesion_duration_days=lesion_set$mean_lesion_duration_days,
                           gamma_shape=lesion_set$gamma_shape,
                           alt_base_prob=lesion_set$alt_base_prob,
                           n_lesions=lesion_set$n_lesions,
                           captured_PVV_id=captured_PVV_id,
                           nsample=tree_ss$ntips-1,
                           tree_ss=tree_ss,
                           captured_lesion_info_PVVs=captured_lesion_info_PVVs)
  return(captured_lesion_set)
})
names(captured_lesion_set_list)<-data_ids

# Final steps: save the output --------------------------------------------
cat("Completed all lesion capture analyses. Saving results.",sep="\n")
saveRDS(captured_lesion_set_list,file = paste0(captured_lesion_dir,"captured_lesions_",mean_LD_years,"_",captured_PVV_id,".Rds"))
cat("Script completed. Saving results.",sep="\n")



