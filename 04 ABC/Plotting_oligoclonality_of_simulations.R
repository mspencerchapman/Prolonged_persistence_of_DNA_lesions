#### OLIGOCLONALITY PLOT

# Import packages and define the custom functions ---------------------------------------------
library(ape)
library(rsimpop)
library(parallel)
library(dplyr)
library(stringr)


#========================================#
# Set the ggplot2 themes for plotting ####
#========================================#

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 7),
                axis.title = element_text(size=8),
                axis.line = element_line(linewidth = 0.4),
                axis.ticks = element_line(linewidth = 0.3),
                legend.text = element_text(size=6),
                legend.title = element_text(size=8),
                strip.text = element_text(size=7),
                strip.background = element_rect(fill="lightgray",linewidth = 0.4),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))



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


get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

# Define directories & import the simulated phylogenies ----------------------------------------
my_working_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/Lesion_Seg/ABC_new","/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_benchmarking")
setwd(my_working_dir)

simpop_dir="simulated_populations/"
lesion_sim_dir="lesion_simulation_files/"
captured_lesion_dir="captured_lesion_files/"

# Import the 40 simulated populations used in the ABC, and downsample to similar size to the data ----
n_sim=40
simpop_files=list.files(simpop_dir,full.names = T)
all_simpops<-lapply(simpop_files[1:n_sim],function(file) readRDS(file))

nsamp=400
all_sim_trees<-lapply(all_simpops,function(simpop) {
  time_tree<-get_elapsed_time_tree(simpop,mutrateperdivision = 1,backgroundrate = 16/365)
  tree_ss<-keep.tip(time_tree,tip = sample(time_tree$tip.label,size=nsamp))
  return(tree_ss)
})

# Get oligoclonality data on the simulated trees ----------------------------------------
sim.expanded.clades<-dplyr::bind_rows(Map(tree=all_sim_trees,exp_ID=paste("sim",1:n_sim,sep="_"),function(tree,exp_ID){
  cat(exp_ID)
  exp_nodes<-get_expanded_clade_nodes(tree,height_cut_off = 100,min_clonal_fraction=0.01) #Clones originating after 100 mutaitons of molecular time
  if(nrow(exp_nodes)>0){
    exp_nodes$exp_ID<-exp_ID
    return(exp_nodes)
  } else {
    return(NULL)
  }
}))

# Read in all the data trees ----------------------------------------

library(dplyr)
library(ape)
source(paste0(root_dir,"/my_programs/Prolonged_persistence_functions.R"))
mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
lesion_seg_input_dir=paste0(root_dir,"/lesion_segregation/input_data")
PVV_blood_table<-mutations%>%
  filter(Type=="PVV" & cat%in%c("Adult_HSPC","Chemo_HSPC","Foetal_HSPC") &!data_set=="NW" & Class=="PASS")

data_sets=c("MSC_BMT","MF","EM","MSC_fetal")
all_trees=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  all_data_set_trees<-lapply(data_set_samples,function(sample_ID) {
    tree<-import_tree_file(SampleID=sample_ID,data_set = data_set,type = "mpboot")
    return(tree)
  })
  names(all_data_set_trees)<-data_set_samples
  return(all_data_set_trees)
})%>%unlist(recursive = F)

# Get oligoclonality on the old trees used for the ABC ----------------------------------------
old_individuals=c("KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01")
data.expanded.clades<-dplyr::bind_rows(Map(tree=all_trees[old_individuals],exp_ID=old_individuals,function(tree,exp_ID){
  exp_nodes<-get_expanded_clade_nodes(tree,height_cut_off = 100,min_clonal_fraction=0.01)
  exp_nodes$exp_ID<-exp_ID
  return(exp_nodes)
}))%>%mutate(type="Data")

individual_order<-dplyr::bind_rows(sim.expanded.clades,data.expanded.clades)%>%
  group_by(exp_ID)%>%
  dplyr::summarise(total_clone=sum(clonal_fraction))%>%
  arrange(total_clone)%>%
  pull(exp_ID)

expanded.clades.plot<-sim.expanded.clades%>%
  mutate(type="Simulation")%>%
  bind_rows(data.expanded.clades)%>%
  arrange(desc(clonal_fraction))%>%
  ggplot(aes(x=factor(exp_ID,levels = individual_order),y=clonal_fraction,fill=type))+
  geom_bar(stat="identity",position="stack",col="black",alpha=0.5,linewidth=0.3)+
  scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(labels=ifelse(grepl("sim",individual_order),"",stringr::str_split(individual_order,pattern="_",simplify=T)[,1]))+ #only show labels for the data
  labs(x="",y="Fractions of\n clonal expansions",fill="")+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust=+0.5))+
  labs(x="")


ggsave(filename=paste0(plots_dir,"Oligoclonality_of_simulations_vs_data.pdf"),expanded.clades.plot,width=4,height=2.5)
