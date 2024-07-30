library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(plotrix)
library(phangorn)
library(RColorBrewer)
library(Quartet)
options(stringsAsFactors = F)

#Define function required later in script for comparing phylogenies
comparePhylo_and_plot=function(tree1,tree2,names){
  plot_comp_tree=function(tree,comp,title,col="red",lwd=1,tree_pos){
    tree_name=deparse(substitute(tree))
    shared_clades=comp[,tree_pos]
    edge_width=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),lwd,2*lwd))
    edge_col=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),"black",col))
    plot(tree,show.tip.label=F,direction="downwards",edge.color=edge_col,edge.width=edge_width,main=title)
  }
  comp<-compare_nodes(tree1,tree2)
  par(mfrow=c(1,2))
  plot_comp_tree(tree1,comp=comp,title=names[1],tree_pos = 1)
  plot_comp_tree(tree2,comp=comp,title=names[2],tree_pos = 2)
}

get_RF_dist=function(tree1,tree2){
  comp<-comparePhylo(tree1,tree2)
  RF_dist=1-(length(comp$NODES[,1])/length(unique(tree$edge[,1])))
  return(RF_dist)
}

#Function to import the different trees, either on the farm, or via the mounted directories
import_tree_file=function(SampleID,data_set,type) {
  if(!type%in%c("mpboot","iqtree","scite")) {stop(print("Type must be either mpboot, iqtree, or scite"))}
  root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
  
  if(type=="mpboot") {
    files=list.files(paste0(root_dir,"/lesion_segregation/input_data/",data_set,"/"),pattern = ".tree",full.names = T)
    relevant_tree_path=grep(paste0("_",SampleID),x=files,value = T)
  } else {
    if(data_set=="MSC_BMT") {
      HSCT_dir=paste0(root_dir,"/Zur_HSCT/final_versions")
      relevant_pattern=ifelse(type=="iqtree",".fa.allocated_muts.tree","scite_tree_allocated_muts_")
      relevant_tree_files=list.files(path=HSCT_dir,pattern=relevant_pattern,full.names = T)
      relevant_tree_path=grep(SampleID,relevant_tree_files,value=T)
      
    } else if(data_set=="MF") {
      relevant_tree_path=ifelse(type=="iqtree",paste0(root_dir,"/Marga_benchmarking/",SampleID,"/iqtree_tree_allocated_muts_",SampleID,".tree"),paste0(root_dir,"/Marga_benchmarking/",SampleID,"/scite_tree_allocated_muts_",SampleID,".tree"))
      tree=ape::read.tree(relevant_tree_path)
    } else if(data_set=="EM") {
      shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
      relevant_tree_path=ifelse(type=="iqtree",paste0(root_dir,"/Emily_benchmarking/",shortID,"/vaf_filtered/iqtree_tree_allocated_muts_",SampleID,".tree"),paste0(root_dir,"/Emily_benchmarking/",shortID,"/vaf_filtered/scite_tree_allocated_muts_",SampleID,".tree"))
    }
  }
  
  #Return the relevant tree if it exists
  if(file.exists(relevant_tree_path)) {
    print(paste("The",type,"tree for",SampleID,"is being read in."))
    tree=read.tree(relevant_tree_path)
    return(tree)
  } else {
    print(paste("The",type,"tree for",SampleID,"is not found at the expected location."))
    return(NULL)
  }
}

import_muts_file=function(SampleID,data_set) {
  root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
  
  data_set_muts_files=list.files(path=paste0(root_dir,"/lesion_segregation/input_data/",data_set),pattern="annotated",full.names = T)
  filepath=grep(SampleID,data_set_muts_files,value=T)
  load(filepath)
  return(filtered_muts)
}

compare_nodes=function(x, y) 
{
  tree1 <- deparse(substitute(x))
  tree2 <- deparse(substitute(y))
  n1 <- Ntip(x)
  n2 <- Ntip(y)
  
  key1 <- makeNodeLabel(x, "md5sum")$node.label
  key2 <- makeNodeLabel(y, "md5sum")$node.label
  mk12 <- match(key1, key2)
  mk21 <- match(key2, key1)
  if (any(tmp <- is.na(mk12))) {
    nk <- sum(tmp)
  }
  if (any(tmp <- is.na(mk21))) {
    nk <- sum(tmp)
  }
  nodes1 <- which(!is.na(mk12))
  nodes2 <- mk12[!is.na(mk12)]
  
  NODES <- data.frame(nodes1 + n1,nodes2 + n2)
  names(NODES) <- c(tree1, tree2)
  return(NODES)
}

my_theme=theme_classic(base_family="Helvetica")+theme(text=element_text(size=7,family="Helvetica"),
                                                              axis.text=element_text(size=7,family="Helvetica"),
                                                              strip.text = element_text(size=7,family="Helvetica"),
                                                              legend.text = element_text(size=7,family="Helvetica"),
                                                              axis.title=element_text(size=7,family="Helvetica"),
                                                  axis.title.y=element_text(size=7,family="Helvetica"),
                                                  axis.title.x=element_text(size=7,family="Helvetica"))


#Source functions needed for the script
my_working_directory<-getwd()
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)


#Settings for this individual
for(SampleID in data_set_samples) {
  cat(SampleID,sep="\n")
  exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1]
  filtering="vaf_filtered"
  filtering_settings="standard_rho01"
  
  #Define the file paths for the data files
  my_working_directory=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/tree_bootstraps")
  tree_file_path=paste0("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/EM/tree_",SampleID,"_standard_rho01.tree")
  file_annot=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/annotated_mut_set_",SampleID,"_",filtering_settings)
  bootstrapped_tree_stats_file_path=paste0("Bootstrapped_tree_stats_",SampleID)
  
  #Source functions needed for the script
  R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
  treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
  sapply(R_function_files[-2],source)
  if(!dir.exists(my_working_directory)) {next}
  setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
  
  #Read count bootstraps
  bootstrapped_trees_dir=paste0(my_working_directory,"/final_trees/")
  all_bootstrapped_trees =paste0(my_working_directory,"/all_trees")
  #bootstrap_tree_files=list.files(path=bootstrapped_trees_dir,full.names = T,pattern=".tree")
  
  if(!file.exists(all_bootstrapped_trees)) {
    system("cat final_trees/Tree_*>all_trees")
  }
  
  all_boots=read.tree(all_bootstrapped_trees)
  all_boots<-di2multi(all_boots)
  
  #Load up the main tree
  tree=di2multi(read.tree(tree_file_path))
  ROOT=1+length(tree$tip.label)
  tree_node_numbers=unique(tree$edge[,1]) #store list of the tree node numbers (not including the tips)
  
  if(!all(all_boots[[1]]$tip.label%in%tree$tip.label)) {
    cat("Dropping tips from the bootstrapped trees that aren't in the data.")
    all_boots<-lapply(all_boots,function(boot_tree) drop.tip(boot_tree,tree$tip.label))
  }
  
  
  #Do comparisons for both sets of bootstraps - Quartet status
  boots_comparisons_quartets=Quartet::QuartetStatus(all_boots,cf=tree)
  Quartet_divergence_plot<-as.data.frame(SimilarityMetrics(boots_comparisons_quartets))%>%
    dplyr::select(QuartetDivergence,SymmetricDifference)%>%
    gather(key="metric",value="value")%>%
    filter(metric=="QuartetDivergence")%>%
    ggplot(aes(x=value,fill=metric))+
    geom_density(linewidth=0.2)+
    facet_grid(cols=vars(metric),scales="free_y")+
    scale_x_continuous(limits=c(0.9,1))+
    my_theme + theme(legend.position = "none") +
    labs(x="Similarity score",y="Density")
  
  as.data.frame(SimilarityMetrics(boots_comparisons_quartets))%>%
    dplyr::select(QuartetDivergence,SymmetricDifference)%>%
    gather(key="metric",value="value")%>%
    filter(metric=="QuartetDivergence")%>%
    summarise(median=median(value))
  
  boots_comparisons_splits=Quartet::SplitStatus(all_boots,cf=tree)
  RF_dist_plot<-data.frame(RF_dist=RawSymmetricDifference(boots_comparisons_splits) / boots_comparisons_splits[, 'N'])%>%
    ggplot(aes(x=RF_dist))+
    geom_density(size=0.2)+
    scale_x_continuous(limits=c(0,0.2))+
    my_theme + theme(legend.position = "none") +
    labs(x="Robinson-Foulds distance",y="Density")
  
  data.frame(RF_dist=RawSymmetricDifference(boots_comparisons_splits) / boots_comparisons_splits[, 'N'])%>%
    
    
    RF_sim_plot<-data.frame(SimilarityMetrics(boots_comparisons_splits))%>%
    gather(key="metric",value="value")%>%
    dplyr::filter(metric=="SymmetricDifference")%>%
    ggplot(aes(x=value,fill=metric))+
    geom_density(size=0.2)+
    facet_grid(cols = vars(metric))+
    scale_x_continuous(limits=c(0.9,1))+
    my_theme + theme(legend.position = "none") +
    labs(x="Robinson-Foulds distance",y="Density")
  
  data.frame(SimilarityMetrics(boots_comparisons_splits))%>%
    gather(key="metric",value="value")%>%
    dplyr::filter(metric=="SymmetricDifference")%>%
    summarise(median=median(value))
  
  comb_plots=arrangeGrob(Quartet_divergence_plot,RF_sim_plot,nrow=2)
  plot(comb_plots)
  ggsave(comb_plots,filename = paste0("Bootstrap_score_comparisons_",exp_ID,"_",filtering,".pdf"),device="pdf",width = 3,height=4)
  
  #Create list of comparison stats to compare each of the bootstrap trees
  boot_stats_list=list(nmuts=numeric(),
                       mean_mut_burden=numeric(),
                       nnodes_boot=numeric(),
                       shared_clades=list(),
                       new_clades=list())
  
  #Fill up the stats list from the data
  for(i in 1:length(all_boots)) {
    print(i)
    tree_boot=all_boots[[i]]
    comp=compare_nodes(di2multi(tree),tree_boot) #Run the tree comparison function (hacked from the ape "comparePhylo" function)
    
    #Extract parameters of interest
    tree_boot_node_numbers=unique(tree_boot$edge[,1])
    tree_boot_lost_clades=tree_boot_node_numbers[which(!tree_boot_node_numbers%in%comp$tree_boot)]
    tree_boot_new_clades=lapply(tree_boot_lost_clades,function(node) {getTips(tree_boot,node)})
    
    #Store the output
    boot_stats_list$nmuts[i]<-sum(tree_boot$edge.length)
    boot_stats_list$mean_mut_burden[i]<-mean(get_mut_burden(tree_boot))
    boot_stats_list$nnodes_boot[i]<-length(tree_boot_node_numbers)
    boot_stats_list$shared_clades[[i]]<-as.numeric(comp[,1])
    boot_stats_list$tree_boot_new_clades[[i]]<-tree_boot_new_clades
  }
  
  
  # save(boot_stats_list,file = bootstrapped_tree_stats_file_path)
  # 
  # load(bootstrapped_tree_stats_file_path)
  
  retained_node_read_bootstraps=sapply(tree_node_numbers,function(node) {
    sum(unlist(lapply(boot_stats_list$shared_clades,function(x) node%in%as.numeric(x))))/length(boot_stats_list$nmuts)
  })
  
  RF_dist_read_bootstraps=sapply(boot_stats_list$shared_clades,function(x) 1-(length(x)/length(tree_node_numbers)))
  
  #Save the data as to the proportion of bootstraps in which the node is retained
  node_confidence_df<-data.frame(node=tree_node_numbers,retained=retained_node_read_bootstraps)%>%
    filter(node!=ROOT)%>%
    mutate(Sample_ID=exp_ID,.before=1)
  
  write.table(node_confidence_df,file=paste0("bootstrap_retained_",SampleID,".tsv"),quote=F,sep="\t",row.names = F)
}
