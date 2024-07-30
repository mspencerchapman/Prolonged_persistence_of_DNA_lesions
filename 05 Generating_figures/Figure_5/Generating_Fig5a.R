#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("abc","ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr")
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
# Set the ggplot2 themes for plotting ####
#========================================#

my_theme<-theme_classic()+
  theme(text = element_text(family="Helvetica"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size=8),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8),
        strip.text = element_text(size=8),
        strip.background = element_rect(fill=NA,linewidth = 0.4),
        legend.spacing = unit(1,"mm"),
        legend.key.size= unit(5,"mm"))

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

#Set the sample IDs used for the ABC
ABC_Sample_IDs=c("KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01")

#========================================#
# Define custom functions for script ####
#========================================#
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

# Function 'Trim' the TLS base order, to remove bases that would not be detected in the phylogeny
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

#========================================#
# COMPARE OLIGOCLONALITY OF SIMULATIONS to DATA ####
#========================================#

##NB. The simulated population files are quite large & are therefore not available on github
# However, you can easily generate similar simulated populations using the 'generate_simpops.R' script

ABC_dir=paste0(root_dir,"Data/ABC_results/")
lesion_sim_dir=paste0(ABC_dir,"lesion_simulation_files/")
simpop_dir=paste0(ABC_dir,"simulated_populations")
captured_lesion_dir=paste0(ABC_dir,"captured_lesion_files/")
simpop_files=list.files(simpop_dir,full.names = T)

# my_working_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/Lesion_Seg/ABC_new","/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/ABC_benchmarking")
# setwd(my_working_dir)

## Import the 40 simulated populations used in the ABC, and downsample to similar size to the data ----
n_sim=40
simpop_files=list.files(simpop_dir,full.names = T)
all_simpops<-lapply(simpop_files[1:n_sim],function(file) readRDS(file))

nsamp=400
all_sim_trees<-lapply(all_simpops,function(simpop) {
  time_tree<-get_elapsed_time_tree(simpop,mutrateperdivision = 1,backgroundrate = 16/365)
  tree_ss<-keep.tip(time_tree,tip = sample(time_tree$tip.label,size=nsamp))
  return(tree_ss)
})

## Get oligoclonality data on the simulated trees ----------------------------------------
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


data_sets=c("EM")
all_trees=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  all_data_set_trees<-lapply(data_set_samples,function(sample_ID) {
    tree<-read.tree(get_file_paths_and_project(dataset=data_set,Sample_ID=sample_ID,input_data_dir = lesion_seg_input_dir)$tree_file_path)
    return(tree)
  })
  names(all_data_set_trees)<-data_set_samples
  return(all_data_set_trees)
})%>%unlist(recursive = F)

## Get oligoclonality on the old trees used for the ABC -----
data.expanded.clades<-dplyr::bind_rows(Map(tree=all_trees[ABC_Sample_IDs],exp_ID=ABC_Sample_IDs,function(tree,exp_ID){
  exp_nodes<-get_expanded_clade_nodes(tree,height_cut_off = 100,min_clonal_fraction=0.01)
  exp_nodes$exp_ID<-exp_ID
  return(exp_nodes)
}))%>%mutate(type="Data")

## Generate Extended Data Fig. 11a ----
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
  guides(fill=guide_legend(position="inside"))+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust=+0.5),legend.position.inside = c(0.2,0.75))+
  labs(x="")

ggsave(filename=paste0(plots_dir,"ExtDatFig11a.pdf"),expanded.clades.plot,width=3,height=2.5)

#========================================#
# Get simulation results in format for extract summary statistics ####
#========================================#

##--------Test the total branch lengths of the different simpops: are they comparable in terms of total 'lineage time'--------
root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
plots_dir=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/Prolonged_persistence_of_DNA_lesions/Rebuttal_plots/","/lustre/scratch126/casm/team154pc/ms56/")
ABC_dir=paste0(root_dir,"Data/ABC_results/")
lesion_sim_dir=paste0(ABC_dir,"lesion_simulation_files/")
simpop_dir=paste0(ABC_dir,"simulated_populations")
simpop_files=list.files(simpop_dir,full.names = T)

lineage_time_file=paste0(ABC_dir,"lineage_time.tsv")
if(file.exists(lineage_time_file)) {
  lt_df<-read.delim(lineage_time_file)
} else {
  lt_vec<-sapply(simpop_files,function(file) {
    simpop<-readRDS(file); tree<-get_elapsed_time_tree(simpop)
    sum(tree$edge.length)
  })
  names(lt_vec)<-list.files(simpop_dir,full.names = F)
  lt_df<-data.frame(file=names(lt_vec),lineage_time=lt_vec)%>%
    tibble::remove_rownames()%>%
    mutate(pop_number=gsub("simulated_populations_","",gsub(".Rds","",file)),.before=1)%>%
    mutate(pop_number=as.integer(pop_number))%>%
    arrange(lt_df)
  write.table(lt_df,file = lineage_time_file,quote = F,sep="\t",row.names=F)
}

#Extract information on the total number of PVVs for each run - needed to infer the number of lesions per cell
total_PVV_info_file=paste0(ABC_dir,"total_PVV_info.tsv")
if(file.exists(total_PVV_info_file)) {
  total_PVV_df<-read.delim(total_PVV_info_file)
} else {
  #If you have rerun the simulations, then can extract the information from them here
  lesion_sim_files=list.files(lesion_sim_dir,full.names = T)
  total_PVV_df<-lapply(lesion_sim_files,function(file) {
    les_sim<-readRDS(file)
    data.frame(pop_number=les_sim$pop_number,mean_lesion_duration_years=les_sim$mean_lesion_duration_years,total_PVVs_in_set=length(les_sim$lesion_info_PVVs))
  })%>%dplyr::bind_rows()
  write.table(total_PVV_df,file = total_PVV_info_file,quote = F,sep="\t",row.names=F)
}


## Import the captured lesion files and extract the lesion & parameter data from each -----
setwd(ABC_dir)

reextract=F
if(reextract) {
  ##Can run this code to extract information from the raw captured lesion files if you have re-run the ABC
  cat("Re-extracting info from all captured lesion files",sep="\n")
  files=list.files(path=paste0(ABC_dir,"captured_lesion_files"),full.names = T)
  uids=sapply(files,function(filename) {x<-gsub(".Rds","",filename); gsub("^.*_", "", x)})
  
  cat(paste("Extracting information from",length(files),"files"),sep="\n")
  sim_res<-lapply(files,function(file) {
    cat(file,sep="\n")
    x<-readRDS(file)
    temp=lapply(x, function (individual) {
      simpop<-readRDS(paste0(simpop_dir,"/simulated_populations_",individual$pop_number,".Rds"))
      tree_ss_raw<-get_subsampled_tree2(simpop,tips=which(simpop$tip.label%in%individual$tree_ss$tip.label))
      tree_mut<-get_elapsed_time_tree(tree_ss_raw,mutrateperdivision = 1,backgroundrate = 16/365)
      
      lesions<-individual$captured_lesion_info_PVVs
      this_sim_df<-lapply(1:length(lesions),function(i) {
        this_lesion_info<-lesions[[i]]
        LN<-get_ancestor_node(this_lesion_info$captured_pure_subclades[1],tree=tree_mut)
        LRN<-get_ancestor_node(tail(this_lesion_info$captured_pure_subclades,n=1),tree=tree_mut)
        LN_timing=nodeheight(tree_mut,LN)
        LRN_timing=nodeheight(tree_mut,LRN)
        MMLD=LRN_timing-LN_timing
        df<-data.frame(lesion_duration=this_lesion_info$lesion_duration,
                       selected_branch=this_lesion_info$selected_branch,
                       LN=LN,
                       LRN=LRN,
                       pure_subclade_nodes=paste(this_lesion_info$pure_subclades,collapse="-"),
                       pure_subclade_bases=paste(this_lesion_info$pure_subclade_bases,collapse = "-"),
                       LN_timing=LN_timing,
                       LRN_timing=LRN_timing,
                       MLD_days=this_lesion_info$MMLD,
                       MLD_MT=MMLD)
      })%>%dplyr::bind_rows()
      
      if(is.null(individual$n_PVVs_assessed)) {
        individual_info=list(lesion_simulation_id=paste(individual$lesion_simulation_id,collapse=","),
                             captured_PVV_id=individual$captured_PVV_id,
                             pop_number=individual$pop_number,
                             total_lesions=individual$n_lesions,
                             lineage_time=lt_vec[paste0("simulated_populations_",individual$pop_number,".Rds")],
                             mean_lesion_duration_years=individual$mean_lesion_duration_years)
      } else {
        
        individual$total_PVVs_in_set<-total_PVV_df%>%filter(pop_number==individual$pop_number & mean_lesion_duration_years==individual$mean_lesion_duration_years)%>%pull(total_PVVs_in_set)
        individual_info=list(lesion_simulation_id=paste(individual$lesion_simulation_id,collapse=","),
                             captured_PVV_id=individual$captured_PVV_id,
                             pop_number=individual$pop_number,
                             n_PVVs_assessed=individual$n_PVVs_assessed,
                             total_PVVs=individual$total_PVVs_in_set,
                             prop_PVVs_assessed=individual$n_PVVs_assessed/individual$total_PVVs_in_set,
                             total_lesions=individual$n_lesions,
                             implied_total_lesions=individual$n_lesions*individual$n_PVVs_assessed/individual$total_PVVs_in_set,
                             lineage_time=(lt_df%>%filter(pop_number==individual$pop_number)%>%pull(lineage_time)),
                             lesion_density=individual$n_lesions*individual$n_PVVs_assessed/individual$total_PVVs_in_set/(lt_df%>%filter(pop_number==individual$pop_number)%>%pull(lineage_time)),
                             mean_lesion_duration_years=individual$mean_lesion_duration_years)
      }
      
      return(list(df=this_sim_df,tree=tree_mut,individual_info=individual_info))
      
    })
    
    duration_list<-lapply(temp,function(list) {list$df$MLD_MT})
    individual_info_list<-lapply(temp,function(list) {list$individual_info})
    
    return(list(info=list(captured_PVV_id=x[[1]]$captured_PVV_id,mean_lesion_duration_days=x[[1]]$mean_lesion_duration_days),duration_list=duration_list,individual_info=individual_info_list))
  })
  
  saveRDS(sim_res,file = paste0(ABC_dir,"extracted_lesion_durations.Rds"))
  
} else {
  cat("Importing previously saved info from captured lesion files",sep="\n")
  sim_res<-readRDS(paste0(ABC_dir,"extracted_lesion_durations.Rds"))
}

#========================================#
# Get summary stats from simulations ####
#========================================#

##Get the mean/ dispersion coeffecients from these mutation sets ----
all_res<-lapply(sim_res,function(list) {
  
  nPVVs_in_data=c(27,71,8,23)
  duration_list<-Map(nPVV=nPVVs_in_data,durations=list$duration_list,function(nPVV,durations) {
    durations_filt<-durations[durations<200 & durations>0]
    if(length(durations_filt)>nPVV) {durations_filt<-sample(durations_filt,size = nPVV)}
    return(durations_filt)
  })
  
  m1 <- glm(unlist(duration_list) ~ 1, family = Gamma(link = "identity"))
  mean=coefficients(m1)
  disp=summary(m1)$dispersion
  data.frame(captured_PVV_id=list$info$captured_PVV_id,mean_lesion_duration_days=list$info$mean_lesion_duration_days,nPVV=length(unlist(duration_list)),mean=mean,disp=disp)
})%>%bind_rows()%>%
  filter(nPVV>20) #a small number of runs (~100) have failed and have very small numbers of PVVs - exclude these

##How many simulations had inadequate numbers of PVVs to match the data?
#The data has 129 PVVs - compare this to the total from each simulation
all_res%>%
  ggplot(aes(x=nPVV))+
  geom_histogram(binwidth = 3,fill="lightblue")+
  geom_vline(xintercept=129,linetype=2)+
  my_theme+
  scale_x_continuous(limits=c(0,135))+
  labs(x="Number of PVVs in ABC simulation",y="Count")

#========================================#
# Get summary statistics from the data ####
#========================================#

PVV_blood_table<-mutations%>%
  filter(Type=="PVV" & Sample_ID %in% ABC_Sample_IDs & Class=="PASS" & Sub1=="C>T" & lesion_duration<200)

#Confirm the number of mutations per sample ----
table(PVV_blood_table$Sample_ID)

#Estimate the parameters using the glm ----
m1 <- glm(unlist(PVV_blood_table$lesion_duration) ~ 1, family = Gamma(link = "identity"))
mean=coefficients(m1)
disp=summary(m1)$dispersion

#========================================#
# Perform the ABC ####
#========================================#

select_vec=1:2 #Can run the ABC with only the mean if desired
res<-abc(target=c(mean=mean,disp=disp)[select_vec],param = as.matrix(all_res[,-1])[,c('mean_lesion_duration_days'),drop=F]/365,sumstat = as.matrix(all_res[,-1])[,c('mean','disp'),drop=F][,select_vec,drop=F],tol = 0.05,method = "neuralnet")

#Write the posterior distribution for use in the PPC ----
writeLines(as.character(res$unadj.values),paste0(ABC_dir,"posterior_means2.txt"))

#Generate Fig. 5a ----
quantile(x=res$unadj.values,c(0.025,0.5,0.975))
p.posterior<-data.frame(mean_lesion_duration_years=res$unadj.values)%>%
  ggplot(aes(x=mean_lesion_duration_years))+
  geom_density(data=all_res%>%mutate(mean_lesion_duration_years=mean_lesion_duration_days/365),col="black",fill="lightblue",alpha=0.3,linewidth=0.5,aes(x=mean_lesion_duration_years))+
  geom_density(col="black",linewidth=0.5,fill="red",alpha=0.5)+
  scale_x_continuous(limits=c(0.5,5))+
  theme_classic()+
  my_theme+
  labs(x="Mean lesion duration (years)",y="Density")
ggsave(filename=paste0(plots_dir,"Fig5a.pdf"),p.posterior,width=7, height=2.5)

#Generate Extended Data Fig. 11b ----
plot2<-all_res%>%
  ggplot(aes(x=mean_lesion_duration_days/365,y=mean))+
  geom_bin2d(bins = 44) +
  scale_fill_continuous(type = "viridis") +
  guides(fill=guide_legend(position="inside",title = "Count"))+
  theme_bw()+
  my_theme+
  theme(legend.position.inside = c(0.8,0.25),legend.key.size = unit(3,"mm"))+
  labs(x="True mean lesion duration\n(Years)",y="Mean MMLD of captured PVVs")

ggsave(filename=paste0(plots_dir,"ExtDatFig11b.pdf"),plot2,width=2, height=2.5)

#Generate Extended Data Fig. 11c ----
plot3<-all_res%>%
  ggplot(aes(x=mean_lesion_duration_days/365,y=disp))+
  geom_bin2d(bins = 44) +
  scale_fill_continuous(type = "viridis") +
  guides(fill=guide_legend(position="inside",title = "Count"))+
  theme_bw()+
  my_theme+
  theme(legend.position.inside = c(0.8,0.8),legend.key.size = unit(3,"mm"))+
  labs(x="True mean lesion duration\n(Years)",y="Estimated dispersion of MMLD\ndistribution")

ggsave(filename=paste0(plots_dir,"ExtDatFig11c.pdf"),plot3,width=2, height=2.5)

#Generate Extended Data Fig. 11f----
plot4<-all_res%>%
  ggplot(aes(x=mean,y=disp))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
  geom_point(aes(x=mean,y=disp),col="red",data=data.frame(mean=mean,disp=disp))+
  theme_classic()+
  my_theme+
  labs(x="Mean MMLD", y = "Estimated dispersion of MMLD distribution")

ggsave(filename=paste0(plots_dir,"ABC_diagnostic_plots.pdf"),comb_plots,width=7, height=2.5)

#========================================#
# Infer info regarding the lesion density, and the number of lesions per cell ####
#========================================#

posterior_ids=all_res$captured_PVV_id[res$region]
uids<-sapply(sim_res,function(list) list$info$captured_PVV_id)
posterior_sim_res=sim_res[uids%in%posterior_ids]
all_posterior_info<-lapply(posterior_sim_res,function(list) {
  all_individual_info=lapply(list$individual_info,function(list2) {
    as.data.frame(list2)
  })%>%dplyr::bind_rows()
})%>%dplyr::bind_rows()

all_posterior_info_filt<-all_posterior_info%>%filter(!is.na(lesion_density))

implied_lesions_distribution=sapply(1:nrow(all_posterior_info_filt), function(k) {
  lesion_density=all_posterior_info_filt$lesion_density[k] #remember the lesion density is the number of lesions per lineage days (not years)
  mean_duration_days=all_posterior_info_filt$mean_lesion_duration_years[k]*365
  
  #Convert this into a metric for the average number of lesions in a lineage at any time
  #1. Imagine 100,000 lineage days; there will be 100,000*lesions_density lesions introduced during this period
  n_lineage_days=1e5
  n_lesions=lesion_density*n_lineage_days
  
  #2. calculate the total number of 'lesion days' from this number of lesions
  lesion_durations=rgamma(n=round(n_lesions),rate = 1/mean_duration_days,shape=1)
  
  #3. Calculate the average number of lesions over the 10000 lineage days
  n_lesions_per_lineage=sum(lesion_durations)/n_lineage_days
  
  return(n_lesions_per_lineage)
})

#Generate Extended Data Fig. 11h----
inferred_lesion_plot<-data.frame(mean_lesion_duration_years=all_posterior_info_filt$mean_lesion_duration_years,nlesions=implied_lesions_distribution)%>%
  ggplot(aes(x=nlesions))+
  geom_density(fill="lightblue")+
  scale_x_continuous(limits=c(0,25))+
  theme_classic()+
  my_theme+
  labs(x="Number of lesions per cell",y="Density")
ggsave(filename=paste0(plots_dir,"ExtDatFig11h.pdf"),inferred_lesion_plot,width=2, height=2)

quantile(x=implied_lesions_distribution,c(0.025,0.5,0.975))


#========================================#
# Performance of ABC - ability to capture 'GROUND TRUTH' ####
#========================================#

#Pre-saved version of extracted data from simulations is available on the github
GT_sim_res_file<-paste0(ABC_dir,"extracted_lesion_durations_GT.Rds")

if(file.exists(GT_sim_res_file)) {
  GT_sim_res<-readRDS(GT_sim_res_file)
} else {
  #If want to re-extract info this is the code
  GT_files=list.files(path=paste0(ABC_dir,"captured_lesions_performance_assessment"),full.names = T)
  cat(paste("Extracting information from",length(GT_files),"files"),sep="\n")
  GT_sim_res<-lapply(GT_files,function(file) {
    cat(file,sep="\n")
    x<-readRDS(file)
    temp=lapply(x, function (individual) {
      simpop<-readRDS(paste0(simpop_dir,"/simulated_populations_",individual$pop_number,".Rds"))
      tree_ss_raw<-get_subsampled_tree2(simpop,tips=which(simpop$tip.label%in%individual$tree_ss$tip.label))
      tree_mut<-get_elapsed_time_tree(tree_ss_raw,mutrateperdivision = 1,backgroundrate = 16/365)
      
      lesions<-individual$captured_lesion_info_PVVs
      this_sim_df<-lapply(1:length(lesions),function(i) {
        this_lesion_info<-lesions[[i]]
        LN<-get_ancestor_node(this_lesion_info$captured_pure_subclades[1],tree=tree_mut)
        LRN<-get_ancestor_node(tail(this_lesion_info$captured_pure_subclades,n=1),tree=tree_mut)
        LN_timing=nodeheight(tree_mut,LN)
        LRN_timing=nodeheight(tree_mut,LRN)
        MMLD=LRN_timing-LN_timing
        df<-data.frame(lesion_duration=this_lesion_info$lesion_duration,
                       selected_branch=this_lesion_info$selected_branch,
                       LN=LN,
                       LRN=LRN,
                       pure_subclade_nodes=paste(this_lesion_info$pure_subclades,collapse="-"),
                       pure_subclade_bases=paste(this_lesion_info$pure_subclade_bases,collapse = "-"),
                       LN_timing=LN_timing,
                       LRN_timing=LRN_timing,
                       MLD_days=this_lesion_info$MMLD,
                       MLD_MT=MMLD)
      })%>%dplyr::bind_rows()
      
      if(is.null(individual$n_PVVs_assessed)) {
        individual_info=list(lesion_simulation_id=paste(individual$lesion_simulation_id,collapse=","),
                             captured_PVV_id=individual$captured_PVV_id,
                             pop_number=individual$pop_number,
                             total_lesions=individual$n_lesions,
                             lineage_time=lt_vec[paste0("simulated_populations_",individual$pop_number,".Rds")],
                             mean_lesion_duration_years=individual$mean_lesion_duration_years)
      } else {
        
        individual$total_PVVs_in_set<-total_PVV_df%>%filter(pop_number==individual$pop_number & mean_lesion_duration_years==individual$mean_lesion_duration_years)%>%pull(total_PVVs_in_set)
        individual_info=list(lesion_simulation_id=paste(individual$lesion_simulation_id,collapse=","),
                             captured_PVV_id=individual$captured_PVV_id,
                             pop_number=individual$pop_number,
                             n_PVVs_assessed=individual$n_PVVs_assessed,
                             total_PVVs=individual$total_PVVs_in_set,
                             prop_PVVs_assessed=individual$n_PVVs_assessed/individual$total_PVVs_in_set,
                             total_lesions=individual$n_lesions,
                             implied_total_lesions=individual$n_lesions*individual$n_PVVs_assessed/individual$total_PVVs_in_set,
                             lineage_time=(lt_df%>%filter(pop_number==individual$pop_number)%>%pull(lineage_time)),
                             lesion_density=individual$n_lesions*individual$n_PVVs_assessed/individual$total_PVVs_in_set/(lt_df%>%filter(pop_number==individual$pop_number)%>%pull(lineage_time)),
                             mean_lesion_duration_years=individual$mean_lesion_duration_years)
      }
      
      return(list(df=this_sim_df,tree=tree_mut,individual_info=individual_info))
      
    })
    
    duration_list<-lapply(temp,function(list) {list$df$MLD_MT})
    individual_info_list<-lapply(temp,function(list) {list$individual_info})
    
    return(list(info=list(captured_PVV_id=x[[1]]$captured_PVV_id,mean_lesion_duration_days=x[[1]]$mean_lesion_duration_days),duration_list=duration_list,individual_info=individual_info_list))
  })
  
  saveRDS(GT_sim_res,file = GT_sim_res_file)
}

##Get the mean/ dispersion coeffecients from these mutation sets ----
all_GT_res<-lapply(GT_sim_res,function(list) {
  
  nPVVs_in_data=c(27,71,8,23)
  duration_list<-Map(nPVV=nPVVs_in_data,durations=list$duration_list,function(nPVV,durations) {
    durations_filt<-durations[durations<200 & durations>0]
    if(length(durations_filt)>nPVV) {durations_filt<-sample(durations_filt,size = nPVV)}
    return(durations_filt)
  })
  
  m1 <- glm(unlist(duration_list) ~ 1, family = Gamma(link = "identity"))
  mean=coefficients(m1)
  disp=summary(m1)$dispersion
  data.frame(captured_PVV_id=list$info$captured_PVV_id,mean_lesion_duration_days=list$info$mean_lesion_duration_days,mean=mean,disp=disp)
})%>%bind_rows()


## Run the ABC on each 'ground truth' run ---------------------------------

GT_ABC_res<-lapply(1:nrow(all_GT_res),function(j) {
  select_vec=1:2
  res<-abc(target=c(mean=all_GT_res$mean[j],disp=all_GT_res$disp[j])[select_vec],param = as.matrix(all_res[,-1])[,c('mean_lesion_duration_days'),drop=F]/365,sumstat = as.matrix(all_res[,-1])[,c('mean','disp'),drop=F][,select_vec,drop=F],tol = 0.05,method = "neuralnet")
  return(data.frame(mean_lesion_duration_years=all_GT_res$mean_lesion_duration_days[j]/365,posterior=res$unadj.values,posterior.nn=res$adj.values))
})

selected_levels=seq(0.5,5,0.5)
level_names=paste("Mean LD =\n",selected_levels,"yrs")

#Generate Extended Data Fig. 11d----
p.ABC_performance1<-dplyr::bind_rows(GT_ABC_res)%>%
  mutate(ground_truth_mean_LD=paste("Mean LD =\n",mean_lesion_duration_years,"yrs"))%>%
  filter(ground_truth_mean_LD%in%level_names)%>%
  ggplot(aes(x=posterior))+
  geom_density(fill="gray")+
  facet_grid(rows=vars(ground_truth_mean_LD))+
  scale_x_continuous(limits=c(0,5.1))+
  geom_vline(aes(xintercept = gt),col="red",data = data.frame(ground_truth_mean_LD=level_names,gt=selected_levels))+
  cowplot::theme_minimal_grid()+
  my_theme+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y = element_blank(),strip.text.y=element_text(angle=0))+
  labs(x="Estimated mean lesion duration (years)",y="Density")

#Generate Extended Data Fig. 11e----
p.ABC_performance2<-lapply(GT_ABC_res,function(df) {
  quantiles=quantile(x=df$posterior,c(0.025,0.5,0.975))
  return(data.frame(ground_truth_mean_LD=df$mean_lesion_duration_years[1],lowerCI=quantiles[1],median=quantiles[2],upperCI=quantiles[3]))
})%>%dplyr::bind_rows()%>%
  mutate(within_CI=ground_truth_mean_LD>=lowerCI & ground_truth_mean_LD<=upperCI)%>%
  ggplot(aes(x=ground_truth_mean_LD,y=median,ymin=lowerCI,ymax=upperCI))+
  geom_errorbar()+
  geom_abline(slope = 1,linetype = 2,col="red")+
  theme_classic()+my_theme+
  labs(x="Actual mean lesion duration (years)",y="Estimated mean lesion duration (years)")


comb_plots<-gridExtra::arrangeGrob(grobs = list(p.ABC_performance1,p.ABC_performance2),nrow=1,widths = c(1.3,1))
plot(comb_plots)

ggsave(filename=paste0(plots_dir,"ExtDatFig11d-e.pdf"),comb_plots,width=6, height=3.5)

