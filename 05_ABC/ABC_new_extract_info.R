##--------Import packages and define custom functions-----------
library(ape)
library(rsimpop)
library(dplyr)
library(ggplot2)

# Set the ggplot2 theme and palettes for plotting-----
palette2=c("#E30613","#1D71B8")
palette6=c("#72e5ef", "#074d65", "#3d99ce", "#c257d3", "#d1add5", "#61356e")
palette8=c("#2EBAED" ,"#000000" ,"#DE1C14", "#E98C7B", "#D4D2D2" ,"#ADCC54" ,"#F0D0CE","blue")
my_theme=theme_classic(base_family="Helvetica")+theme(text=element_text(size=7,family="Helvetica"),
                                                      axis.text=element_text(size=5,family="Helvetica"),
                                                      strip.text = element_text(size=6,family="Helvetica"),
                                                      legend.key.height = unit(0.25, 'cm'),
                                                      legend.key.width = unit(0.25, 'cm'),
                                                      legend.title = element_text(size=6),
                                                      legend.text = element_text(size=5,family="Helvetica"))

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

estimate_gamma_params=function(value_vec,log_rate_range=c(-2,1),shape_range=c(1,5)) {
  require(tidyr)
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rate_vec = 10^(seq(log_rate_range[1],log_rate_range[2],by=0.1)) # rho will be bounded within 1e-6 and 0.89
  shape_vec=seq(shape_range[1],shape_range[2],0.05)
  params_grid=expand_grid(rate_vec,shape_vec)
  ll = sapply(1:nrow(params_grid), function(i) {shape=params_grid$shape_vec[i]; rate=params_grid$rate_vec[i];sum(dgamma(x=value_vec, shape=shape,rate=rate,log = T))})
  return(params_grid[which.max(ll),])
}

##--------Test the total branch lengths of the different simpops: are they comparable in terms of total 'lineage time'--------
root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
plots_dir=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/Prolonged_persistence_of_DNA_lesions/Rebuttal_plots/","/lustre/scratch126/casm/team154pc/ms56/")
ABC_dir=paste0(root_dir,"/lesion_segregation/ABC_new/")
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


total_PVV_info_file=paste0(ABC_dir,"total_PVV_info.tsv")
if(file.exists(total_PVV_info_file)) {
  total_PVV_df<-read.delim(total_PVV_info_file)
} else {
  lesion_sim_files=list.files(lesion_sim_dir,full.names = T)
  total_PVV_df<-lapply(lesion_sim_files,function(file) {
    les_sim<-readRDS(file)
    data.frame(pop_number=les_sim$pop_number,mean_lesion_duration_years=les_sim$mean_lesion_duration_years,total_PVVs_in_set=length(les_sim$lesion_info_PVVs))
  })%>%dplyr::bind_rows()
  write.table(total_PVV_df,file = total_PVV_info_file,quote = F,sep="\t",row.names=F)
}


##--------Import the captured lesion files and extract the lesion & parameter data from each ---------------------------------
setwd(ABC_dir)
files=list.files(path=paste0(ABC_dir,"captured_lesion_files"),full.names = T)
uids=sapply(files,function(filename) {x<-gsub(".Rds","",filename); gsub("^.*_", "", x)})
previous_sim_res<-readRDS("extracted_lesion_durations.Rds")
already_extracted_ids<-sapply(previous_sim_res, function(list) list$info$captured_PVV_id)

reextract=F
if(reextract) {
  cat("Re-extracting info from all captured lesion files",sep="\n")
  new_files<-files
} else {
  cat("Importing info only from new captured lesion files",sep="\n")
  new_files<-files[!uids%in%already_extracted_ids]
}

if(length(new_files)>0) {
  cat(paste("Extracting information from",length(new_files),"files"),sep="\n")
  new_sim_res<-lapply(new_files,function(file) {
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
  
  if(reextract) {
    sim_res<-new_sim_res
  } else {
    sim_res<-c(previous_sim_res,new_sim_res)
  }
  
  saveRDS(sim_res,file = "extracted_lesion_durations.Rds")
} else {
  sim_res<-previous_sim_res
}

##How many simulations had inadequate numbers of PVVs to match the data?

##Get the mean/ dispersion coeffecients from these mutation sets
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
  filter(nPVV>20) #a small number of runs have failed and have very small numbers of PVVs

all_res%>%
  ggplot(aes(x=nPVV))+
  geom_histogram(binwidth = 3,fill="lightblue")+
  geom_vline(xintercept=129,linetype=2)+
  my_theme+
  scale_x_continuous(limits=c(0,135))+
  labs(x="Number of PVVs in ABC simulation",y="Count")

if(!reextract) {
  library(abc)
  library(ggplot2)
  #Read in the relevant data
  mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
  ABC_Sample_IDs=c("KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01")
  PVV_blood_table<-mutations%>%
    filter(Type=="PVV" & cat%in%c("Adult_HSPC","Chemo_HSPC","Foetal_HSPC") & Sample_ID %in% ABC_Sample_IDs & Class=="PASS" & Sub1=="C>T" & lesion_duration<200)
  table(PVV_blood_table$Sample_ID)
  m1 <- glm(unlist(PVV_blood_table$lesion_duration) ~ 1, family = Gamma(link = "identity"))
  mean=coefficients(m1)
  disp=summary(m1)$dispersion

  select_vec=1:2
  res<-abc(target=c(mean=mean,disp=disp)[select_vec],param = as.matrix(all_res[,-1])[,c('mean_lesion_duration_days'),drop=F]/365,sumstat = as.matrix(all_res[,-1])[,c('mean','disp'),drop=F][,select_vec,drop=F],tol = 0.05,method = "neuralnet")
  
  writeLines(as.character(res$unadj.values),paste0(ABC_dir,"posterior_means2.txt"))
  
  quantile(x=res$unadj.values,c(0.025,0.5,0.975))
  plot1<-data.frame(mean_lesion_duration_years=res$unadj.values)%>%
    ggplot(aes(x=mean_lesion_duration_years))+
    geom_density(data=all_res%>%mutate(mean_lesion_duration_years=mean_lesion_duration_days/365),col="black",fill="lightblue",alpha=0.3,linewidth=0.5,aes(x=mean_lesion_duration_years))+
    geom_density(col="black",linewidth=0.5,fill="red",alpha=0.5)+
    scale_x_continuous(limits=c(0.5,5))+
    theme_classic()+
    my_theme+
    labs(x="Mean lesion duration (years)",y="Density")
  
  plot2<-all_res%>%
    ggplot(aes(x=mean_lesion_duration_days/365,y=mean))+
    geom_bin2d(bins = 44) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+
    my_theme+
    labs(x="True mean lesion duration (Years)",y="Mean MMLD of captured PVVs")
  
  plot3<-all_res%>%
    ggplot(aes(x=mean_lesion_duration_days/365,y=disp))+
    geom_bin2d(bins = 44) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+
    my_theme+
    labs(x="True mean lesion duration (Years)",y="Estimated dispersion of MMLD distribution")
  
  plot4<-all_res%>%
    ggplot(aes(x=mean,y=disp))+
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
    geom_point(aes(x=mean,y=disp),col="red",data=data.frame(mean=mean,disp=disp))+
    theme_classic()+
    my_theme+
    labs(x="Mean MMLD", y = "Estimated dispersion of MMLD distribution")
  
  comb_plots<-gridExtra::arrangeGrob(grobs = list(plot2,plot3,plot1),nrow=1)
  plot(comb_plots)
  
  ggsave(filename=paste0(plots_dir,"ABC_diagnostic_plots.pdf"),comb_plots,width=7, height=2.5)
  
}

##Infer info regarding the lesion density, and the number of lesions per cell--------

posterior_ids=all_res$captured_PVV_id[res$region]
uids<-sapply(sim_res,function(list) list$info$captured_PVV_id)
posterior_sim_res=sim_res[uids%in%posterior_ids]
all_posterior_info<-lapply(posterior_sim_res,function(list) {
  all_individual_info=lapply(list$individual_info,function(list2) {
    as.data.frame(list2)
  })%>%dplyr::bind_rows()
})%>%dplyr::bind_rows()

all_posterior_info_filt<-all_posterior_info%>%filter(!is.na(lesion_density))

hist(all_posterior_info_filt$mean_lesion_duration_years)

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
quantile(implied_lesions_distribution,c(0.025,0.5,0.975))

inferred_lesion_plot<-data.frame(mean_lesion_duration_years=all_posterior_info_filt$mean_lesion_duration_years,nlesions=implied_lesions_distribution)%>%
  ggplot(aes(x=nlesions))+
  geom_density(fill="lightblue")+
  scale_x_continuous(limits=c(0,25))+
  theme_classic()+
  my_theme+
  labs(x="Number of lesions per cell",y="Density")
ggsave(filename=paste0(plots_dir,"inferred_lesion_plot.pdf"),inferred_lesion_plot,width=2, height=2)

quantile(x=implied_lesions_distribution,c(0.025,0.5,0.975))
