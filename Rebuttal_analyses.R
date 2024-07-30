####------ANALYSES FOR REBUTTAL---------
# Set the ggplot2 theme and palettes for plotting-----
library(dplyr)
library(ggplot2)
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

# Set the root directory and read in the necessary files-----

options(stringsAsFactors = FALSE)
# root_dir="~/R_work/Prolonged_persistence_of_DNA_lesions/"
# data_dir=paste0(root_dir,"Data/")
# plots_dir=paste0(root_dir,"Plots/")
# source(paste0(data_dir,"Prolonged_persistence_functions.R")) #Source functions needed for the script

root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
source(paste0(root_dir,"/my_programs/Prolonged_persistence_functions.R"))
lesion_seg_input_dir=paste0(root_dir,"/lesion_segregation/input_data")
plots_dir="~/R_work/Prolonged_persistence_of_DNA_lesions/Rebuttal_plots/"

#Create a reference table of sample names and datasets (helpful for some functions)
sample_files_full<-list.files(lesion_seg_input_dir,pattern=".txt",full.names = T)
datasets<-gsub("_samples.txt",replacement="",list.files(lesion_seg_input_dir,pattern=".txt",full.names = F))
sample_ref_df<-Map(file=sample_files_full,dataset=datasets,function(file,dataset) {
  samples<-readLines(file)
  return(data.frame(data_set=dataset,Sample_ID=samples))
})%>%dplyr::bind_rows()

###------------ASSESSING THE ROBUSTNESS OF THE TREES------------------------------------------------------------
##-------------Compare the trees with alternative phylogeny-reconstruction algorithms---------------------------
library(dplyr)
library(ape)
root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
source(paste0(root_dir,"/my_programs/Prolonged_persistence_functions.R"))
mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
sample_ref=read.csv(paste0(root_dir,"/lesion_segregation/Individual_ref.csv"))

#IMPORT THE TREES OF THE DIFFERENT TYPES
data_sets=c("MSC_BMT","MF","EM")
combined_tree_lists=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(root_dir,"/lesion_segregation/input_data/",data_set,"_samples.txt"))
  data_set_trees<-lapply(data_set_samples,function(SampleID) {
    sample_trees<-lapply(c("mpboot","iqtree","scite"),function(tree_type) {
      tree<-import_tree_file(SampleID=SampleID,data_set = data_set,type = tree_type)
      if(!is.null(tree)) {tree<-di2multi(tree)} #ensure trees are multifurcating, not dichotomous
      return(tree)
    })
    names(sample_trees)<-c("tree","iqtree_tree","scite_tree")
    return(sample_trees)
  })
  names(data_set_trees)<-data_set_samples
  return(data_set_trees)
})%>%unlist(recursive = F)


pdf("All_phylo_comparisons.pdf",height=60,width=10)
par(mfrow=c(27,4),mar=c(2, 0.5, 2, 0.5))
temp=Map(comb_trees=combined_tree_lists,SampleID=names(combined_tree_lists),function(comb_trees,SampleID) {
  tree_mod<-comb_trees$tree;tree_mod$edge.length[tree_mod$edge.length<10]<-10
  SampleID<-stringr::str_split(SampleID,pattern = "_",simplify=T)[,1]
  
  if(!is.null(comb_trees$iqtree_tree)) {
    iqtree_tree_mod<-comb_trees$iqtree_tree; iqtree_tree_mod$edge.length[iqtree_tree_mod$edge.length<10]<-10
    tree_mod<-drop.tip(tree_mod,tree_mod$tip.label[!tree_mod$tip.label%in%iqtree_tree_mod$tip.label])
    iqtree_tree_mod<-drop.tip(iqtree_tree_mod,iqtree_tree_mod$tip.label[!iqtree_tree_mod$tip.label%in%tree_mod$tip.label])
    comparePhylo_and_plot(tree_mod,iqtree_tree_mod,names=c(paste0(SampleID,": MPBoot phylogeny"),paste0(SampleID,": IQTREE phylogeny")),lwd=0.3)
  } else {
    plot.new();plot.new()
  }
  
  if(!is.null(comb_trees$scite_tree)) {
    scite_tree_mod<-comb_trees$scite_tree; scite_tree_mod$edge.length[scite_tree_mod$edge.length<10]<-10
    tree_mod<-drop.tip(tree_mod,tree_mod$tip.label[!tree_mod$tip.label%in%scite_tree_mod$tip.label])
    scite_tree_mod<-drop.tip(scite_tree_mod,scite_tree_mod$tip.label[!scite_tree_mod$tip.label%in%tree_mod$tip.label])
    comparePhylo_and_plot(tree_mod,scite_tree_mod,names=c(paste0(SampleID,": MPBoot phylogeny"),paste0(SampleID,": SCITE phylogeny")),lwd=0.3)
  } else {
    plot.new();plot.new()
  }
})
dev.off()

#Now get the Robinson-Foulds and Quartet statistics for all comparisons with the mpboot tree
comp_stats=Map(comb_trees=combined_tree_lists,SampleID=names(combined_tree_lists),function(comb_trees,SampleID) {
  cat(SampleID,sep="\n")
  iq_stats=get_tree_comp_stats(comb_trees$iqtree_tree,comb_trees$tree)
  scite_stats=get_tree_comp_stats(comb_trees$scite_tree,comb_trees$tree)
  
  colnames(iq_stats)<-paste("iqtree",colnames(iq_stats),sep="_")
  colnames(scite_stats)<-paste("scite",colnames(scite_stats),sep="_")
  
  return(cbind(data.frame(SampleID=SampleID),iq_stats,scite_stats))
  
})%>%dplyr::bind_rows()%>%tibble::remove_rownames()

readr::write_csv(comp_stats,file=paste0(root_dir,"/Rebuttal_plots/tree_comparison_stats.csv"))

#--------------MPBoot bootstrap confidence for whole tree-------------------------------------------------------
library(dplyr)
library(ape)
root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
source(paste0(root_dir,"/my_programs/Prolonged_persistence_functions.R"))
mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
sample_ref=read.csv(paste0(root_dir,"/lesion_segregation/Individual_ref.csv"))

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

df<-Map(tree=all_trees,Sample_ID=names(all_trees),function(tree,Sample_ID) {
  mpboot_conf<-as.numeric(tree$node.label)
  n_private_nodes=length(tree$tip.label)
  node_number=1:length(mpboot_conf)+n_private_nodes
  node_heights=data.frame(node=tree$edge[,1],height=nodeHeights(tree)[,1])
  mpboot_conf_df<-data.frame(Sample_ID=Sample_ID,node=node_number,bc=mpboot_conf)%>%
    left_join(node_heights,by="node")%>%
    filter(!is.na(bc))%>%
    arrange(desc(bc))%>%
    mutate(idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.))
})%>%dplyr::bind_rows()

overall_mpboot_bc<-df%>%
  ggplot(aes(x=xpos,y=bc,col=Sample_ID))+
  geom_line(alpha=0.2)+
  geom_point(size=0.2)+
  theme_classic()+
  my_theme+
  theme(legend.position = "left")+
  labs(title="All nodes",
       x="Nodes\n(by decreasing bootstrap support)",
       y="MPBoot bootstrap support value")

df_embryonic<-lapply(unique(df$Sample_ID),function(sampleID) {
  this_sample_df<-df%>%
    filter(Sample_ID==sampleID)%>%
    filter(height<50)
  if(nrow(this_sample_df)>0) {return(this_sample_df%>%mutate(idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.)))} else {return(NULL)}
})%>%dplyr::bind_rows()

df_postembryonic<-lapply(unique(df$Sample_ID),function(sampleID) {
  this_sample_df<-df%>%
    filter(Sample_ID==sampleID)%>%
    filter(height>=50)
  if(nrow(this_sample_df)>0) {return(this_sample_df%>%mutate(idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.)))} else {return(NULL)}
})%>%dplyr::bind_rows()

embryonic_mpboot_bc<-df_embryonic%>%
  ggplot(aes(x=xpos,y=bc,col=Sample_ID))+
  geom_line(alpha=0.2)+
  geom_point(size=0.2)+
  theme_classic()+
  my_theme+
  theme(legend.position = "none")+
  labs(title="Embryonic nodes",
       x="Nodes\n(by decreasing bootstrap support)",
       y="MPBoot bootstrap support value")

postembryonic_mpboot_bc<-df_postembryonic%>%
  ggplot(aes(x=xpos,y=bc,col=Sample_ID))+
  geom_line(alpha=0.2)+
  geom_point(size=0.2)+
  theme_classic()+
  my_theme+
  theme(legend.position = "none")+
  labs(title="Post-embryonic nodes",
       x="Nodes\n(by decreasing bootstrap support)",
       y="MPBoot bootstrap support value")



comb_plots<-gridExtra::arrangeGrob(grobs=list(overall_mpboot_bc,embryonic_mpboot_bc,postembryonic_mpboot_bc),nrow=1,widths = c(2,1,1))
plot(comb_plots)
ggsave(paste0(plots_dir,"mpboot_bootstrap_FULLTREE_plot.pdf"),plot=comb_plots,width=8,height=2.5)

#-------------------Overall read count boostrap confidence across all nodes-------------
#--------------Read count bootstrap confidence for trees where read count bootstraps have been performed ---------------------

# root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
# source(paste0(root_dir,"/my_programs/Prolonged_persistence_functions.R"))

root_dir="~/R_work/Prolonged_persistence_of_DNA_lesions/"
source(paste0(root_dir,"/Data/Prolonged_persistence_functions.R"))
lesion_seg_input_dir=paste0(root_dir,"/Data/input_data")
library(dplyr)
library(ape)

mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
PVV_blood_table<-mutations%>%
  filter(Type=="PVV" & cat%in%c("Adult_HSPC","Chemo_HSPC","Foetal_HSPC") &!data_set=="NW" & Class=="PASS")

#Compile the read count boostrap tables
data_sets=c("MF","EM")
combined_boostrap_confidence=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(root_dir,"/lesion_segregation/input_data/",data_set,"_samples.txt"))
  data_set_boostraps<-lapply(data_set_samples,function(SampleID) {
    
    cat(SampleID,sep="\n")
    exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1]
    filtering="vaf_filtered"
    filtering_settings="standard_rho01"
    
    #Define the file paths for the data files
    my_working_directory=ifelse(data_set=="EM",paste0(root_dir,"/Emily_benchmarking/",exp_ID,"/",filtering,"/tree_bootstraps"),paste0(root_dir,"/Marga_benchmarking/",SampleID,"/tree_bootstraps"))
    bootstrapped_tree_stats_file_path=paste0("Bootstrapped_tree_stats_",SampleID)
    
    if(!dir.exists(my_working_directory)) {stop(return(NULL))}
    setwd(my_working_directory)
    
    node_confidence_file=paste0("bootstrap_retained_",SampleID,".tsv")
    
    if(file.exists(node_confidence_file)) {
      cat("Reading in existing node_confidence_df file",sep="\n")
      node_confidence_df<-read.delim(node_confidence_file)
    } else {
      cat("Attempting to generate node confidence table from the bootstrap trees",sep="\n")
      
      #Read count bootstraps
      bootstrapped_trees_dir=paste0(my_working_directory,"/final_trees/")
      all_bootstrapped_trees =paste0(my_working_directory,"/all_trees")
      
      tree_files=list.files(path=paste0(my_working_directory,"/final_trees/"),pattern=".tree")
      if(length(tree_files)<100) {cat("Insufficient trees generated.",sep="\n");stop(return(NULL))}
      
      if(!file.exists(all_bootstrapped_trees)) {system("cat final_trees/Tree_*>all_trees")}
      
      all_boots=read.tree(all_bootstrapped_trees)
      all_boots<-di2multi(all_boots)
      
      #Load up the main tree
      tree=di2multi(import_tree_file(SampleID = SampleID,data_set = data_set,type="mpboot"))
      ROOT=1+length(tree$tip.label)
      tree_node_numbers=unique(tree$edge[,1]) #store list of the tree node numbers (not including the tips)
      
      if(!all(all_boots[[1]]$tip.label%in%tree$tip.label)) {
        cat("Dropping tips from the bootstrapped trees that aren't in the data.",sep="\n")
        cat(all_boots[[1]]$tip.label[!all_boots[[1]]$tip.label%in%tree$tip.label],sep="\n")
        all_boots<-lapply(all_boots,function(boot_tree) drop.tip(boot_tree,boot_tree$tip.label[!boot_tree$tip.label%in%tree$tip.label]))
      }
      
      boot_stats_list=list(nmuts=numeric(),
                           mean_mut_burden=numeric(),
                           nnodes_boot=numeric(),
                           shared_clades=list(),
                           new_clades=list())
      
      #Fill up the stats list from the data
      for(i in 1:length(all_boots)) {
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
      
      retained_node_read_bootstraps=sapply(tree_node_numbers,function(node) {
        sum(unlist(lapply(boot_stats_list$shared_clades,function(x) node%in%as.numeric(x))))/length(boot_stats_list$nmuts)
      })
      
      #Save the data as to the proportion of bootstraps in which the node is retained
      node_confidence_df<-data.frame(node=tree_node_numbers,retained=retained_node_read_bootstraps)%>%
        filter(node!=ROOT)%>%
        mutate(SampleID=SampleID,.before=1)
      
      write.table(node_confidence_df,file=node_confidence_file,quote=F,sep="\t",row.names = F)
    }
    
    return(node_confidence_df)
    
  })
  names(data_set_boostraps)<-data_set_samples
  return(data_set_boostraps)
})%>%unlist(recursive = F)

# df<-Map(tree=all_trees,Sample_ID=names(all_trees),function(tree,Sample_ID) {
#   mpboot_conf<-as.numeric(tree$node.label)
#   n_private_nodes=length(tree$tip.label)
#   node_number=1:length(mpboot_conf)+n_private_nodes
#   node_heights=data.frame(node=tree$edge[,1],height=nodeHeights(tree)[,1])
#   mpboot_conf_df<-data.frame(Sample_ID=Sample_ID,node=node_number,bc=mpboot_conf)%>%
#     left_join(node_heights,by="node")%>%
#     filter(!is.na(bc))%>%
#     arrange(desc(bc))%>%
#     mutate(idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.))
# })%>%dplyr::bind_rows()

rc_bootstrap_df<-Map(df=combined_boostrap_confidence,Sample_ID=names(combined_boostrap_confidence),function(df,Sample_ID) {
  if(is.null(df)) {stop(return(NULL))}
  tree<-all_trees[[Sample_ID]]
  node_heights=data.frame(node=tree$edge[,1],height=nodeHeights(tree)[,1])%>%
    filter(!duplicated(.))
  df%>%
    left_join(node_heights,by="node")%>%
    arrange(desc(retained))%>%
    mutate(Sample_ID=Sample_ID,idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.))%>%
    dplyr::rename("bc"=retained)
})%>%dplyr::bind_rows()

overall_rc_bc<-rc_bootstrap_df%>%
  ggplot(aes(x=xpos,y=bc,col=Sample_ID))+
  geom_line(alpha=0.2)+
  geom_point(size=0.2)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position = "left")+
  labs(title="All nodes",
       x="Nodes\n(by decreasing bootstrap support)",
       y="Read count bootstrap support")

df_embryonic<-lapply(unique(rc_bootstrap_df$Sample_ID),function(sampleID) {
  this_sample_df<-rc_bootstrap_df%>%
    filter(Sample_ID==sampleID)%>%
    filter(height<50)
  if(nrow(this_sample_df)>0) {return(this_sample_df%>%mutate(idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.)))} else {return(NULL)}
})%>%dplyr::bind_rows()

df_postembryonic<-lapply(unique(rc_bootstrap_df$Sample_ID),function(sampleID) {
  this_sample_df<-rc_bootstrap_df%>%
    filter(Sample_ID==sampleID)%>%
    filter(height>=50)
  if(nrow(this_sample_df)>0) {return(this_sample_df%>%mutate(idx=1:nrow(.))%>%mutate(xpos=idx/nrow(.)))} else {return(NULL)}
})%>%dplyr::bind_rows()

embryonic_rc_bc<-df_embryonic%>%
  ggplot(aes(x=xpos,y=bc,col=Sample_ID))+
  geom_line(alpha=0.2)+
  geom_point(size=0.2)+
  theme_classic()+
  my_theme+
  theme(legend.position = "none")+
  labs(title="Embryonic nodes",
       x="Nodes\n(by decreasing bootstrap support)",
       y="Read count bootstrap support value")

postembryonic_rc_bc<-df_postembryonic%>%
  ggplot(aes(x=xpos,y=bc,col=Sample_ID))+
  geom_line(alpha=0.2)+
  geom_point(size=0.2)+
  theme_classic()+
  my_theme+
  theme(legend.position = "none")+
  labs(title="Post-embryonic nodes",
       x="Nodes\n(by decreasing bootstrap support)",
       y="Read count bootstrap support value")

comb_rc_bc_plots<-gridExtra::arrangeGrob(grobs=list(overall_rc_bc,embryonic_rc_bc,postembryonic_rc_bc),nrow=1,widths = c(1.5,1,1))
plot(comb_rc_bc_plots)
ggsave(paste0(plots_dir,"rc_bootstrap_FULLTREE_plot.pdf"),plot=comb_rc_bc_plots,width=8,height=2.5)

### Look at genotype maps of the blood phylogenies----------------------------
library(pheatmap)

#Convenience function
get_highly_shared_muts=function(details,tree,min_samples=3) {
  require(dplyr)
  highly_shared_logical<-sapply(tree$edge[,2],function(node) length(getTips(tree,node))>min_samples)
  highly_shared_nodes<-tree$edge[,2][highly_shared_logical]
  details%>%filter(node%in%highly_shared_nodes)%>%pull(mut_ref)
}

data_sets=c("MSC_BMT","MF","EM")
combined_gt_maps=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(root_dir,"Data/input_data/",data_set,"_samples.txt"))
  data_set_gt<-lapply(data_set_samples,function(SampleID) {
    cat(SampleID,sep="\n")
    load(get_file_paths_and_project(Sample_ID=SampleID,dataset=data_set,input_data_dir = lesion_seg_input_dir)$filtered_muts_path)
    tree<-read.tree(get_file_paths_and_project(Sample_ID=SampleID,dataset=data_set,input_data_dir = lesion_seg_input_dir)$tree_file_path)
    details<-filtered_muts$COMB_mats.tree.build$mat
    gt<-filtered_muts$Genotype_shared_bin
    highly_shared_muts<-get_highly_shared_muts(details,tree,min_samples = 5)
    mut_select<-highly_shared_muts[highly_shared_muts%in%rownames(gt)]
    gt<-gt[mut_select,]
    
    #Reduce the size to a maximum of 1000 mutations size
    max_mut=1000
    if(dim(gt)[1]>max_mut) {gt<-gt[sample(1:nrow(gt),size=max_mut,replace = F),]}
    
    return(gt)
  })
  names(data_set_gt)<-data_set_samples
  return(data_set_gt)
})%>%unlist(recursive = F)

pdf("temp.pdf")
selected_individuals=c(3,4,11,13,14,16,18,24,27)
gt_plots=Map(gt=combined_gt_maps[selected_individuals],SampleID=names(combined_gt_maps)[selected_individuals],function(gt,SampleID) {
  p<-ggplotify::as.ggplot(pheatmap(gt,
                                   main=SampleID,
                                   legend=T,
                                   show_colnames=F,
                                   show_rownames=F,
                                   legend_breaks = c(0.8,0.5,0.2),
                                   legend_labels = c("+","?","-"),
                                   color = c("#08519C","light grey","#CB181D")))
  return(p)
})
dev.off()

comb_plots=cowplot::plot_grid(plotlist = gt_plots, ncol = 3)
ggsave(paste0(plots_dir,"selectedl_shared_genotypes.pdf"),comb_plots,width = 10,height=12)

#--------------Running the disagreement scores------------------------------------------------------------------

disagreement_score=function(binary_mut_mat,max_mut=5000) {
  if(dim(binary_mut_mat)[1]>max_mut) {
    binary_mut_mat<-binary_mut_mat[sample(1:nrow(binary_mut_mat),size=max_mut,replace = F),]
  }
  c=matrix(as.logical(binary_mut_mat),ncol=ncol(binary_mut_mat)) #c is a binary matrix with 0.5's rounded up to 1's
  d<-binary_mut_mat;d[d==0.5]<-0;d<-matrix(as.logical(d),ncol=ncol(d)) #d is a binary matrix with 0.5's rounded down to 0's
  e_exc=(d)%*%t(!c)
  e_inc=(d)%*%t(d) #Need to use the version where 0.5 is cooerced to 0
  min_exc=matrix(mapply("min",e_exc,t(e_exc)),ncol=ncol(e_exc))
  out=matrix(mapply("min",min_exc,e_inc),ncol=ncol(e_exc))
  mean(out)
}

disagreement_scores_file=paste0(root_dir,"/all_disagreement_scores.Rds")
if(file.exists(disagreement_scores_file)) {
  disagreement_scores<-readRDS(disagreement_scores_file)
} else {
  data_sets=c("MSC_BMT","MF","EM","MSC_fetal")
  diagreement_scores=lapply(data_sets,function(data_set) {
    data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
    data_set_disagreement_scores<-lapply(data_set_samples,function(sample_ID) {
      cat(sample_ID,sep="\n")
      tree<-import_tree_file(SampleID=sample_ID,data_set = data_set,type = "mpboot")
      tree<-drop.tip(tree,tip = "Ancestral")
      cat("Importing the filtered muts files")
      filtered_muts=import_muts_file(SampleID=sample_ID,data_set=data_set)
      gt<-filtered_muts$Genotype_shared_bin[,tree$tip.label]
      cat("Calculating disagreement score of the data")
      dat=disagreement_score(gt)
      nrand=20; max_mut=5000
      cat("Calculating disagreement score of random shuffles of the data")
      rand=sapply(1:nrand,function(i) {print(i);gt_sub<-gt[sample(1:nrow(gt),size=max_mut,replace = F),];gt_shuffled<-apply(gt_sub,1,FUN = base::sample);disagreement_score(gt_shuffled,max_mut = max_mut)})
      sample_results=list(data=dat,rand=rand)
      return(sample_results)
    })
    names(data_set_disagreement_scores)<-data_set_samples
    return(data_set_disagreement_scores)
  })%>%unlist(recursive = F)
  
  #Save the output, as this takes a while
  saveRDS(diagreement_scores,file = disagreement_scores_file)
}

all_ds_combined=Map(list=disagreement_scores,Sample_ID=names(disagreement_scores),function(list,Sample_ID) {
  bind_rows(data.frame(ds=list$data,type="data"),
            data.frame(ds=list$rand,type="random shuffles"))%>%
    mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1],.before=1)
  })%>%bind_rows()

breaks <- 10^(-2:2)
minor_breaks <- NULL
all_ds_plot<-all_ds_combined%>%
  filter(type=="random shuffles")%>%
  ggplot(aes(x=Sample_ID,y=ds,col=type))+
  #geom_point(alpha=0.5)+
  geom_jitter(size=0.15,width = 0.1,alpha=0.3)+
  geom_point(data=all_ds_combined%>%filter(type=="data"),size=0.3)+
  scale_y_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),panel.grid.major.x = element_line())+
  labs(x="Individual",y="Disagreement score",col="Type")+
  annotation_logticks(sides="l",short=unit(0.05,"cm"),mid=unit(0.05,"cm"),long=unit(0.1,"cm"))
  
ggsave(filename=paste0(plots_dir,"/all_disagreement_scores.pdf"),width=7,height=2.5)

###------------ASSESSING THE ROBUSTNESS OF THE PVV NODES-----------------------------------------------------
#---------------MPBoot bootstrap confidence for all the PVV nodes---------------------
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

#Use the original version of the mpboot tree
all_trees$KX003_5_01<-read.tree(paste0(lesion_seg_input_dir,"/EM/tree_KX003_5_01_standard_rho01_ORIGINAL.tree"))

mpboot_bc_df=lapply(unique(PVV_blood_table$Sample_ID),function(sample_ID) {
  cat(sample_ID,sep="\n")
  this_sample_PVVs<-PVV_blood_table%>%filter(Sample_ID==sample_ID)
  this_tree<-all_trees[[sample_ID]]
  n_private_nodes=length(this_tree$tip.label) #Need this to work out the index of the node label relating to the relevant node
  
  sample_bc<-lapply(1:nrow(this_sample_PVVs),function(i) {
    LN<-this_sample_PVVs$lesion_node[i]
    LRN<-this_sample_PVVs$lesion_repair_node[i]
    
    data.frame(Sample_ID=sample_ID,mut_ref=this_sample_PVVs$mut_ref1[i],LN_bc=as.numeric(this_tree$node.label[LN-n_private_nodes]),LRN_bc=as.numeric(this_tree$node.label[LRN-n_private_nodes]))
    
  })%>%dplyr::bind_rows()
  
  return(sample_bc)
})%>%dplyr::bind_rows()

#Some summary proportion regarding the bootstrap support
sum(mpboot_bc_df$LN_bc>95 & mpboot_bc_df$LRN_bc>95)/nrow(mpboot_bc_df)
sum(mpboot_bc_df$LN_bc<80 | mpboot_bc_df$LRN_bc<80)/nrow(mpboot_bc_df)

#Some lower confidence PVVs defined by <80 bootstrap support for either the lesion node, or lesion repair node
mpboot_flagged_mut_refs<-mpboot_bc_df%>%
  filter(LN_bc<80|LRN_bc<80)%>%
  pull(mut_ref)

Sample_cols<-c("#88e99a", "#8e0049", "#24a475", "#115e41", "#64d4fd",
               "#34466d", "#d0a8f9", "#7430a3", "#b1e632", "#673d17",
               "#51f310", "#dc5dd8", "#799d10", "#3f16f9", "#cdd9b8",
               "#e72525", "#21a708", "#1288da", "#f4d403", "#f4327e",
               "#f7931e", "#809776", "#fd9f9f", "#fa1bfc")

names(Sample_cols)<-stringr::str_split(unique(mpboot_bc_df$Sample_ID),pattern="_",simplify=T)[,1]
bootstrap_confidence_plot<-mpboot_bc_df%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=LN_bc,y=LRN_bc,col=Sample_ID))+
  scale_x_continuous(limits=c(0,105))+
  scale_y_continuous(limits=c(0,105))+
  scale_color_manual(values=Sample_cols)+
  labs(x="Lesion node\nbootstrap support",y="Lesion repair node\nbootstrap support")+
  geom_jitter(alpha=0.3,width=2,height=2,size=0.5)+
  theme_classic()+
  my_theme+
  theme(legend.position="left")

ggsave("bootstrap_confidence_plot.pdf",plot=bootstrap_confidence_plot,width=4.5,height=2.5)

##Visualize the flagged up PVVs
unique_node_pair_df<-PVV_blood_table%>%dplyr::filter(mut_ref1%in%mpboot_flagged_mut_refs)%>%dplyr::select(Sample_ID,lesion_node,lesion_repair_node)%>%unique()

pdf("PVVs_poor_mpboot_bs_support.pdf",width=12,height=20)
par(mfrow=c(5,2))
temp=sapply(1:nrow(unique_node_pair_df),function(i) {
  cat(i,sep="\n")
  SampleID=unique_node_pair_df$Sample_ID[i]
  tree<-all_trees[[SampleID]]
  LN=unique_node_pair_df$lesion_node[i]; LRN=unique_node_pair_df$lesion_repair_node[i]
  
  info=get_file_paths_and_project(dataset=sample_ref_df%>%filter(Sample_ID==SampleID)%>%pull(data_set),Sample_ID=SampleID)
  load(info$filtered_muts_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  
  ##Plot the original PVV on the mpboot tree
  sub_tree<-extract.clade(phy=tree,node=LN)
  sub_tree<-squash_tree(sub_tree,cut_off = min(200,20+max(nodeHeights(sub_tree)[sub_tree$edge[1,]])))
  PVVs<-mutations%>%filter(Sample_ID==SampleID & lesion_node==LN & lesion_repair_node==LRN)%>%pull(mut_ref1)
  cat(PVVs,sep="\n")
  #confirmatory_mutations<-details%>%filter(node%in%c(LRN,402,403,397,get_direct_daughters(LN,tree)))%>%arrange(node)%>%pull(mut_ref)
  confirmatory_mutations<-details%>%filter(node==LRN)%>%pull(mut_ref)
  hm=create_mutation_heatmap(mut_refs1=PVVs,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=sub_tree)
  
  print(filtered_muts$COMB_mats.tree.build$NR[c(PVVs,confirmatory_mutations),sub_tree$tip.label]) #confirm the coverage
  
  #The plots
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=1,title = paste(SampleID,":",paste(PVVs,collapse = ",")))
  add_mut_heatmap_PVV(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.6)
})
dev.off()

#-------------------Readcount bootstrap read confidence for PVV nodes-------------
combined_bootstrap_confidence_df<-combined_boostrap_confidence%>%dplyr::bind_rows()

#From this dataframe, create the bootstrap confidence dataframe for each PVV
rc_bc_df<-left_join(PVV_blood_table%>%dplyr::select(Sample_ID,mut_ref1,lesion_node,lesion_repair_node),combined_bootstrap_confidence_df,by=c("Sample_ID"="SampleID","lesion_node"="node"))%>%
  dplyr::rename("LN_bc"=retained)%>%
  left_join(combined_bootstrap_confidence_df,by=c("Sample_ID"="SampleID","lesion_repair_node"="node"))%>%
  dplyr::rename("LRN_bc"=retained)%>%
  mutate(LN_bc=100*LN_bc,LRN_bc=100*LRN_bc)%>%
  filter(!is.na(LN_bc))

sum(rc_bc_df$LN_bc==100 & rc_bc_df$LRN_bc==100)/nrow(rc_bc_df)
sum(rc_bc_df$LN_bc<80 | rc_bc_df$LRN_bc<80)/nrow(rc_bc_df)

##List of the lower confidence PVVs defined by <80 bootstrap support for either the lesion node, or lesion repair node
rc_flagged_mut_refs<-rc_bc_df%>%
  filter(LN_bc<80|LRN_bc<80)%>%
  pull(mut_ref1)

readcount_bootstrap_confidence_plot<-rc_bc_df%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=LN_bc,y=LRN_bc,col=Sample_ID))+
  scale_x_continuous(limits=c(0,105))+
  scale_y_continuous(limits=c(0,105))+
  scale_color_manual(values=Sample_cols)+
  labs(x="Lesion node\nbootstrap support",y="Lesion repair node\nbootstrap support")+
  geom_jitter(alpha=0.2,width=2,height=2,size=0.5)+
  theme_classic()+
  my_theme+
  theme(legend.position="none")

ggsave("readcount_bootstrap_confidence_plot.pdf",plot=readcount_bootstrap_confidence_plot,width=3.5,height=2.5)

bootstrap_confidence_combined_line_plot<-rc_bc_df%>%
  mutate(min_bc=pmin(LN_bc,LRN_bc))%>%
  arrange(desc(min_bc))%>%
  mutate(pos=1:nrow(.),type="Read count")%>%
  mutate(xpos=pos/nrow(.))%>%
  bind_rows(mpboot_bc_df%>%
              mutate(min_bc=pmin(LN_bc,LRN_bc))%>%
              arrange(desc(min_bc))%>%
              mutate(pos=1:nrow(.),type="MPBoot")%>%
              mutate(xpos=pos/nrow(.)))%>%
  ggplot(aes(x=xpos,y=min_bc,col=type))+
  geom_point(size=0.1,alpha=0.5)+
  geom_line(alpha=0.5)+
  theme_classic()+
  my_theme+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  labs(x="PVVs (arranged by\ndecreasing bootstrap support)",
       y="Bootstrap support for lesion/ lesion repair node\n(Minimum of the two nodes)",
       col="Bootstrap\ntype")

ggsave("bootstrap_confidence_combined_line_plot.pdf",plot=bootstrap_confidence_combined_line_plot,width=3.5,height=2.5)

combined_PVV_bs_support<-gridExtra::arrangeGrob(grobs =list(bootstrap_confidence_plot,
                                                            readcount_bootstrap_confidence_plot,
                                                            bootstrap_confidence_combined_line_plot),
                                                 nrow=1,
                                                 widths=c(1.5,1,1.2))

ggsave(paste0(plots_dir,"bootstrap_confidence_all_combined.pdf"),plot=combined_PVV_bs_support,width=8,height=2.5)


unique_node_pair_df<-PVV_blood_table%>%dplyr::filter(mut_ref1%in%rc_flagged_mut_refs)%>%dplyr::select(Sample_ID,lesion_node,lesion_repair_node)%>%unique()

pdf("PVVs_poor_readcount_bs_support.pdf",width=12,height=10)
par(mfrow=c(2,2))
temp=sapply(1:nrow(unique_node_pair_df),function(i) {
  cat(i,sep="\n")
  SampleID=unique_node_pair_df$Sample_ID[i]
  tree<-all_trees[[SampleID]]
  LN=unique_node_pair_df$lesion_node[i]; LRN=unique_node_pair_df$lesion_repair_node[i]
  
  info=get_file_paths_and_project(dataset=sample_ref_df%>%filter(Sample_ID==SampleID)%>%pull(data_set),Sample_ID=SampleID)
  load(info$filtered_muts_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  
  ##Plot the original PVV on the mpboot tree
  sub_tree<-extract.clade(phy=tree,node=LN)
  sub_tree<-squash_tree(sub_tree,cut_off = min(200,20+max(nodeHeights(sub_tree)[sub_tree$edge[1,]])))
  PVVs<-mutations%>%filter(Sample_ID==SampleID & lesion_node==LN & lesion_repair_node==LRN)%>%pull(mut_ref1)
  cat(PVVs,sep="\n")
  #confirmatory_mutations<-details%>%filter(node%in%c(LRN,402,403,397,get_direct_daughters(LN,tree)))%>%arrange(node)%>%pull(mut_ref)
  confirmatory_mutations<-details%>%filter(node==LRN)%>%pull(mut_ref)
  hm=create_mutation_heatmap(mut_refs1=PVVs,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=sub_tree)
  
  print(filtered_muts$COMB_mats.tree.build$NR[c(PVVs,confirmatory_mutations),sub_tree$tip.label]) #confirm the coverage
  
  #The plots
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=1,title = paste(SampleID,":",paste(PVVs,collapse = ",")))
  add_mut_heatmap_PVV(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.6)
})
dev.off()

#---------------Now do a similar analysis but looking at whether clades are identical within IQTree and SCITE trees---------------------

###Test if the clades defined by the LN and LRN are identical between tree-building methods----------

retained_in_alternative_tree_building<-lapply(names(combined_tree_lists),function(SampleID) {
  cat(SampleID,sep="\n")
  #Load in the trees
  mpboot_tree<-combined_tree_lists[[SampleID]]$tree
  iqtree_tree<-combined_tree_lists[[SampleID]]$iqtree_tree
  scite_tree<-combined_tree_lists[[SampleID]]$scite_tree
  
  #Do node comparisons
  if(!is.null(iqtree_tree)) {iqtree_tree<-keep.tip(iqtree_tree,iqtree_tree$tip.label[iqtree_tree$tip.label%in%mpboot_tree$tip.label]);iq_comp<-compare_nodes(mpboot_tree,iqtree_tree); iq_shared_clades<-as.numeric(iq_comp[,1])} else {iq_shared_clades<-NA}
  if(!is.null(scite_tree)) {scite_tree<-keep.tip(scite_tree,scite_tree$tip.label[scite_tree$tip.label%in%mpboot_tree$tip.label]);scite_comp<-compare_nodes(mpboot_tree,scite_tree); scite_shared_clades<-as.numeric(scite_comp[,1])} else {scite_shared_clades<-NA}
  
  #Extract parameters of interest
  this_sample_PVVs<-mutations%>%filter(Sample_ID==SampleID & Type == "PVV" & Class=="PASS")
  retained_df<-this_sample_PVVs%>%
    mutate(LN_retained_iq=lesion_node%in%iq_shared_clades,LRN_retained_iq=lesion_repair_node%in%iq_shared_clades)%>%
    mutate(LN_retained_scite=lesion_node%in%scite_shared_clades,LRN_retained_scite=lesion_repair_node%in%scite_shared_clades)%>%
    dplyr::select(Sample_ID,Chrom_pos,mut_ref1,lesion_node,lesion_repair_node,LN_retained_iq,LRN_retained_iq,LN_retained_scite,LRN_retained_scite)
  
  if(is.null(iqtree_tree)) {retained_df$LN_retained_iq<-NA;retained_df$LRN_retained_iq<-NA}
  if(is.null(scite_tree)) {retained_df$LN_retained_scite<-NA;retained_df$LRN_retained_scite<-NA}
  
  return(retained_df)
  
})%>%dplyr::bind_rows()

nPVV_assessed<-mutations%>%filter(data_set%in%data_sets & Type == "PVV" & Class=="PASS")%>%nrow()
imperfect_match_df<-retained_in_alternative_tree_building%>%filter(!LN_retained_iq|!LRN_retained_iq)

#Compare with the set of PVVs called in these samples by calling on the IQTREE set---------------

library(stringr)
library(dplyr)
library(ggplot2)
library(Rsamtools)

root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
lesion_seg_input_dir=paste0(root_dir,"/lesion_segregation/input_data")
lesion_seg_output_dir=paste0(root_dir,"/lesion_segregation/output_iqtree")

#Set data file paths
iqtree_data_sets=c("EMiq","MSC_BMTiq","MFiq")
out_list=lapply(iqtree_data_sets,function(data_set) {
  print(data_set)
  files=list.files(paste0(lesion_seg_output_dir,"/",data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
  print(files)
  mut_tables_list=lapply(files,readr::read_delim,delim="\t",col_types=c("cccccccciiicciiiiiiccccc"))
  mut_tables_list=lapply(mut_tables_list,function(df) {
    df$colonies_per_negative_subclade<-as.character(df$colonies_per_negative_subclade)
    df$depth_per_negative_subclade<-as.character(df$depth_per_negative_subclade)
    df$Ref<-as.character(df$Ref);df$Alt1<-as.character(df$Alt1);df$Alt2<-as.character(df$Alt2)
    return(df)})
  mut_tables_df=dplyr::bind_rows(mut_tables_list)
  mut_tables_df$data_set=data_set
  return(mut_tables_df)
})
iqtree_mutations=dplyr::bind_rows(out_list)


imperfect_match_df<-imperfect_match_df%>%
  mutate(called_in_iq=ifelse(mut_ref1%in%iqtree_mutations$mut_ref1,T,F))%>%
  mutate(called_in_iq=ifelse(Sample_ID%in%iqtree_mutations$Sample_ID,called_in_iq,NA))

unique_node_pair_df<-imperfect_match_df%>%filter(!LN_retained_iq|!LRN_retained_iq)%>%dplyr::filter(!called_in_iq)%>%dplyr::select(Sample_ID,lesion_node,lesion_repair_node)%>%unique()

alt_treebuilding_flagged_mut_refs<-imperfect_match_df%>%dplyr::filter(!called_in_iq)%>%pull(mut_ref1)
all_low_confidence_df<-data.frame(mut_ref=unique(c(mpboot_flagged_mut_refs,rc_flagged_mut_refs,alt_treebuilding_flagged_mut_refs)))%>%left_join(PVV_blood_table%>%select(Sample_ID,mut_ref1),by=c("mut_ref"="mut_ref1"))
all_low_confidence_df$mpboot<-all_low_confidence_df$mut_ref%in%mpboot_flagged_mut_refs
all_low_confidence_df$rc<-all_low_confidence_df$mut_ref%in%rc_flagged_mut_refs
all_low_confidence_df$alt_phylo<-all_low_confidence_df$mut_ref%in%alt_treebuilding_flagged_mut_refs
readr::write_csv(all_low_confidence_df,file=paste0(root_dir,"/Rebuttal_plots/low_confidence_PVVs.csv"))

#Visualize these mutations - could they be explained by an alternative phylogeny??
pdf("Inconsistent_PVVs.pdf",width=5,height=10)
par(mfrow=c(8,2))
temp=sapply(1:nrow(unique_node_pair_df),function(i) {
  cat(i,sep="\n")
  SampleID=unique_node_pair_df$Sample_ID[i]
  
  tree<-combined_tree_lists[[SampleID]]$tree
  LN=unique_node_pair_df$lesion_node[i]; LRN=unique_node_pair_df$lesion_repair_node[i]
  LRN_ancestors<-get_ancestral_nodes(LRN,tree$edge)
  LN_ancestors<-get_ancestral_nodes(LN,tree$edge)
  LP_nodes=LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
  
  info=get_file_paths_and_project(dataset=sample_ref_df%>%filter(Sample_ID==SampleID)%>%pull(data_set),Sample_ID=SampleID)
  load(info$filtered_muts_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  
  ##Plot the original PVV on the mpboot tree
  sub_tree<-extract.clade(phy=tree,node=LN)
  sub_tree<-squash_tree(sub_tree,cut_off = min(200,20+max(nodeHeights(sub_tree)[sub_tree$edge[1,]])))
  PVVs<-mutations%>%filter(Sample_ID==SampleID & lesion_node==LN & lesion_repair_node==LRN)%>%pull(mut_ref1)
  cat(PVVs,sep="\n")
  #confirmatory_mutations<-details%>%filter(node%in%c(LRN,402,403,397,get_direct_daughters(LN,tree)))%>%arrange(node)%>%pull(mut_ref)
  confirmatory_mutations<-details%>%filter(node%in%LP_nodes)%>%pull(mut_ref)
  hm=create_mutation_heatmap(mut_refs1=PVVs,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=sub_tree)
  
  print(filtered_muts$COMB_mats.tree.build$NR[c(PVVs,confirmatory_mutations),sub_tree$tip.label]) #confirm the coverage
  
  #The plots
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=1,title = paste(SampleID,":",paste(PVVs,collapse = ",")))
  add_mut_heatmap_PVV(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.6)
  
  ##Now on the altered iqtree tree
  iqtree_tree<-combined_tree_lists[[SampleID]]$iqtree_tree
  iqtree_tree<-keep.tip(iqtree_tree,iqtree_tree$tip.label[iqtree_tree$tip.label%in%tree$tip.label])
  
  pos_samples=getTips(tree,LN); pos_samples<-pos_samples[pos_samples%in%iqtree_tree$tip.label]
  iqtree_LN<-find_latest_acquisition_node(iqtree_tree,pos_samples=pos_samples)
  iq_sub_tree<-extract.clade(phy=iqtree_tree,node=iqtree_LN)
  iq_sub_tree<-squash_tree(iq_sub_tree,cut_off = min(200,20+max(nodeHeights(iq_sub_tree)[iq_sub_tree$edge[1,]])))
  print(filtered_muts$COMB_mats.tree.build$NR[c(PVVs,confirmatory_mutations),iq_sub_tree$tip.label])
  hm=create_mutation_heatmap(mut_refs1=PVVs,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=iq_sub_tree)
  
  #The plots
  iq_sub_tree=plot_tree(iq_sub_tree,cex.label=0,vspace.reserve=1,title = paste(SampleID,":",paste(PVVs,collapse = ",")))
  add_mut_heatmap_PVV(tree=iq_sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.6)
  
})
dev.off()

#Separately plot the mutations highlighted by the bootstrap support but not by the tree-building
low_bs_unique_node_pair_df<-PVV_blood_table%>%
  dplyr::filter(mut_ref1%in%(all_low_confidence_df%>%filter(!alt_phylo)%>%pull(mut_ref)))%>%
  dplyr::select(Sample_ID,lesion_node,lesion_repair_node)%>%
  unique()

pdf(paste0("low_bs_PVVs.pdf"),width=5,height=5)
par(mfrow=c(2,2))
temp=sapply(1:nrow(low_bs_unique_node_pair_df),function(i) {
  cat(i,sep="\n")
  SampleID=low_bs_unique_node_pair_df$Sample_ID[i]
  LN=low_bs_unique_node_pair_df$lesion_node[i]
  LRN=low_bs_unique_node_pair_df$lesion_repair_node[i]
  
  info=get_file_paths_and_project(dataset=sample_ref_df%>%filter(Sample_ID==SampleID)%>%pull(data_set),Sample_ID=SampleID)
  tree<-read.tree(info$tree_file_path)
  load(info$filtered_muts_path); details<-filtered_muts$COMB_mats.tree.build$mat
  sub_tree<-extract.clade(phy=tree,node=LN)
  sub_tree<-squash_tree(sub_tree,cut_off = min(200,20+max(nodeHeights(sub_tree)[sub_tree$edge[1,]])))
  PVVs<-mutations%>%filter(Sample_ID==SampleID & lesion_node==LN & lesion_repair_node==LRN)%>%pull(mut_ref1)
  if(!"node"%in% colnames(details)) {stop(return(NULL))}
  LRN_ancestors<-get_ancestral_nodes(LRN,tree$edge)
  LN_ancestors<-get_ancestral_nodes(LN,tree$edge)
  LP_nodes<-LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
  confirmatory_mutations<-details%>%filter(node%in%LP_nodes)%>%arrange(node)%>%pull(mut_ref)
  combined_muts=c(PVVs,confirmatory_mutations)
  print(filtered_muts$COMB_mats.tree.build$NR[combined_muts,sub_tree$tip.label]) #confirm the coverage
  
  ##Plot the original PVV on the mpboot tre
  hm=create_mutation_heatmap(mut_refs1=PVVs,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=sub_tree)
  
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=1,title = paste(SampleID,":",paste(PVVs,collapse = ",")))
  add_mut_heatmap_PVV(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.3)
  
})
dev.off()


##Manually confirm some of the PVVs with low numbers of phylogeny-confirming mutations---------------

low_conf_PVVs<-PVV_blood_table%>%
  filter(Type=="PVV" & cat=="Adult_HSPC" & Class=="PASS" &lesion_duration<5)

#11/19 PVVs in this category were already flagged in previous approaches
sum(low_conf_PVVs$mut_ref1%in%all_low_confidence_df$mut_ref)

low_conf_unique_node_pair_df<-low_conf_PVVs%>%
  dplyr::filter(!mut_ref1%in%all_low_confidence_df$mut_ref)%>%
  dplyr::select(Sample_ID,lesion_node,lesion_repair_node)%>%
  unique()

pdf(paste0("low_number_of_confirmatory_mutation_PVVs.pdf"),width=5,height=10)
par(mfrow=c(4,2))
temp=sapply(1:nrow(low_conf_unique_node_pair_df),function(i) {
  cat(i,sep="\n")
  SampleID=low_conf_unique_node_pair_df$Sample_ID[i]
  LN=low_conf_unique_node_pair_df$lesion_node[i]
  LRN=low_conf_unique_node_pair_df$lesion_repair_node[i]
  
  info=get_file_paths_and_project(dataset=sample_ref_df%>%filter(Sample_ID==SampleID)%>%pull(data_set),Sample_ID=SampleID)
  tree<-read.tree(info$tree_file_path)
  load(info$filtered_muts_path); details<-filtered_muts$COMB_mats.tree.build$mat
  sub_tree<-extract.clade(phy=tree,node=LN)
  sub_tree<-squash_tree(sub_tree,cut_off = min(200,20+max(nodeHeights(sub_tree)[sub_tree$edge[1,]])))
  PVVs<-mutations%>%filter(Sample_ID==SampleID & lesion_node==LN & lesion_repair_node==LRN)%>%pull(mut_ref1)
  if(!"node"%in% colnames(details)) {stop(return(NULL))}
  LRN_ancestors<-get_ancestral_nodes(LRN,tree$edge)
  LN_ancestors<-get_ancestral_nodes(LN,tree$edge)
  LP_nodes<-LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
  confirmatory_mutations<-details%>%filter(node%in%LP_nodes)%>%pull(mut_ref)
  combined_muts=c(PVVs,confirmatory_mutations)
  print(filtered_muts$COMB_mats.tree.build$NR[combined_muts,sub_tree$tip.label]) #confirm the coverage
  
  ##Plot the original PVV on the mpboot tre
  hm=create_mutation_heatmap(mut_refs1=PVVs,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=sub_tree)
  
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=1,title = paste(SampleID,":",paste(PVVs,collapse = ",")))
  add_mut_heatmap_PVV(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.3)
  
})
dev.off()

###------------IMPORTING THE LESION SEGREGATION DATA TO LOOK AT AFFECTED BRANCHES------------

out_list=lapply(summary_table_df$Sample_ID,function(sample_ID) {
  binom<-read.delim(paste0(data_dir,"lesion_segregation_analysis/",sample_ID,"_binom.tsv"))
  rl20<-read.delim(paste0(data_dir,"lesion_segregation_analysis/",sample_ID,"_rl20.tsv"))
  wald<-read.delim(paste0(data_dir,"lesion_segregation_analysis/",sample_ID,"_wald.tsv"))
  
  left_join(binom%>%dplyr::select(sample_name,p.value)%>%dplyr::rename("binom"=p.value),
            wald%>%dplyr::select(sample_name,p.value)%>%dplyr::rename("wald"=p.value),by="sample_name")%>%
    left_join(rl20%>%dplyr::select(sample_name,rl20),by="sample_name")%>%
    mutate(Sample_ID=sample_ID,.before=1)
})

dplyr::bind_rows(out_list)%>%
  filter(rl20>=7)

#--------------Compare multiple occurrences of mutations *within* versus *between* individuals i.e. mutation hotspots vs lesion persistence---------

##Import the full mutation datasets across the blood phylogenies
root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")

data_sets=c("MSC_BMT","MF","EM")
lesion_seg_input_dir=paste0(root_dir,"/lesion_segregation/input_data")
all_muts_df=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  all_data_set_mat<-lapply(data_set_samples,function(sample_ID) {
    filtered_muts<-import_muts_file(SampleID = sample_ID,data_set=data_set)
    mutation_mat<-filtered_muts$COMB_mats.tree.build$mat%>%dplyr::mutate(Sample_ID=sample_ID,.before=1)
    return(mutation_mat)
  })%>%dplyr::bind_rows()
  return(all_data_set_mat)
})%>%dplyr::bind_rows()

all_trees=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  all_data_set_trees<-lapply(data_set_samples,function(sample_ID) {
    tree<-import_tree_file(SampleID=sample_ID,data_set = data_set,type = "mpboot")
    return(tree)
  })
  names(all_data_set_trees)<-data_set_samples
  return(all_data_set_trees)
})%>%unlist(all_trees,recursive = F)

#Summarise the number of mutations per individual
summary_df=data.frame(Sample_ID=unique(all_muts_df$Sample_ID))
summary_df$n_muts<-sapply(summary_df$Sample_ID,function(sample_ID) {nrow(all_muts_df%>%filter(Sample_ID==sample_ID))})

#Find all mutations that are present more than once in the full dataset
dup_muts<-all_muts_df$mut_ref[duplicated(all_muts_df$mut_ref)]

#Filter the dataframe to include only these mutation
dup_muts_df<-all_muts_df%>%
  filter(mut_ref%in%dup_muts)%>%
  arrange(mut_ref)

#Annotate if mutations are on private or shared branches
dup_muts_df$branch_type<-sapply(1:nrow(dup_muts_df),function(i) {ifelse(dup_muts_df$node[i] <= length(all_trees[[dup_muts_df$Sample_ID[i]]]$tip.label),"private","shared")})

#Give a sense of how many doublets/ triplets etc there are
table(table(dup_muts_df$mut_ref))

#Create a matrix of these mutations and which samples they are present in
dup_mat<-matrix(0,nrow = length(dup_muts),ncol=length(unique(all_muts_df$Sample_ID)))
rownames(dup_mat)<-dup_muts
colnames(dup_mat)<-unique(all_muts_df$Sample_ID)

for(sample_ID in unique(all_muts_df$Sample_ID)) {
  Sample_dups<-dup_muts_df%>%filter(Sample_ID==sample_ID)%>%pull(mut_ref)
  dup_mat[Sample_dups,sample_ID]<-1
}
head(dup_mat[order(rowSums(dup_mat),decreasing = F),])

#Now do a similar matrix, but one which shows the number of overlaps, relative to the number of opportunities for overlaps (n_mut1 * nmut2)
sample_by_sample_dups<-matrix(NA,nrow=length(unique(all_muts_df$Sample_ID)),ncol=length(unique(all_muts_df$Sample_ID)),dimnames = list(unique(all_muts_df$Sample_ID),(unique(all_muts_df$Sample_ID))))

for(i in 1:nrow(sample_by_sample_dups)) {
  for(j in 1:ncol(sample_by_sample_dups)) {
    Sample1=rownames(sample_by_sample_dups)[i];Sample2=colnames(sample_by_sample_dups)[j]
    nmut_1=summary_df$n_muts[summary_df$Sample_ID==Sample1]; nmut_2=summary_df$n_muts[summary_df$Sample_ID==Sample2]
    sample_by_sample_dups[i,j] <- sum(dup_mat[,i]==1 & dup_mat[,j]==1)/prod(nmut_1,nmut_2)
  }
  #Set the same sample comparisons to 'NA'
  sample_by_sample_dups[i,i]<-NA
}

#Show this as a heatmap
library(pheatmap)
pheatmap(sample_by_sample_dups,cluster_cols=F,cluster_rows = F)

###Look at the mutation signatures of these duplicates
mutations<-tidyr::separate(data.frame(mut_ref=rownames(dup_mat)),col = mut_ref,into = c("chr","pos","ref","mut"),sep = "-",remove = F)%>%mutate(SampleID="BS_duplicates",.before=1)
mutations$pos=as.numeric(mutations$pos)
library("GenomicRanges")
library("Rsamtools")
library("MASS")
genomeFile=lustre_ref="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                   as.numeric(mutations$pos)+1))))
ntcomp = c(T="A",G="C",C="G",A="T")
mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
mutations$trinuc_ref_py = mutations$trinuc_ref
for (j in 1:nrow(mutations)) {
  if (mutations$ref[j] %in% c("A","G")) { # Purine base
    mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
    mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
  }
}

freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec

plot_96profile(mutations,colnames=c("SampleID","chr","pos","ref","mut"),genomeFile = lustre_ref,savefile="recurrent_mutation_mutsig")


#---------------Deeper analysis of the 'phylogeny-confirming' mutations---------------------
library(dplyr)
library(ape)
library(ggplot2)
mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
head(mutations)
PVV_blood_table<-mutations%>%
  filter(Type=="PVV" & cat=="Adult_HSPC" & Class=="PASS")
number_of_confirming_muts.plot<-PVV_blood_table%>%
  filter(Type=="PVV" & cat=="Adult_HSPC" & Class=="PASS")%>%
  ggplot(aes(x=lesion_duration))+
  #scale_x_continuous(limits=c(0,50))+
  geom_histogram(binwidth=1,breaks=seq(0.5,100.5,1),fill="lightblue",col="black",linewidth=0.1)+
  scale_x_continuous(breaks=seq(0,100,10))+
  #geom_histogram(binwidth=1,fill="lightblue",col="black",linewidth=0.1)+
  labs(x="Number of mutations confirming\nthe consensus phylogeny")+
  theme_classic()+
  my_theme

ggsave(paste0(plots_dir,"/confirmatory_muts_histogram.pdf"),width=3, height=2.5)

#Now look at mean coverage and signatures of the confirmatory mutations
data_sets="EM"
confirmatory_mutations_df=lapply(data_sets,function(dataset) {
  cat(dataset,sep="\n")
  dataset_samples=unique(PVV_blood_table%>%filter(data_set==dataset)%>%pull(Sample_ID))
  
  dataset_LRN_muts=lapply(dataset_samples,function(sampleID) {
    cat(sampleID,sep="\n")
    sample_mutations=PVV_blood_table%>%filter(Sample_ID==sampleID)
    info=get_file_paths_and_project(dataset=dataset,Sample_ID = sampleID)
    tree=read.tree(info$tree_file_path)
    load(info$filtered_muts_path)
    details=filtered_muts$COMB_mats.tree.build$mat
    individual_LRN_muts<-lapply(sample_mutations$mut_ref1, function(this_PVV) {
      LN<-sample_mutations%>%filter(mut_ref1==this_PVV)%>%pull(lesion_node)
      LRN<-sample_mutations%>%filter(mut_ref1==this_PVV)%>%pull(lesion_repair_node)
      LN_ancestors<-get_ancestral_nodes(LN,tree$edge)
      LRN_ancestors<-get_ancestral_nodes(LRN,tree$edge)
      LP_nodes<-LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
      cat(LP_nodes,sep = ",")
      all_LP_muts<-details%>%filter(node%in%LP_nodes)%>%mutate(associated_PVV=this_PVV,.before=1)
      all_LP_muts$mean_depth<-rowMeans(filtered_muts$COMB_mats.tree.build$NR[all_LP_muts$mut_ref,,drop=F])
      return(all_LP_muts)
    })%>%dplyr::bind_rows()%>%
      mutate(Sample_ID=sampleID)
    return(individual_LRN_muts)
  })%>%dplyr::bind_rows()%>%
    mutate(data_set=dataset,.before=1)
  return(dataset_LRN_muts)
})%>%dplyr::bind_rows()

write.table(confirmatory_mutations_df,file="confirmatory_mutations.tsv",sep="\t",row.names=F,quote=F)

plot_96profile(confirmatory_mutations_df,colnames=c("Sample_ID","Chrom","Pos","Ref","Alt"),genomeFile = lustre_ref,savefile = "branch_creating_mutations_signature")


mean_cov_df=lapply(data_sets,function(dataset) {
  cat(dataset,sep="\n")
  dataset_samples=unique(PVV_blood_table%>%filter(data_set==dataset)%>%pull(Sample_ID))
  dataset_LRN_muts=lapply(dataset_samples,function(sampleID) {
    cat(sampleID,sep="\n")
    sample_mutations=PVV_blood_table%>%filter(Sample_ID==sampleID)
    info=get_file_paths_and_project(dataset=dataset,Sample_ID = sampleID)
    load(info$filtered_muts_path)
    mean_coverage_by_mut<-rowMeans(filtered_muts$COMB_mats.tree.build$NR)
    return(data.frame(Sample_ID=sampleID,mean_cov=mean_coverage_by_mut))
  })%>%dplyr::bind_rows()%>%
    mutate(data_set=dataset,.before=1)
  return(dataset_LRN_muts)
})%>%dplyr::bind_rows()


coverage_df<-bind_rows(mean_cov_df%>%mutate(Type="Full mutation set"),
          confirmatory_mutations_df%>%dplyr::select(Sample_ID,"mean_cov"=mean_depth)%>%mutate(Type="Branch-creating\nmutations"))

plot<-coverage_df%>%
  ggplot(aes(x=mean_cov))+
  geom_histogram()+
  scale_x_continuous(limits=c(0,25))+
  facet_grid(Type~Sample_ID,scales = "free")+
  theme_classic()+
  my_theme+
  labs(x="Mean coverage (x)")

kruskal.test(x=coverage_df%>%filter(Sample_ID=="KX009_1_01")%>%pull(mean_cov),g=coverage_df%>%filter(Sample_ID=="KX009_1_01")%>%pull(Type))

#coverage_df%>%group_by(Sample_ID,Type)%>%summarise(coverage_means=mean(mean_cov))%>%tidyr::pivot_wider(names_from = "Type",values_from = "coverage_means")
ggsave("branch_creating_mutations_coverage.pdf",plot=plot,width=10,height=4)

#Look for any correlation between age at which lesion was present, and lesion duration---------
#Is there a suggestion of changes is DNA repair efficiency through life?

filtered_mutations_file=paste0(data_dir,"mutations_filtered.tsv")
library(ggpubr)
mutations<-read.delim(filtered_mutations_file)
lesion_duration_vs_age<-mutations%>%
  filter(!(Class=="FAIL"&Type=="MAV") & cat=="Adult_HSPC")%>%
  mutate(lesion_age_midpoint=0.5*(lesion_timing+lesion_repair_timing))%>% #Add in a statistic with the 'midpoint' of age at lesion persistence i.e. halfway along the lesion path
  ggplot(aes(x=lesion_age_midpoint, y=lesion_duration, col = Type))+
  geom_point(alpha=0.5,size=0.5)+
  guides(col=guide_legend(override.aes = list(alpha=1,size=1)))+
  #scale_y_continuous(limits=c(0,200))+
  theme_classic()+
  #stat_cor(aes(x=lesion_repair_timing, y=lesion_duration),method = "pearson", label.x = 50, label.y = 150)+
  my_theme+
  labs(x="Time of lesion persistence\n(Molecular time)",y="Minimum duration of lesion\n(Molecular time)")

ggsave(paste0(plots_dir,"/lesion_duration_vs_age.pdf"),lesion_duration_vs_age,width=4, height=2.5)

#Do a simple linear regression
x<-lm(lesion_duration~lesion_age_midpoint,data=mutations%>%
        filter(!(Class=="FAIL"&Type=="MAV") & cat=="Adult_HSPC")%>%
        mutate(lesion_age_midpoint=0.5*(lesion_timing+lesion_repair_timing)))
summary(x)


###---PLOT A SCALE BARS------

create_color_scale=function(tick_positions,tick_values,color.scale) {
  require(dichromat)
  colfunc = colorRampPalette(color.scale)
  lut=colfunc(101)
  scale = (length(lut)-1)
  plot(c(0,10), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
  axis(side=2,pos=1,las=1,at=tick_positions,labels=as.character(tick_values))
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale
    rect(1,y,2,y+1/scale, col=lut[i], border=NA)
  }
}

cols=c("#FF7F00","#377EB8")
pdf(paste0(plots_dir,"/vaf_scale.pdf"),width = 2,height=3)
create_color_scale(tick_positions = seq(0,1,0.2),tick_values=seq(0,1,0.2),color.scale=c("white",cols[1]))
create_color_scale(tick_positions = seq(0,1,0.2),tick_values=seq(0,1,0.2),color.scale=c("white",cols[2]))
dev.off()





