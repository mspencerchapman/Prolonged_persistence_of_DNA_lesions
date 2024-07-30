library(dplyr)
library(ape)


#Source functions needed for the script
my_working_directory<-getwd()
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

source("/lustre/scratch126/casm/team154pc/ms56/my_programs/Prolonged_persistence_functions.R")
lesion_seg_input_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data"


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

#-------------------------------------------------------------------------------------------
#----------------- EMILY SAMPLES-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

###FIND AND RUN SCITE FOR ANY SAMPLES THAT DON'T ALREADY HAVE IT
data_set="EM"
data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))

SCITE_comp_list=vector(mode="list",length = length(data_set_samples))
names(SCITE_comp_list)<-data_set_samples

for(SampleID in data_set_samples) {
  cat(SampleID,sep="\n"); shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
  my_working_directory=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered/")
  output_dir=my_working_directory; setwd(my_working_directory)
  
  filtered_muts_file=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered/annotated_mut_set_",SampleID,"_standard_rho01")
  tree_file_path=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered/tree_",SampleID,"_standard_rho01.tree")
  
  #Set SCITE paths
  scite_input_file_path=paste0("scite_input_",SampleID)
  scite_path="/lustre/scratch126/casm/team154pc/ms56/programs/SCITE/scite"
  scite_output_tree_path=paste0(scite_input_file_path,"_standard_rho01_ml0.newick")
  scite_output_tree_path_option2=paste0(scite_input_file_path,"_ml0.newick")
  scite_allocated_muts_tree_path=paste0("scite_tree_allocated_muts_",SampleID,".tree")
  
  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  tree_mod<-di2multi(tree)
  tree_mod$edge.length[tree_mod$edge.length<10]<-10
  
  #2. SCITE
  #SCITE has a binary input format where 0 is absent, 1 is present, and 3 is "missing data".
  if(!file.exists(scite_output_tree_path) & !file.exists(scite_output_tree_path_option2)){
    cat("No scite tree found - setting off the SCITE algorithm",sep="\n")
    
    cat("Reading in the filtered muts file",sep="\n")
    if(!file.exists(filtered_muts_file)) {next}
    load(filtered_muts_file)
    gt<-filtered_muts$Genotype_shared_bin; NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR); details<-filtered_muts$COMB_mats.tree.build$mat
    
    #Write file in required format for SCITE
    gt_rows=apply(gt,1,paste,collapse=" ") #Collapse into strings
    gt_rows=sapply(gt_rows, function(x) gsub("0.5","3", x)) #Replace the 0.5's with 3's as per the SCITE input
    writeLines(gt_rows,scite_input_file_path)
    
    #Write the SCITE command
    scite_command=paste0(scite_path," -i ",scite_input_file_path," -n ",nrow(gt)," -m ",ncol(gt)," -r 1 -l 1000000 -fd 0.001 -ad 0.005 -transpose")
    #system(paste0("bsub -o $PWD/log.%J -e $PWD/err.%J -q basement -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J scite_",SampleID," ",scite_command)) #SCITE takes a long time: days to weeks
    
    SCITE_comp_list[[SampleID]]<-data.frame(RF_diff=NA,Quartet_diff=NA)
    
  } else if(!file.exists(scite_allocated_muts_tree_path)){
    cat("SCITE tree is found, but without the allocated mutations. Running this algorithm.",sep="\n")
    
    if(!file.exists(scite_output_tree_path) & file.exists(scite_output_tree_path_option2)) {
      scite_output_tree_path<-scite_output_tree_path_option2
    }
    
    tree_scite=read.tree(text=paste0(readLines(scite_output_tree_path),";"))
    tree_scite$edge.length=rep(1,nrow(tree_scite$edge))
    tree_scite$tip.label<-sapply(tree_scite$tip.label,function(x) colnames(gt)[as.numeric(x)]) #Map back the sample names from the numbers
    
    cat("Reading in the filtered muts file",sep="\n")
    if(!file.exists(filtered_muts_file)) {next}
    load(filtered_muts_file)
    gt<-filtered_muts$Genotype_shared_bin; NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR); details<-filtered_muts$COMB_mats.tree.build$mat
    
    #Drop tips that aren't in the original tree
    tree_scite<-drop.tip(tree_scite,colnames(gt)[!colnames(gt)%in%colnames(NV)])
    
    #Assign mutations back to the tree
    df = reconstruct_genotype_summary(tree_scite) #Define df (data frame) for treeshape
    res = assign_to_tree(NV[,df$samples],NR[,df$samples],df,error_rate = c(rep(0.01, ncol(NR))))
    tree_scite$edge.length <- res$df$df$edge_length #Assign edge lengths from the res object
    tree_scite<-di2multi(tree_scite)
    write.tree(tree_scite,file=scite_allocated_muts_tree_path)
    
  } else if(file.exists(scite_allocated_muts_tree_path)){
    cat("SCITE tree is found with allocated mutations. Comparing this tree to the MPBoot phylogeny.",sep="\n")
    tree_scite=read.tree(scite_allocated_muts_tree_path)
    tree_scite_mod<-tree_scite
    tree_scite_mod$edge.length[tree_scite_mod$edge.length<10]<-10
    
    pdf(paste0(output_dir,"/MPBoot_vs_SCITE_comparisons_",SampleID,"_",".pdf"),width = 15,height=6)
    comparePhylo_and_plot(di2multi(tree),tree_scite,names=c("MPBoot phylogeny","SCITE phylogeny"))
    dev.off()
    
    pdf(paste0(output_dir,"/MPBoot_vs_SCITE_comparisons_mod_",SampleID,"_",".pdf"),width = 15,height=6)
    comparePhylo_and_plot(di2multi(tree_mod),tree_scite_mod,names=c("MPBoot phylogeny","SCITE phylogeny"))
    dev.off()
    
    #Calculate and print the similarity metrics for the scite tree
    tree_scite_splits_comp<-Quartet::SplitStatus(tree_scite,di2multi(tree))
    tree_scite_quartets_comp<-Quartet::QuartetStatus(tree_scite,di2multi(tree))
    RF_diff<-round(Quartet::SimilarityMetrics(tree_scite_splits_comp)[,"SymmetricDifference"],digits=5)
    Quartet_diff<-round(Quartet::SimilarityMetrics(tree_scite_quartets_comp)[,"SymmetricDifference"],digits = 5)
    print(paste("The Robinson-Foulds similarity is",RF_diff))
    print(paste("The Quartet Similarity is",Quartet_diff))
    
    SCITE_comp_list[[SampleID]]<-data.frame(RF_diff=RF_diff,Quartet_diff=Quartet_diff)
  }
}
temp=dplyr::bind_rows(SCITE_comp_list)%>%mutate(SampleID=data_set_samples,.before=1)%>%dplyr::rename("SCITE_RF_diff"=RF_diff,"SCITE_Quartet_diff"=Quartet_diff)

###FIND AND RUN IQTree FOR ANY SAMPLES THAT DON'T ALREADY HAVE IT
IQTREE_comp_list=vector(mode="list",length = length(data_set_samples))
names(IQTREE_comp_list)<-data_set_samples

for(SampleID in data_set_samples) {
  cat(SampleID,sep="\n"); shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
  my_working_directory=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered/")
  output_dir=my_working_directory; setwd(my_working_directory)
  filtered_muts_file=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered/annotated_mut_set_",SampleID,"_standard_rho01")
  tree_file_path=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered/tree_",SampleID,"_standard_rho01.tree")
  
  #Set IQtree paths
  iqtree_path="/lustre/scratch126/casm/team154pc/ms56/programs/IQ-TREE/build/iqtree"
  iqtree_input_file_path="iqtree_fasta.fa"
  iqtree_output_tree_path=paste0(iqtree_input_file_path,".treefile")
  iqtree_allocated_muts_tree_path=paste0("iqtree_tree_allocated_muts_",SampleID,".tree")
  
  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  tree_mod<-di2multi(tree)
  tree_mod$edge.length[tree_mod$edge.length<10]<-10

  if(!file.exists(iqtree_output_tree_path)){
    cat("No iqtree tree found - setting off the IQtree algorithm",sep="\n")
    
    cat("Reading in the filtered muts file",sep="\n")
    if(!file.exists(filtered_muts_file)) {next}
    load(filtered_muts_file)
    gt<-filtered_muts$Genotype_shared_bin; NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR); details<-filtered_muts$COMB_mats.tree.build$mat
    
    #Write file in required format for IQtree
    seqinr::write.fasta(lapply(filtered_muts$dna_strings,function(x) gsub("W","0",gsub("V","1",x))),names=names(filtered_muts$dna_strings),iqtree_input_file_path)
    #system(paste0(iqtree_path," -s ",iqtree_input_file_path," -m JC2 -bb 1000 -czb -redo")) #-czb option allows collapse to polytomies
    IQTREE_comp_list[[SampleID]]<-data.frame(RF_diff=NA,Quartet_diff=NA)
    
  } else if(!file.exists(iqtree_allocated_muts_tree_path)){
    cat("IQtree tree is found, but without the allocated mutations. Running this algorithm.",sep="\n")
    
    cat("Reading in the filtered muts file",sep="\n")
    if(!file.exists(filtered_muts_file)) {next}
    load(filtered_muts_file)
    gt<-filtered_muts$Genotype_shared_bin; NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR); details<-filtered_muts$COMB_mats.tree.build$mat
    
    #Read in the output tree
    tree_iq=read.tree(iqtree_output_tree_path)
    tree_iq <- drop.tip(tree_iq,"Ancestral")
    tree_iq<-di2multi(tree_iq)
    tree_iq$edge.length = rep(1, nrow(tree_iq$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
    
    #Drop tips that aren't in the original tree
    tree_iq<-drop.tip(tree_iq,colnames(gt)[!colnames(gt)%in%colnames(NV)])
    
    #Assign mutations back to the tree
    df = reconstruct_genotype_summary(tree_iq) #Define df (data frame) for treeshape
    res = assign_to_tree(NV[,df$samples],NR[,df$samples],df,error_rate = c(rep(0.01, ncol(NR))))
    tree_iq$edge.length <- res$df$df$edge_length #Assign edge lengths from the res object
    tree_iq<-di2multi(tree_iq)
    write.tree(tree_iq,file=iqtree_allocated_muts_tree_path)
    
  } else if(file.exists(iqtree_allocated_muts_tree_path)){
    cat("IQtree tree is found with allocated mutations. Comparing this tree to the MPBoot phylogeny.",sep="\n")
    tree_iq=read.tree(iqtree_allocated_muts_tree_path)
    tree_iq<-drop.tip(tree_iq,tip = tree_iq$tip.label[!tree_iq$tip.label%in%tree$tip.label]) #In PX001, there are some additional samples in the bootstrapped trees
    
    #Versions of the tree where short branches are lengthened to a 10 to allow visualization of differences
    tree_iq_mod<-tree_iq; tree_iq_mod$edge.length[tree_iq_mod$edge.length<10]<-10
    
    #Calculate and print the similarity metrics for the iqtree tree
    tree_iq_splits_comp<-Quartet::SplitStatus(tree_iq,di2multi(tree))
    RF_diff<-round(Quartet::SimilarityMetrics(tree_iq_splits_comp)[,"SymmetricDifference"],digits=5)
    
    #Quartet is unable to run the quartet status for very large trees (includes KX004)
    if(length(tree_iq$tip.label)>477) {
      Quartet_diff<-NA
    } else {
      tree_iq_quartets_comp<-Quartet::QuartetStatus(tree_iq,di2multi(tree))
      Quartet_diff<-round(Quartet::SimilarityMetrics(tree_iq_quartets_comp)[,"SymmetricDifference"],digits = 5)
    }
   
    print(paste("The Robinson-Foulds similarity is",RF_diff))
    print(paste("The Quartet Similarity is",Quartet_diff))
    
    IQTREE_comp_list[[SampleID]]<-data.frame(RF_diff=RF_diff,Quartet_diff=Quartet_diff)
  }
}
temp2=dplyr::bind_rows(IQTREE_comp_list)%>%mutate(SampleID=data_set_samples,.before=1)%>%dplyr::rename("IQ_RF_diff"=RF_diff,"IQ_Quartet_diff"=Quartet_diff)

#-------------------------------------------------------------------------------------------
#----------------- MARGA SAMPLES-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

###FIND AND RUN SCITE FOR ANY MARGA SAMPLES THAT DON'T ALREADY HAVE IT
data_set="MF"
data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))

SCITE_comp_list_MF=vector(mode="list",length = length(data_set_samples))
names(SCITE_comp_list_MF)<-data_set_samples

for(SampleID in data_set_samples) {
  cat(SampleID,sep="\n"); shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
  my_working_directory=paste0("/lustre/scratch126/casm/team154pc/ms56/Marga_benchmarking/",SampleID)
  output_dir=my_working_directory; setwd(my_working_directory)
  
  all_filtered_muts_files=list.files("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/MF",pattern="annotated_mut_set",full.names = T)
  all_tree_files=list.files("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/MF",pattern="tree_",full.names = T)
  filtered_muts_file=grep(shortID,all_filtered_muts_files,value=T)
  tree_file_path=grep(shortID,all_tree_files,value=T)
  
  #Set SCITE paths
  scite_input_file_path=paste0("scite_input_",SampleID)
  scite_path="/lustre/scratch126/casm/team154pc/ms56/programs/SCITE/scite"
  scite_output_tree_path=paste0(scite_input_file_path,"_ml0.newick")
  scite_allocated_muts_tree_path=paste0("scite_tree_allocated_muts_",SampleID,".tree")
  
  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  cat("Reading in the filtered muts file",sep="\n")
  if(!file.exists(filtered_muts_file)) {next}
  load(filtered_muts_file)
  gt<-filtered_muts$Genotype_shared_bin; NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR); details<-filtered_muts$COMB_mats.tree.build$mat
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  tree_mod<-di2multi(tree)
  tree_mod$edge.length[tree_mod$edge.length<10]<-10
  
  #2. SCITE
  #SCITE has a binary input format where 0 is absent, 1 is present, and 3 is "missing data".
  gt=filtered_muts$Genotype_shared_bin
  if(!file.exists(scite_output_tree_path)){
    cat("No scite tree found - setting off the SCITE algorithm",sep="\n")
    
    #Write file in required format for SCITE
    gt_rows=apply(gt,1,paste,collapse=" ") #Collapse into strings
    gt_rows=sapply(gt_rows, function(x) gsub("0.5","3", x)) #Replace the 0.5's with 3's as per the SCITE input
    writeLines(gt_rows,scite_input_file_path)
    
    #Write the SCITE command
    scite_command=paste0(scite_path," -i ",scite_input_file_path," -n ",nrow(gt)," -m ",ncol(gt)," -r 1 -l 1000000 -fd 0.001 -ad 0.005 -transpose")
    #system(paste0("bsub -o $PWD/log.%J -e $PWD/err.%J -q basement -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J scite_",SampleID," ",scite_command)) #SCITE takes a long time: days to weeks
    
  } else if(!file.exists(scite_allocated_muts_tree_path)){
    cat("SCITE tree is found, but without the allocated mutations. Running this algorithm.",sep="\n")
    
    tree_scite=read.tree(text=paste0(readLines(scite_output_tree_path),";"))
    tree_scite$edge.length=rep(1,nrow(tree_scite$edge))
    tree_scite$tip.label<-sapply(tree_scite$tip.label,function(x) colnames(gt)[as.numeric(x)]) #Map back the sample names from the numbers
    
    #Drop tips that aren't in the original tree
    tree_scite<-drop.tip(tree_scite,colnames(gt)[!colnames(gt)%in%colnames(NV)])
    
    #Assign mutations back to the tree
    df = reconstruct_genotype_summary(tree_scite) #Define df (data frame) for treeshape
    res = assign_to_tree(NV[,df$samples],NR[,df$samples],df,error_rate = c(rep(0.01, ncol(NR))))
    tree_scite$edge.length <- res$df$df$edge_length #Assign edge lengths from the res object
    tree_scite<-di2multi(tree_scite)
    write.tree(tree_scite,file=scite_allocated_muts_tree_path)
    
  } else if(file.exists(scite_allocated_muts_tree_path)){
    cat("SCITE tree is found with allocated mutations. Comparing this tree to the MPBoot phylogeny.",sep="\n")
    tree_scite=read.tree(scite_allocated_muts_tree_path)
    tree_scite_mod<-tree_scite
    tree_scite_mod$edge.length[tree_scite_mod$edge.length<10]<-10
    
    pdf(paste0(output_dir,"/MPBoot_vs_SCITE_comparisons_",SampleID,"_",".pdf"),width = 15,height=6)
    comparePhylo_and_plot(di2multi(tree),tree_scite,names=c("MPBoot phylogeny","SCITE phylogeny"))
    dev.off()
    
    pdf(paste0(output_dir,"/MPBoot_vs_SCITE_comparisons_mod_",SampleID,"_",".pdf"),width = 15,height=6)
    comparePhylo_and_plot(di2multi(tree_mod),tree_scite_mod,names=c("MPBoot phylogeny","SCITE phylogeny"))
    dev.off()
    
    #Calculate and print the similarity metrics for the scite tree
    tree_scite_splits_comp<-Quartet::SplitStatus(tree_scite,di2multi(tree))
    tree_scite_quartets_comp<-Quartet::QuartetStatus(tree_scite,di2multi(tree))
    RF_diff<-round(Quartet::SimilarityMetrics(tree_scite_splits_comp)[,"SymmetricDifference"],digits=5)
    Quartet_diff<-round(Quartet::SimilarityMetrics(tree_scite_quartets_comp)[,"SymmetricDifference"],digits = 5)
    print(paste("The Robinson-Foulds similarity is",RF_diff))
    print(paste("The Quartet Similarity is",Quartet_diff))
    
    SCITE_comp_list_MF[[SampleID]]<-data.frame(RF_diff=RF_diff,Quartet_diff=Quartet_diff)
  }
}

###FIND AND RUN IQTree FOR ANY SAMPLES THAT DON'T ALREADY HAVE IT
IQTREE_comp_list_MF=vector(mode="list",length = length(data_set_samples))
names(IQTREE_comp_list_MF)<-data_set_samples

for(SampleID in data_set_samples) {
  cat(SampleID,sep="\n"); shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
  my_working_directory=paste0("/lustre/scratch126/casm/team154pc/ms56/Marga_benchmarking/",SampleID)
  output_dir=my_working_directory; setwd(my_working_directory)
  
  all_filtered_muts_files=list.files("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/MF",pattern="annotated_mut_set",full.names = T)
  all_tree_files=list.files("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/MF",pattern="tree_",full.names = T)
  filtered_muts_file=grep(shortID,all_filtered_muts_files,value=T)
  tree_file_path=grep(shortID,all_tree_files,value=T)
  
  #Set IQtree paths
  iqtree_path="/lustre/scratch126/casm/team154pc/ms56/programs/IQ-TREE/build/iqtree"
  iqtree_input_file_path="iqtree_fasta.fa"
  iqtree_output_tree_path=paste0(iqtree_input_file_path,".treefile")
  iqtree_allocated_muts_tree_path=paste0("iqtree_tree_allocated_muts_",SampleID,".tree")
  
  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  cat("Reading in the filtered muts file",sep="\n")
  if(!file.exists(filtered_muts_file)) {next}
  load(filtered_muts_file)
  gt<-filtered_muts$Genotype_shared_bin; NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR); details<-filtered_muts$COMB_mats.tree.build$mat
  tree<-drop.tip(read.tree(tree_file_path),"Ancestral")
  tree_mod<-di2multi(tree)
  tree_mod$edge.length[tree_mod$edge.length<10]<-10
  
  gt=filtered_muts$Genotype_shared_bin
  if(!file.exists(iqtree_output_tree_path)){
    cat("No iqtree tree found - setting off the IQtree algorithm",sep="\n")
    
    #Write file in required format for IQtree
    seqinr::write.fasta(lapply(filtered_muts$dna_strings,function(x) gsub("W","0",gsub("V","1",x))),names=names(filtered_muts$dna_strings),iqtree_input_file_path)
    #system(paste0(iqtree_path," -s ",iqtree_input_file_path," -m JC2 -bb 1000 -czb -redo")) #-czb option allows collapse to polytomies
    
  } else if(!file.exists(iqtree_allocated_muts_tree_path)){
    cat("IQtree tree is found, but without the allocated mutations. Running this algorithm.",sep="\n")
    
    #Read in the output tree
    tree_iq=read.tree(iqtree_output_tree_path)
    tree_iq <- drop.tip(tree_iq,"Ancestral")
    tree_iq<-di2multi(tree_iq)
    tree_iq$edge.length = rep(1, nrow(tree_iq$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
    
    #Drop tips that aren't in the original tree
    tree_iq<-drop.tip(tree_iq,colnames(gt)[!colnames(gt)%in%colnames(NV)])
    
    #Assign mutations back to the tree
    df = reconstruct_genotype_summary(tree_iq) #Define df (data frame) for treeshape
    res = assign_to_tree(NV[,df$samples],NR[,df$samples],df,error_rate = c(rep(0.01, ncol(NR))))
    tree_iq$edge.length <- res$df$df$edge_length #Assign edge lengths from the res object
    tree_iq<-di2multi(tree_iq)
    write.tree(tree_iq,file=iqtree_allocated_muts_tree_path)
    
  } else if(file.exists(iqtree_allocated_muts_tree_path)){
    cat("IQtree tree is found with allocated mutations. Comparing this tree to the MPBoot phylogeny.",sep="\n")
    tree_iq=read.tree(iqtree_allocated_muts_tree_path)
    #Versions of the tree where short branches are lengthened to a 10 to allow visualization of differences
    tree_iq_mod<-tree_iq
    tree_iq_mod$edge.length[tree_iq_mod$edge.length<10]<-10
    
    
    #Calculate and print the similarity metrics for the iqtree tree
    tree_iq_splits_comp<-Quartet::SplitStatus(tree_iq,di2multi(tree))
    tree_iq_quartets_comp<-Quartet::QuartetStatus(tree_iq,di2multi(tree))
    RF_diff<-round(Quartet::SimilarityMetrics(tree_iq_splits_comp)[,"SymmetricDifference"],digits=5)
    Quartet_diff<-round(Quartet::SimilarityMetrics(tree_iq_quartets_comp)[,"SymmetricDifference"],digits = 5)
    print(paste("The Robinson-Foulds similarity is",RF_diff))
    print(paste("The Quartet Similarity is",Quartet_diff))
    
    IQTREE_comp_list_MF[[SampleID]]<-data.frame(RF_diff=RF_diff,Quartet_diff=Quartet_diff)
  }
}


# Import only----------------------

SCITE_trees<-lapply(c("EM","MF"),function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  
  data_set_trees<-lapply(data_set_samples,function(SampleID) {
    if(data_set=="EM") {
      cat(SampleID,sep="\n"); shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
      output_dir=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered")
    } else if (data_set=="MF") {
      output_dir=paste0("/lustre/scratch126/casm/team154pc/ms56/Marga_benchmarking/",SampleID)
    }
    scite_allocated_muts_tree_path=paste0(output_dir,"/scite_tree_allocated_muts_",SampleID,".tree")
    if(file.exists(scite_allocated_muts_tree_path)){
      cat("SCITE tree is found with allocated mutations. Comparing this tree to the MPBoot phylogeny.",sep="\n")
      tree_scite=read.tree(scite_allocated_muts_tree_path)
      return(tree_scite)
    } else {
      return(NULL)
    }
  })
  names(data_set_trees)<-data_set_samples
  return(data_set_trees)
})

IQTREE_trees<-lapply(c("EM","MF"),function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  
  data_set_trees<-lapply(data_set_samples,function(SampleID) {
    if(data_set=="EM") {
      cat(SampleID,sep="\n"); shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
      output_dir=paste0("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/",shortID,"/vaf_filtered")
    } else if (data_set=="MF") {
      cat(SampleID,sep="\n")
      output_dir=paste0("/lustre/scratch126/casm/team154pc/ms56/Marga_benchmarking/",SampleID)
    }
    iqtree_allocated_muts_tree_path=paste0(output_dir,"/iqtree_tree_allocated_muts_",SampleID,".tree")
    if(file.exists(iqtree_allocated_muts_tree_path)){
      cat("IQTREE tree is found with allocated mutations. Comparing this tree to the MPBoot phylogeny.",sep="\n")
      tree_iq=read.tree(iqtree_allocated_muts_tree_path)
      return(tree_iq)
    } else {
      return(NULL)
    }
  })
  names(data_set_trees)<-data_set_samples
  return(data_set_trees)
})


