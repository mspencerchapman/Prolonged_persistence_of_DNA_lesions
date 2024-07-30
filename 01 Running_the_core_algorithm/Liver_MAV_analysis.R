library(stringr)
library(ape)
library(seqinr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(phangorn)
library(optparse)
library(parallel)
library("GenomicRanges")
library("Rsamtools")
library("MASS")
options(stringsAsFactors = F)

opt=list(w="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation",
         s="PD48372i",
         t="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/tree_PD48372g.tree",
         f="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/Filtered_muts_PD48372g",
         o="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/output2/SN",
         a=F,
         d=F,
         r=F,
         c=1
)

if(!is.null(opt$w)) {my_working_directory = opt$w} else {my_working_directory<-getwd()}
sampleID=opt$s
tree_file_path = opt$t
filtered_muts_file = opt$f
output_dir=opt$o
include_ancestral_tip=opt$a  #Does your tree include an ancestral tip that you want to keep in the analysis?
remove_duplicates=opt$d
MC_CORES=opt$c
re_run=opt$r
print(paste(MC_CORES,"cores available"))

R_scripts_dir = "/lustre/scratch119/casm/team154pc/ms56/my_functions"
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

filter_output_df_file=paste0(output_dir,"/",sampleID,"_filter_df.tsv")

#Import functions/ files
setwd(my_working_directory)
R_function_files = list.files(R_scripts_dir,pattern=".R",full.names=TRUE)
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Load up the objects, and pull out the details data frame and the NV and NR matrices
filtered_muts_files=list.files("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN",pattern="Filtered_muts")
samples=gsub("Filtered_muts_","",filtered_muts_files)

plot_MAV_simple=function(tree,details,mut1,mut2) {
  mut1_branches=details$node[details$mut_ref==mut1]
  mut2_branches=details$node[details$mut_ref==mut2]
  
  edge.col=rep("lightgray",nrow(tree$edge))
  edge.width=rep(1,nrow(tree$edge))
  edge.col[tree$edge[,2]%in%mut1_branches]<-"blue"
  edge.col[tree$edge[,2]%in%mut2_branches]<-"red"
  edge.col[tree$edge[,2]%in%mut1_branches & tree$edge[,2]%in%mut2_branches]<-"purple"
  edge.width[edge.col!="lightgray"]<-3
  
  
  plot.phylo(tree,direction="downwards",cex=0.3,show.node.label=T,edge.color = edge.col,edge.width = edge.width)
  text(1,1,pos=4,paste0(mut1," / ",mut2))
  
}
all_results=list()

sample_summary=data.frame(sample=samples)
sample_summary$n_samp=sapply(sample_summary$sample,function(sample){
  print(sample)
  tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  return(length(tree$tip.label))
})
sample_summary$pdn=sapply(sample_summary$sample,function(sample){
  print(sample)
  tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  internal_nodes=sapply(1:max(tree$edge),function(node) sum(tree$edge[,1]==node)>1)
  height_over_50=sapply(1:max(tree$edge),function(node) nodeHeights(tree)[tree$edge==node][1])>50
  pdn=sum(internal_nodes&height_over_50)
  return(pdn)
})
sample_summary$n_mut=sapply(sample_summary$sample,function(sample){
  print(sample)
  tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  return(sum(tree$edge.length))
})

load("input_data/SN/x.snv.liver.RData")

for(sample in samples) {
  print(sample)
  tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/tree_",sample,".tree")
  filtered_muts_file=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/Filtered_muts_",sample)
  
  print("Loading the tree and converting to multifurcating structure")
  tree <- di2multi(read.tree(tree_file_path))
  print("Loading the mutations matrix")
  load(filtered_muts_file)
  
  #CORRECTLY CLASSIFY MNVs
  filtered_muts$COMB_mats.tree.build$mat$Mut_type="SNV"
  # filtered_muts$COMB_mats.tree.build$mat$pval=1
  # filtered_muts$COMB_mats.tree.build=reclassify_MNVs(COMB_mats=filtered_muts$COMB_mats.tree.build,region_size=2,genomeFile = genome_file)
  # 
  #Add a "Chrom_pos" column for identifying multi-allelic variants
  details=filtered_muts$COMB_mats.tree.build$mat
  details$Chrom=as.character(details$Chrom)
  details$Ref=as.character(details$Ref)
  details$Alt=as.character(details$Alt)
  NV = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
  NR = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  
  details$Chrom_pos=paste(details$Chrom,details$Pos,sep = "-") #Create a "chrom_pos" column
  
  #Now look for the different mutations at the same site
  print("Assessing mutations for duplicates at the same site")
  MAV_list=get_multi_allelic_variant_list(details,SNV_only=T)
  MAV_list<-lapply(MAV_list,function(vec) {if(length(unique(vec))==1){return(NULL)} else {return(sort(unique(vec)))}})
  MAV_list<-MAV_list[-which(sapply(MAV_list,is.null))]
  if(any(duplicated(MAV_list))){MAV_list<-MAV_list[-which(duplicated(MAV_list))]}
  print(paste0("There are ",length(MAV_list)," mutations at overlapping sites"))
  
  if(any(sapply(MAV_list,length)>2)) {
    print("Evidence of triallelic variants that cannot be handled automatically. These are being removed and should be analysed manually.")
    print(MAV_list[sapply(MAV_list,length)>2])
    MAV_list[sapply(MAV_list,length)>2]<-NULL
  }
  
  if(length(MAV_list)!=0) {
    multi_allelic_list=lapply(MAV_list,function(muts) {
      
      Ref1=details$Ref[details$mut_ref==muts[1]][1]; Ref2=details$Ref[details$mut_ref==muts[2]][1]
      Alt1<-details$Alt[details$mut_ref==muts[1]][1]; Alt2<-details$Alt[details$mut_ref==muts[2]][1]
      Pos1<-as.numeric(details$Pos[details$mut_ref==muts[1]][1]); Pos2<-as.numeric(details$Pos[details$mut_ref==muts[2]][1])
      
      new_refs=establish_ref_and_alt(Ref1,Ref2,Alt1,Alt2,Pos1,Pos2)
      
      clone1=x.snv%>%filter(pos==Pos1 & grepl(sample,id) &alt==Alt1)%>%pull(id)%>%unique()%>%str_split(pattern = "_",simplify=T)
      clone2=x.snv%>%filter(pos==Pos2 & grepl(sample,id) &alt==Alt2)%>%pull(id)%>%unique()%>%str_split(pattern = "_",simplify=T)
      
      clone1<-clone1[,2]
      clone2<-clone2[,2]
      
      nodes=1:max(tree$edge)
      names(nodes)<-c(tree$tip.label,tree$node.label)
      
      
      df=data.frame(Chrom_pos=details$Chrom_pos[details$mut_ref==muts[1]][1],
                    mut_ref1=muts[1],
                    mut_ref2=muts[2],
                    Ref=new_refs$Ref,
                    Alt1=new_refs$Alt1,
                    Alt2=new_refs$Alt2,
                    Mut_type1=details$Mut_type[details$mut_ref==muts[1]][1],
                    Mut_type2=details$Mut_type[details$mut_ref==muts[2]][1],
                    Node1=nodes[clone1],
                    Node2=nodes[clone2],
                    Type="MAV",
                    Class=NA,
                    Sample_ID=sample,
                    lesion_node=NA,
                    lesion_timing=NA,
                    lesion_repair_node=NA,
                    lesion_repair_timing=NA,
                    lesion_duration=NA,
                    no_of_cell_divisions=NA,
                    base_order=NA,
                    MAV_is_PVV=NA,
                    colonies_per_negative_subclade=NA,
                    depth_per_negative_subclade=NA
      )
      return(df)
    })
    multi_allelic_df=dplyr::bind_rows(multi_allelic_list)
    
    for(i in 1:nrow(multi_allelic_df)){
      details$node[details$mut_ref==multi_allelic_df$mut_ref1[i]]<-multi_allelic_df$Node1[i]
      details$node[details$mut_ref==multi_allelic_df$mut_ref2[i]]<-multi_allelic_df$Node2[i]
    }
    
    
    multi_allelic_df<-multi_allelic_df[!duplicated(multi_allelic_df$Chrom_pos),]
    ROOT=tree$edge[1,1]
    multi_allelic_df$Class=sapply(1:nrow(multi_allelic_df),function(i){
      node1=multi_allelic_df$Node1[i]
      node2=multi_allelic_df$Node2[i]
      parent1=get_ancestor_node(multi_allelic_df$Node1[i],tree = tree)
      parent2=get_ancestor_node(multi_allelic_df$Node2[i],tree = tree)
      
      if(parent1%in%c(parent2,get_all_node_children(parent2,tree))|parent2%in%c(parent1,get_all_node_children(parent1,tree))){
        ROOT=tree$edge[1,1]
        if(parent1==ROOT&!parent2%in%c(node1,get_all_node_children(node1,tree))|parent2==ROOT&!parent1%in%c(node2,get_all_node_children(node2,tree))){
          return("FAIL")
        } else {
          return("PASS")
        }
      } else {
        return("FAIL")
      }
    })
    
    all_results[[sample]]<-multi_allelic_df
  }
}

load("~/Downloads/Supplementary Code/x.structure.RData")
load("~/Downloads/Supplementary Code/all.trees.RData")

load("x.structure.RData")
load("all.trees.RData")

x.structure<-Map(tree=all.trees.updated,structure_df=x.structure,f=function(tree,structure_df) {
  internal_node_labels<-tree$node.label
  internal_node_numbers<-(length(tree$tip.label)+1):max(tree$edge)
  tip_node_numbers<-1:length(tree$tip.label)
  tip_labels<-tree$tip.label
  convert_node_label_vec<-c(tip_node_numbers,internal_node_numbers)
  names(convert_node_label_vec)<-c(tip_labels,internal_node_labels)
  structure_df$node1<-sapply(structure_df$From,function(label) convert_node_label_vec[paste("Cl",label,sep=".")])
  structure_df$node2<-sapply(structure_df$To,function(label) convert_node_label_vec[paste("Cl",label,sep=".")])
  return(structure_df)
})

all_results<-lapply(samples,function(sampleID){
  sample_file=paste0(output_dir,"/",sampleID,"_mut_table.tsv")
  if(!file.exists(sample_file)) {
    return(NULL)
  } else {
    df<-read.delim(paste0(output_dir,"/",sampleID,"_mut_table.tsv"),stringsAsFactors = F)
    return(df)
  }
})
names(all_results)<-samples

dataset="SN"
SN_project_info<-read.csv("/lustre/scratch119/realdata/mdt1/team154/ms56/lesion_segregation/input_data/SN/Samples_project_ref_SN.csv")
phasing_output_dir=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/phasing_output/",dataset)
all_phasing_results<-Map(results=all_results[names(all.trees.updated)],structure_df=x.structure,tree=all.trees.updated,f=function(results,structure_df,tree) {
  if(is.null(results)) {
    return(NULL)
  } else {
    results$Chrom=str_split(results$Chrom_pos,pattern="-",simplify=T)[,1]
    results$Pos=as.numeric(str_split(results$Chrom_pos,pattern="-",simplify=T)[,2])
    all_samples=unique(unlist(str_split(structure_df$microdissection_id,pattern=", ")))
    phasing_results=vector(mode="list",length=nrow(results))
    for(i in 1:nrow(results)){
      MAV1_samples=unlist(str_split(structure_df$microdissection_id[which(structure_df$node2==results$Node1[i])][[1]],pattern=", "))
      MAV2_samples=unlist(str_split(structure_df$microdissection_id[which(structure_df$node2==results$Node2[i])][[1]],pattern=", "))
      
      set.seed(1)
      Ref_sample_set=sample(all_samples[!all_samples%in%c(MAV1_samples,MAV2_samples)],5)
      phasing_list1=get_phasing_list(samples=MAV1_samples,Chrom=results$Chrom[i],Pos=results$Pos[i],tree=tree,project=SN_project_info,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set,use_tree = F,verbose = F)
      phasing_list2=get_phasing_list(samples=MAV2_samples,Chrom=results$Chrom[i],Pos=results$Pos[i],tree=tree,project=SN_project_info,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set,use_tree = F,verbose = F)
      
      base_counts_list1=get_base_counts_list(samples=MAV1_samples,Chrom=results$Chrom[i],Pos=results$Pos[i],tree=tree,project=SN_project_info,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set,use_tree = F,verbose = F)
      base_counts_list2=get_base_counts_list(samples=MAV2_samples,Chrom=results$Chrom[i],Pos=results$Pos[i],tree=tree,project=SN_project_info,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set,use_tree = F,verbose = F)
      
      
      phasing_info_by_subclade=lapply(list(phasing_list1,phasing_list2),extract_phasing_info,Ref=results$Ref[i], Alt=c(results$Alt1[i],results$Alt2[i]))
      phasing_results[[i]]<-list(phasing_info_by_subclade=phasing_info_by_subclade,base_counts=list(base_counts_list1,base_counts_list2))
    }
    return(phasing_results)
  }
})

#Save the results to the output2 directory is the same format as all the other datasets
output_dir="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/output2/SN"
mapply(FUN=function(df,sampleID) {
  write.table(df,file = paste0(output_dir,"/",sampleID,"_mut_table.tsv"),quote = FALSE,sep = "\t",row.names = FALSE)
},df=all_results,sampleID=names(all_results))

