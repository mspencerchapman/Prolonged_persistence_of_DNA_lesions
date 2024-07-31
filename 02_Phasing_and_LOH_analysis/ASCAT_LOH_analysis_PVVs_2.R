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

#Set file paths for the script - WILL NEED TO UPDATE THESE TO RUN LOCALLY
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
source("/lustre/scratch126/casm/team154pc/ms56/my_programs/Prolonged_persistence_functions.R")
lesion_segregation_output_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/output2/"

setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Set data file paths
data_sets=c("NW","MF","EM","KY","MSC_BMT","MSC_fetal")
out_list=lapply(data_sets,function(data_set) {
  files=list.files(paste0(lesion_segregation_output_dir,data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
  mut_tables_list=lapply(files,read.delim)
  mut_tables_df=Reduce(rbind,mut_tables_list)
  mut_tables_df$data_set=data_set
  return(mut_tables_df)
})
mutations=Reduce(rbind,out_list)
mutations$Ref[mutations$Ref=="TRUE"]<-"T"
mutations$Alt1[mutations$Alt1=="TRUE"]<-"T"
mutations$Alt2[mutations$Alt2=="TRUE"]<-"T"

get_coords_of_mutated_bases=function(mut_ref) {
  Pos<-as.numeric(str_split(mut_ref,pattern = "-",simplify=T)[,2])
  Ref<-str_split(mut_ref,pattern = "-",simplify=T)[,3]
  covered_coords<-Pos:(Pos+nchar(Ref)-1)
  return(covered_coords)
}

#CHECK PHASING
####NOW THE ACTUAL CODE STARTS####

ASCAT_LOH_analysis_PVVs=lapply(data_sets, function(dataset) {
  data_set_PVVs<-mutations%>%filter(Type=="PVV" & data_set==dataset)
  Sample_IDs=unique(data_set_PVVs$Sample_ID)
  
  data_set_out=lapply(Sample_IDs,function(sample) {
    print(paste("Starting analysis for sample",sample))
    data_set_sample_PVVs=data_set_PVVs%>%dplyr::filter(Sample_ID==sample)
    
    #Find the relevant project and tree file
    sample_info=get_file_paths_and_project(dataset,Sample_ID=sample)
    tree_curr=read.tree(sample_info$tree_file_path)
    load(sample_info$filtered_muts_path)
    details=filtered_muts$COMB_mats.tree.build$mat
    matrices=list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR)
    
    if(any(data_set_sample_PVVs$Mut_type1=="MNV")){
      print("Reclassifying the MNVs")
      filtered_muts$COMB_mats.tree.build=reclassify_MNVs(COMB_mats = filtered_muts$COMB_mats.tree.build,genomeFile = genome_file)
      details=filtered_muts$COMB_mats.tree.build$mat
      matrices=list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR) 
    }
    
    sample_out=lapply(1:nrow(data_set_sample_PVVs),function(i) {
      print(i)
      
      #Define the mut parameers
      Chrom=str_split(data_set_sample_PVVs$Chrom_pos[i],pattern = "-",simplify=T)[,1]
      mut=data_set_sample_PVVs$mut_ref1[i]
      print(mut)
      Pos=as.numeric(str_split(data_set_sample_PVVs$Chrom_pos[i],pattern = "-",simplify=T)[,2])
      Ref<-data_set_sample_PVVs$Ref[i]; Alt<-data_set_sample_PVVs$Alt1[i]
      
      #Only able to run if a lesion node has been allocated
      if(is.na(data_set_sample_PVVs$lesion_node[i])) {
        stop(return("No lesion node"))
      }
      
      #Adjust the position/ ref/ alt if mutation is an indel
      if(data_set_sample_PVVs$Mut_type1[i]=="INDEL") {
        print("This is an indel - advise checking result manually")
        Pos=Pos+nchar(data_set_sample_PVVs$Ref[i])-1
        Ref=str_split(Ref,"",simplify=T)[,nchar(Ref)]
        Alt="-"
      }
      
      #Display exactly what will be looked for in the phasing scro[t]
      print(paste(sample,Chrom,Pos))
      print(paste("Using Ref =",Ref,", Alt =",Alt,"at position",Pos))
      
      #Define a random set of samples from the tree used for finding heterozgous SNPs in the .jl phasing script
      set.seed(1)
      Ref_sample_set=paste0(sample(tree_curr$tip.label,min(5,length(tree_curr$tip.label))),collapse=",") #This actually gets ignored if the project is supplied as a look-up df
      pure_subclades=get_pure_subclades(mut1 = mut,lesion_node = data_set_sample_PVVs$lesion_node[i],tree=tree_curr,matrices=matrices)
      print(pure_subclades)
      if(sum(names(pure_subclades)=="pure_positive")<2|!"pure_negative"%in%names(pure_subclades)) {stop(return("Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected"))}
      if(any(pure_subclades=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")) {stop(return("More than one mixed subclade"))}
      
      #Get the mean copy number from ASCAT at the position of the mutation
      #This requires access to output from ASCAT - the functions here have fixed file paths based on the Sanger file structure; will need to edit these to run locally
      subclade_minor_allele_cn=sapply(pure_subclades,function(node,tree=tree_curr) {samples=getTips(tree=tree,node=node);mean_cn=get_mean_ASCAT_minor_allele_cn(samples=samples,Chrom=Chrom,Pos=Pos,project=sample_info$project);return(mean_cn)})
      
      neg_subclade_minor_allele_cn=round(subclade_minor_allele_cn[names(subclade_minor_allele_cn)=="pure_negative"])
      pos_subclade_minor_allele_cn=round(subclade_minor_allele_cn[names(subclade_minor_allele_cn)=="pure_positive"])
      if(all(is.na(neg_subclade_minor_allele_cn))) {
        res<-NA
      } else if(all(neg_subclade_minor_allele_cn[!is.na(neg_subclade_minor_allele_cn)]==0) & !any(pos_subclade_minor_allele_cn==0)){
        res<-"LOH"
      } else if(all(neg_subclade_minor_allele_cn[!is.na(neg_subclade_minor_allele_cn)]==0) & any(pos_subclade_minor_allele_cn==0)){
        res<-"Homozygous in all clades"
      } else if(any(neg_subclade_minor_allele_cn==1) & !any(neg_subclade_minor_allele_cn[!is.na(neg_subclade_minor_allele_cn)]==0)) {
        res<-"No LOH"
      } else {
        res<-"Conflicting"
      }
      print(res)
      return(res)
    })
    return(sample_out)
  })
  return(data_set_out)
})

#Unlist x 2 so that results are just a list by MAV (not grouped by dataset or sample ID)
ASCAT_LOH_analysis_PVVs=unlist(ASCAT_LOH_analysis_PVVs,recursive=F)
ASCAT_LOH_analysis_PVVs=unlist(ASCAT_LOH_analysis_PVVs,recursive=F)
names(ASCAT_LOH_analysis_PVVs)<-mutations%>%filter(Type=="PVV")%>%pull(Chrom_pos)
save(ASCAT_LOH_analysis_PVVs,file = "output2/ASCAT_LOH_analysis_PVVs_all")



