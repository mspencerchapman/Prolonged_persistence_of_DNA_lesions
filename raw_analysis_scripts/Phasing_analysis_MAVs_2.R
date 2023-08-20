#!/software/R-3.6.1/bin/Rscript
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

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Set data file paths
data_sets=c("MSC_chemo","NW","MF","EM","KY","PR","MSC_BMT","MSC_fetal")#,"SN")
out_list=lapply(data_sets,function(data_set) {
  files=list.files(paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/output2/",data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
  mut_tables_list=lapply(files,read.delim)
  mut_tables_df=Reduce(rbind,mut_tables_list)
  mut_tables_df$data_set=data_set
  return(mut_tables_df)
})
mutations=Reduce(rbind,out_list)
mutations$Ref[mutations$Ref=="TRUE"]<-"T"
mutations$Alt1[mutations$Alt1=="TRUE"]<-"T"
mutations$Alt2[mutations$Alt2=="TRUE"]<-"T"

genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

get_coords_of_mutated_bases=function(mut_ref) {
  Pos<-as.numeric(str_split(mut_ref,pattern = "-",simplify=T)[,2])
  Ref<-str_split(mut_ref,pattern = "-",simplify=T)[,3]
  covered_coords<-Pos:(Pos+nchar(Ref)-1)
  return(covered_coords)
}

#CHECK PHASING
####NOW THE ACTUAL CODE STARTS####

Phasing_results_MAVs=lapply(data_sets, function(dataset) {
  phasing_output_dir=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/phasing_output/",dataset)
  data_set_MAVs<-mutations%>%filter(Type=="MAV" & data_set==dataset)
  Sample_IDs=unique(data_set_MAVs$Sample_ID)
  data_set_out=lapply(Sample_IDs,function(sample) {
    print(paste("Starting analysis for sample",sample))
    data_set_sample_MAVs=data_set_MAVs%>%dplyr::filter(Sample_ID==sample)
    
    #Find the relevant project and tree file
    sample_info=get_file_paths_and_project(dataset,Sample_ID=sample)
    tree_curr=read.tree(sample_info$tree_file_path)
    load(sample_info$filtered_muts_path)
    details=filtered_muts$COMB_mats.tree.build$mat
    matrices=list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR)
    
    if(any(data_set_sample_MAVs$Mut_type1=="MNV")|any(data_set_sample_MAVs$Mut_type2=="MNV")){
      print("Reclassifying the MNVs")
      filtered_muts$COMB_mats.tree.build=reclassify_MNVs(COMB_mats = filtered_muts$COMB_mats.tree.build,genomeFile = genome_file)
      details=filtered_muts$COMB_mats.tree.build$mat
      matrices=list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR) 
    }
    
    sample_out=lapply(1:nrow(data_set_sample_MAVs),function(i) {
      print(i)
      Chrom=str_split(data_set_sample_MAVs$Chrom_pos[i],pattern = "-",simplify=T)[,1]
      mut1=data_set_sample_MAVs$mut_ref1[i]
      mut2=data_set_sample_MAVs$mut_ref2[i]
      
      print(paste("Mut 1 is:",mut1,", Mut 2 is:", mut2))
      
      Ref=data_set_sample_MAVs$Ref[i]
      Alt1<-data_set_sample_MAVs$Alt1[i]
      Alt2<-data_set_sample_MAVs$Alt2[i]
      
      #Decide which position to phase in all samples
      if(data_set_sample_MAVs$Mut_type1[i]=="SNV"&data_set_sample_MAVs$Mut_type2[i]=="SNV") {
        Pos=as.numeric(str_split(data_set_sample_MAVs$Chrom_pos[i],pattern = "-",simplify=T)[,2])
      } else {
        covered_1=get_coords_of_mutated_bases(mut1)
        covered_2=get_coords_of_mutated_bases(mut2)
        covered_both=intersect(covered_1,covered_2)
        Pos<-max(covered_both) #Choose the max of the intersect as INDELs are annotated including a preceding base that actually remains unchanged
        shift<-Pos-as.numeric(str_split(data_set_sample_MAVs$Chrom_pos[i],pattern = "-",simplify=T)[,2])+1
        Alt1<-unlist(strsplit(Alt1,""))[shift]
        Alt2<-unlist(strsplit(Alt2,""))[shift]
        Ref<-unlist(strsplit(Ref,""))[shift]
        if(data_set_sample_MAVs$Mut_type1[i]=="INDEL") {
          Alt1<-"-"
        }
        if(data_set_sample_MAVs$Mut_type2[i]=="INDEL") {
          Alt2<-"-"
        }
      }
      
      #Display exactly what will be looked for in the phasing script
      print(paste(sample,Chrom,Pos))
      print(paste("Using Ref =",Ref,", Alt1 =",Alt1,", Alt2 =",Alt2,"at position",Pos))
      
      if(dataset=="SN") {
        mut1_pos_samples=colnames(matrices$NV[mut1,matrices$NV[mut1,]>=1])
        mut2_pos_samples=colnames(matrices$NV[mut2,matrices$NV[mut2,]>=1])
        set.seed(1)
        Ref_sample_set=sample(size=5,x=colnames(matrices$NV))
        sample_sets=list(mut1_pos_samples,mut2_pos_samples)
        phasing_lists=lapply(sample_sets,function(samples) {get_phasing_list(samples=samples,Chrom=Chrom,Pos=as.integer(Pos),project=sample_info$project,tree=tree_curr,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set,use_tree=F)})
        phasing_summaries=Map(f=function(phasing_list,Alt){extract_phasing_info(phasing_list,Ref=Ref,Alt=Alt)},phasing_list=phasing_lists,Alt=c(Alt1,Alt2))
        res=assess_phasing_non_clonal(phasing_summaries=phasing_summaries,sample_sets=sample_sets,Chrom=Chrom,Pos=as.integer(Pos),project=sample_info$project,tree=tree_curr,output_dir=phasing_output_dir,ref_sample_set = Ref_sample_set,use_tree=F)
        print(res)
        phasing_info_by_subclade=phasing_summaries; positive_subclade_res=res; negative_subclade_res=NA
      } else {
        #Define a random set of samples from the tree used for finding heterozgous SNPs in the .jl phasing script
        set.seed(1)
        Ref_sample_set=paste0(sample(tree_curr$tip.label,min(5,length(tree_curr$tip.label))),collapse=",") #This actually gets ignored if the project is supplied as a look-up df
        
        #Find which samples should have each MAV from the tree & the node in the mutations table
        if(is.na(data_set_sample_MAVs$lesion_node[i])|data_set_sample_MAVs$Class[i]=="simple") {
          
          pure_subclade_phasing_list=lapply(c(data_set_sample_MAVs$Node1[i],data_set_sample_MAVs$Node2[i]),function(node,tree=tree_curr) {samples=getTips(tree=tree,node=node);phasing_list=get_phasing_list(samples=samples,Chrom=Chrom,Pos=as.integer(Pos),project=sample_info$project,tree=tree,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set);return(phasing_list)})
          phasing_info_by_subclade=lapply(pure_subclade_phasing_list,extract_phasing_info,Ref=Ref, Alt=c(Alt1,Alt2))
          print(phasing_info_by_subclade)
          
          #Compare the alt phasing of all the positive subclades
          positive_subclade_res=check_matching_phasing(phasing_info_by_subclade[[1]],phasing_info_by_subclade[[2]])
          
          if(any(grepl("suggested",unlist(positive_subclade_res)))) {
            #Get aggregated base counts from apparent heterozygous SNPs in each of the subclades
            aggregated_subclade_base_counts_list=get_clade_base_counts(c(data_set_sample_MAVs$Node1[i],data_set_sample_MAVs$Node2[i]),tree=tree_curr,Chrom=Chrom,Pos=Pos,project=sample_info$project,ref_sample_set=Ref_sample_set,phasing_output_dir=phasing_output_dir)
            #Confirm those that are truly heterozygous in the positive subclades - this may include the mutation itself
            het_positions<-return_heterozygous_SNPs(base_counts_list=aggregated_subclade_base_counts_list)
            het_positions<-het_positions[het_positions!=Pos] #Exclude the position of the actual mutation
            positive_subclade_res=check_matching_phasing(phasing_info_by_subclade[[1]],phasing_info_by_subclade[[2]],het_positions=het_positions)
          }
          print(positive_subclade_res)
          negative_subclade_res=NA
        } else {
          #Here is the PVV assessment code
          pure_subclades=get_pure_subclades(mut1 = mut1, mut2 = mut2,lesion_node = data_set_sample_MAVs$lesion_node[i],tree=tree_curr,matrices=matrices)
          print(pure_subclades)
          if(!"pure_mut1"%in%names(pure_subclades)|!"pure_mut2"%in%names(pure_subclades)) {stop(return("Need both MAVs to be within the lesion node: one or both not detected"))}
          if(any(pure_subclades=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")) {stop(return("More than one mixed subclade"))}
          pure_subclade_phasing_list=lapply(pure_subclades,function(node,tree=tree_curr) {samples=getTips(tree=tree,node=node);phasing_list=get_phasing_list(samples=samples,Chrom=Chrom,Pos=Pos,project=sample_info$project,tree=tree,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set);return(phasing_list)})
          phasing_info_by_subclade=lapply(pure_subclade_phasing_list,extract_phasing_info,Ref=Ref, Alt=c(Alt1,Alt2))
          
          print(phasing_info_by_subclade)
          
          #Compare the alt phasing of all the positive subclades
          positive_subclades=which(names(phasing_info_by_subclade)%in%c("pure_mut1","pure_mut2"))
          positive_subclade_res=lapply(2:length(positive_subclades),function(j) {
            result=check_matching_phasing(phasing_info_by_subclade[[positive_subclades[1]]],phasing_info_by_subclade[[positive_subclades[j]]])
            names(result)<-paste0(pure_subclades[positive_subclades][1]," (",names(pure_subclades[positive_subclades])[1],") vs ",pure_subclades[positive_subclades][j]," (",names(pure_subclades[positive_subclades])[j],")")
            print(result)
            return(result)
          })
          
          if(any(grepl("suggested",unlist(positive_subclade_res)))) {
            #Get aggregated base counts from apparent heterozygous SNPs in each of the subclades
            aggregated_subclade_base_counts_list=get_clade_base_counts(pure_subclades,tree=tree_curr,Chrom=Chrom,Pos=Pos,project=sample_info$project,ref_sample_set=Ref_sample_set,phasing_output_dir=phasing_output_dir)
            #Confirm those that are truly heterozygous in the positive subclades - this may include the mutation itself
            het_positions<-return_heterozygous_SNPs(base_counts_list=aggregated_subclade_base_counts_list[positive_subclades])
            het_positions<-het_positions[het_positions!=Pos] #Exclude the position of the actual mutation
            
            #Are any of the SNPs used to confirm the phasing in the set of "confirmed heterozygous SNPs"? Update the results accordingly
            positive_subclade_res=lapply(2:length(positive_subclades),function(j) {
              result=check_matching_phasing(phasing_info_by_subclade[[positive_subclades[1]]],phasing_info_by_subclade[[positive_subclades[j]]],het_positions=het_positions)
              names(result)<-paste0(pure_subclades[positive_subclades][1]," (",names(pure_subclades[positive_subclades])[1],") vs ",pure_subclades[positive_subclades][j]," (",names(pure_subclades[positive_subclades])[j],")")
              print(result)
              return(result)
            })
          }
          
          if(any(names(phasing_info_by_subclade)=="pure_negative")) {
            negative_subclade_res=lapply(which(names(phasing_info_by_subclade)=="pure_negative"),function(k) {
              result<-check_for_both_alleles_confirming_ref(phasing_info_by_subclade[[k]])
              print(result)
              return(result)
            }) 
          } else {
            negative_subclade_res=NA
          }
        }
      }
      return(list(phasing_info_by_subclade=phasing_info_by_subclade,positive_subclade_res=positive_subclade_res,negative_subclade_res=negative_subclade_res))
    })
    return(sample_out)
  })
  return(data_set_out)
})

#Unlist x 2 so that results are just a list by MAV (not grouped by dataset or sample ID)
Phasing_results_MAVs=unlist(Phasing_results_MAVs,recursive=F)
Phasing_results_MAVs=unlist(Phasing_results_MAVs,recursive=F)

names(Phasing_results_MAVs)<-mutations%>%filter(Type=="MAV")%>%pull(Chrom_pos)

save(Phasing_results_MAVs,file = "output2/Phasing_results_MAVs_all")
