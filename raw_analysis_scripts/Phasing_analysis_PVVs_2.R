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
data_sets=c("MSC_chemo","NW","MF","EM","KY","PR","MSC_BMT","MSC_fetal")
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

Phasing_results_PVVs=lapply(data_sets, function(dataset) {
  phasing_output_dir=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/phasing_output/",dataset)
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
      
      #Dfine the mut parameers
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
      
      #Here is the PVV assessment code
      pure_subclades=get_pure_subclades(mut1 = mut,lesion_node = data_set_sample_PVVs$lesion_node[i],tree=tree_curr,matrices=matrices)
      print(pure_subclades)
      if(sum(names(pure_subclades)=="pure_positive")<2|!"pure_negative"%in%names(pure_subclades)) {stop(return("Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected"))}
      if(any(pure_subclades=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")) {stop(return("More than one mixed subclade"))}
      pure_subclade_phasing_list=lapply(pure_subclades,function(node,tree=tree_curr) {samples=getTips(tree=tree,node=node);phasing_list=get_phasing_list(samples=samples,Chrom=Chrom,Pos=Pos,project=sample_info$project,tree=tree,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set);return(phasing_list)})
      phasing_info_by_subclade=lapply(pure_subclade_phasing_list,extract_phasing_info,Ref=Ref, Alt=Alt)
      print(phasing_info_by_subclade)
      
      #Compare the alt/ ref phasing of all the positive subclades
      positive_subclades=which(names(phasing_info_by_subclade)=="pure_positive")
      positive_subclade_res=lapply(2:length(positive_subclades),function(j) {
        result=check_matching_phasing(phasing_info_by_subclade[[positive_subclades[1]]],phasing_info_by_subclade[[positive_subclades[j]]])
        names(result)<-paste0(pure_subclades[positive_subclades][1]," (",names(pure_subclades[positive_subclades])[1],") vs ",pure_subclades[positive_subclades][j]," (",names(pure_subclades[positive_subclades])[j],")")
        print(result)
        return(result)
      })
      
      #If the phasing result is only "suggested" need to confirm the heterozygosity of the SNPs used to confirm phasing
      if(any(grepl("suggested",unlist(positive_subclade_res)))) {
        #Get aggregated base counts from apparent heterozygous SNPs in each of the subclades
        aggregated_subclade_base_counts_list=get_clade_base_counts(pure_subclades,tree=tree_curr,Chrom=Chrom,Pos=Pos,project=sample_info$project,ref_sample_set=Ref_sample_set,phasing_output_dir=phasing_output_dir)
        #Confirm those that are truly heterozygous in the positive subclades - this may include the mutation itself
        het_positions<-return_heterozygous_SNPs(base_counts_list=aggregated_subclade_base_counts_list[positive_subclades])
        het_positions<-het_positions[het_positions!=Pos] #Exclude the position of the actual mutation
        
        if(length(het_positions)==0){
          positive_subclade_res<-lapply(positive_subclade_res,function(list) return("Unable to confirm phasing"))
        } else {
          #Are any of the SNPs used to confirm the phasing in the set of "confirmed heterozygous SNPs"? Update the results accordingly
          positive_subclade_res=lapply(2:length(positive_subclades),function(j) {
            result=check_matching_phasing(phasing_info_by_subclade[[positive_subclades[1]]],phasing_info_by_subclade[[positive_subclades[j]]],het_positions=het_positions)
            names(result)<-paste0(pure_subclades[positive_subclades][1]," (",names(pure_subclades[positive_subclades])[1],") vs ",pure_subclades[positive_subclades][j]," (",names(pure_subclades[positive_subclades])[j],")")
            print(result)
            return(result)
          })
        }
      }
      
      #Now check phasing of the negative subclades
      if(any(names(phasing_info_by_subclade)=="pure_negative")) {
        negative_subclade_res=lapply(which(names(phasing_info_by_subclade)=="pure_negative"),function(k) {
          result<-check_for_both_alleles_confirming_ref(phasing_info_by_subclade[[k]])
          print(result)
          return(result)
        }) 
      } else {
        negative_subclade_res=NA
      }
      
      return(list(phasing_info_by_subclade=phasing_info_by_subclade,positive_subclade_res=positive_subclade_res,negative_subclade_res=negative_subclade_res))
    })
    return(sample_out)
  })
  return(data_set_out)
})

#Unlist x 2 so that results are just a list by MAV (not grouped by dataset or sample ID)
Phasing_results_PVVs=unlist(Phasing_results_PVVs,recursive=F)
Phasing_results_PVVs=unlist(Phasing_results_PVVs,recursive=F)

names(Phasing_results_PVVs)<-mutations%>%filter(Type=="PVV")%>%pull(Chrom_pos)

save(Phasing_results_PVVs,file="output2/Phasing_results_PVVs_all")



