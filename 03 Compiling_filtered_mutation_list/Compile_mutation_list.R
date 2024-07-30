#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr")
bioconductor_packages=c("GenomicRanges","IRanges","MutationalPatterns","MASS","Rsamtools")

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

palette2=c("#E30613","#1D71B8")
palette6=c("#72e5ef", "#074d65", "#3d99ce", "#c257d3", "#d1add5", "#61356e")
palette8=c("#2EBAED" ,"#000000" ,"#DE1C14", "#E98C7B", "#D4D2D2" ,"#ADCC54" ,"#F0D0CE","blue")
my_theme<-theme_classic()+
  theme(text = element_text(family="Helvetica"),
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


#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir="~/R_work/Prolonged_persistence_of_DNA_lesions/"
plots_dir=paste0(root_dir,"plots/")
data_dir=paste0(root_dir,"Data/")
output_dir=paste0(root_dir,"/output/")
source(paste0(root_dir,"/Data/Prolonged_persistence_functions.R"))
ref_table=read.csv(paste0(root_dir,"/Data/metadata/Individual_ref.csv"))
genome_file=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/reference_files/genome.fa","/nfs/cancer_ref02/human/GRCh37d5/genome.fa")##Set to the location of GRCh37 genome file

options(stringsAsFactors = FALSE)


Phasing_MAV_file_path=paste0(data_dir,"phasing_results/Phasing_results_MAVs_all")
Phasing_PVV_file_path=paste0(data_dir,"phasing_results/Phasing_results_PVVs_all")
ASCAT_PVV_file_path=paste0(data_dir,"ASCAT_LOH_analysis/ASCAT_LOH_analysis_PVVs_all")
ASCAT_MAV_file_path=paste0(data_dir,"ASCAT_LOH_analysis/ASCAT_LOH_analysis_MAVs_all")
SN_phasing_file_path=paste0(data_dir,"phasing_results/Phasing_results_MAVs_SN.csv")

lesion_seg_input_dir=paste0(data_dir,"/input_data")
lesion_seg_output_dir=output_dir

mutations_file=paste0(data_dir,"mutations.tsv")
summary_df_file=paste0(data_dir,"summary_df.tsv")

if(file.exists(mutations_file)&file.exists(summary_df_file)) {
  
  mutations<-read.delim(mutations_file)
  summary_table_df<-read.delim(summary_df_file)
  sample_ref=read.csv(paste0(root_dir,"/Data/metadata/Individual_ref.csv"))
  
} else {
  #Set data file paths
  data_sets=c("NW","MF","EM","KY","MSC_BMT","MSC_fetal","SN")
  out_list=lapply(data_sets,function(data_set) {
    print(data_set)
    files=list.files(paste0(output_dir,"/",data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
    mut_tables_list=lapply(files,readr::read_delim,delim="\t",col_types=c("cccccccciiccciiiiiiccccc"))
    mut_tables_list=lapply(mut_tables_list,function(df) {
      df$colonies_per_negative_subclade<-as.character(df$colonies_per_negative_subclade)
      df$depth_per_negative_subclade<-as.character(df$depth_per_negative_subclade)
      df$Ref<-as.character(df$Ref);df$Alt1<-as.character(df$Alt1);df$Alt2<-as.character(df$Alt2)
      return(df)})
    mut_tables_df=dplyr::bind_rows(mut_tables_list)
    mut_tables_df$data_set=data_set
    return(mut_tables_df)
  })
  mutations=dplyr::bind_rows(out_list)
  
  #ADD THE TRINUCLEOTIDE REFERENCE
  mutations$Chrom=str_split(mutations$Chrom_pos,pattern="-",simplify=TRUE)[,1]
  mutations$Pos=as.numeric(str_split(mutations$Chrom_pos,pattern="-",simplify=TRUE)[,2])
  mutations$trinuc_ref = as.vector(scanFa(genome_file, GRanges(mutations$Chrom, IRanges(mutations$Pos-1, mutations$Pos+1))))
  
  #Convert to the format of having the pyrmidine reference base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$Sub1 = paste(mutations$Ref,mutations$Alt1,sep=">")
  mutations$Sub2 = paste(mutations$Ref,mutations$Alt2,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$Ref[j] %in% c("A","G")) { # Purine base
      mutations$Sub1[j] = paste(ntcomp[mutations$Ref[j]],ntcomp[mutations$Alt1[j]],sep=">")
      mutations$Sub2[j] = paste(ntcomp[mutations$Ref[j]],ntcomp[mutations$Alt2[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  #Create the 96 profile (needed later)
  pre_vec=rep(c("A","C","G","T"),each=4,times=6)
  mid_vec=rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16)
  post_vec=rep(c("A","C","G","T"),times=24)
  full_vec=paste0(pre_vec,"[",mid_vec,"]",post_vec)
  mutations$mut_profile_1=paste0(substr(mutations$trinuc_ref_py,1,1),"[",mutations$Sub1,"]",substr(mutations$trinuc_ref_py,3,3))
  mutations$mut_profile_2=paste0(substr(mutations$trinuc_ref_py,1,1),"[",mutations$Sub2,"]",substr(mutations$trinuc_ref_py,3,3))
  
  #Add individual-level metadata
  sample_ref=read.csv(paste0(data_dir,"metadata/Individual_ref.csv"))
  mutations$cat=sapply(mutations$Sample_ID,function(sample) return(sample_ref$Category[sample_ref$Sample_ID==sample]))
  mutations<-mutations%>%filter(data_set!="PR")
  
  write.table(mutations,file=mutations_file,quote=F,sep="\t",row.names = F)
  
  #Summary table
  summary_table_list=lapply(data_sets,function(data_set) {
    print(data_set)
    data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
    data_set_df_list=lapply(data_set_samples,function(sample){
      sample_info=get_file_paths_and_project(data_set,Sample_ID=sample,input_data_dir = lesion_seg_input_dir)
      tree=read.tree(sample_info$tree_file_path)
      df=data.frame(Sample_ID=sample,
                    n_sample=length(tree$tip.label),
                    n_mutations=sum(tree$edge.length),
                    sharedness=calculate_sharedness_stat_2(tree),
                    development_nodes=count_internal_nodes(tree,cut_off=50,type="below"),
                    post_development_nodes=count_internal_nodes(tree,cut_off=50))
      return(df)
    })
    
    data_set_df=dplyr::bind_rows(data_set_df_list)
    data_set_df$data_set=data_set
    return(data_set_df)
  })
  summary_table_df=dplyr::bind_rows(summary_table_list)%>%filter(data_set!="PR")
  summary_table_df$cat=sapply(summary_table_df$Sample_ID,function(sample) return(sample_ref$Category[sample_ref$Sample_ID==sample]))
  summary_table_df$smoking_status=sapply(summary_table_df$Sample_ID,function(sample) return(sample_ref$Smoking_status[sample_ref$Sample_ID==sample]))
  summary_table_df$smoking_status<-ifelse(summary_table_df$smoking_status=="","unknown",summary_table_df$smoking_status)
  
  write.table(summary_table_df,file=summary_df_file,quote=F,sep="\t",row.names = F)
}

#========================================#
# Generate filtered MAV list ####
#========================================#

MAV_mutations=mutations%>%filter(Type=="MAV")
MAV_mutations<-MAV_mutations%>%mutate(Original_Class=Class,.before = "Class") #Record the 'original outcome' from the calling algorithm
dim(MAV_mutations)

## Add info on the number of negative subclades between the two MAV clades----
# As seen later, if >1 negative subclades, these are more likely to be independently-acquired mutations
MAV_mutations$n_neg=sapply(MAV_mutations$colonies_per_negative_subclade,function(x) length(unlist(strsplit(x[!is.na(x)],","))))

## Load the phasing information ----
#This is generated by the script in 02 Phasing_and_LOH_analysis/Phasing_analysis_MAVs_2.R
load(Phasing_MAV_file_path)
Phasing_results_MAVs<-Phasing_results_MAVs[-which(duplicated(names(Phasing_results_MAVs)))] #remove any duplicates

## Extract the basic phasing info to add to the master dataframe ----
MAV_phasing_summary=sapply(Phasing_results_MAVs,extract_MAV_pos_clade_phasing_summary)
names(MAV_phasing_summary)<-str_split(names(MAV_phasing_summary),pattern="\\.",simplify = T)[,1]
MAV_phasing_summary<-MAV_phasing_summary[!duplicated(names(MAV_phasing_summary))]
MAV_mutations<-left_join(MAV_mutations,data.frame(Chrom_pos=names(MAV_phasing_summary),phasing_summary=MAV_phasing_summary),by="Chrom_pos")

## Add the phasing information from the liver samples ----
# Note that these are analyzed separately as they are non-clonal samples so have a different framework
SN_phasing_summary<-read.csv(SN_phasing_file_path,stringsAsFactors=F)
MAV_mutations$phasing_summary<-sapply(1:nrow(MAV_mutations),function(i) {
  if(MAV_mutations$Sample_ID[i]%in%SN_phasing_summary$exp_ID) {
    return(SN_phasing_summary%>%filter(Chrom_pos==MAV_mutations$Chrom_pos[i])%>%pull(res))
  } else {
    return(MAV_mutations$phasing_summary[i])
  }
})

## Create simplified phasing metric, with only 3 categories ----
MAV_mutations<-MAV_mutations%>%
  mutate(phasing_summary=ifelse(phasing_summary%in%c("Same phasing confirmed","Non-matching phasing confirmed"),phasing_summary,"Unable to confirm phasing"))%>%
  mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

#Now extract info on the negative subclades (where there are any)
MAV_both_alleles_summary=sapply(Phasing_results_MAVs,extract_PVV_neg_clade_phasing_summary)
MAV_both_alleles_summary<-MAV_both_alleles_summary[!duplicated(names(MAV_both_alleles_summary))]
MAV_neg_clades_phasing_df=data.frame(Chrom_pos=names(MAV_both_alleles_summary),result=sapply(MAV_both_alleles_summary,function(x) return(x[1])))
MAV_neg_clades_phasing_df$basic_result=sapply(MAV_neg_clades_phasing_df$result,function(result) {
  if(is.na(result)) {
    return(NA)
  } else if(result%in%c("Alt allele reads present","Alt allele reads present in at least one subclade","Both alleles confirmed with reference allele","Both alleles confirmed with reference allele in at least one subclade")){
    return("No LOH")
  } else if(grepl("with maximum",result)){
    n_reads<-as.numeric(sapply(str_extract_all(result,"\\d"),paste,collapse=""))
    if(n_reads<=4) {return("Unable to confirm")} else if(n_reads>4) {return("LOH")}
  } else if(result=="Positive clades have non-matching phasing"){return(NA)} else {return("Unable to confirm")}
})
MAV_mutations=left_join(MAV_mutations,MAV_neg_clades_phasing_df[,c("Chrom_pos","result","basic_result")],by="Chrom_pos") #Join this to the mutations dataframe
MAV_mutations%>%filter(Class!="FAIL" & (is.na(basic_result)|basic_result!="LOH"))%>%pull(phasing_summary)%>%table()

## Review the 7 non-matching phasing mutations ----
MAV_mutations%>%
  mutate(phasing_summary=factor(phasing_summary))%>%filter(phasing_summary%in%c("Non-matching phasing confirmed")&Class!="FAIL")

## Re-annotate the 7 non-matching phasing samples as 'FAIL' ----
MAV_mutations<-MAV_mutations%>%mutate(Class=ifelse(phasing_summary=="Non-matching phasing confirmed","FAIL",Class))

## Given the high rate of non-matching phasing if > 1 negative subclade, re-annotate all of these as "FAIL"----
MAV_mutations%>%group_by(phasing_summary)%>%summarise(n_exclude=sum(n_neg>1 & Class=="removed"))
MAV_mutations<-MAV_mutations%>%mutate(Class=ifelse(Class=="removed"&n_neg>1,"FAIL",Class))


dim(MAV_mutations) #No mutations should have been lost from the list, but the annotations/ filter results are added

#========================================#
# Generate filtered PVV list ####
#========================================#

PVV_mutations<-mutations%>%filter(Type=="PVV")
PVV_mutations<-PVV_mutations%>%mutate(Original_Class=Class,.before = "Class") #Record the 'original outcome' from the calling algorithm
dim(PVV_mutations)

## Load the ASCAT data to confirm no LOH ----
#This is generated by the script in 02 Phasing_and_LOH_analysis/ASCAT_LOH_analysis_PVVs_2.R
load(ASCAT_PVV_file_path)
ASCAT_PVV_df<-data.frame(Chrom_pos=names(ASCAT_LOH_analysis_PVVs),ASCAT_result=unlist(ASCAT_LOH_analysis_PVVs))
PVV_mutations<-left_join(PVV_mutations,ASCAT_PVV_df,by="Chrom_pos")

## Extract the phasing info regarding the positive subclades ----
#This is generated by the script in 02 Phasing_and_LOH_analysis/Phasing_analysis_PVVs_2.R
load(Phasing_PVV_file_path)
Phasing_results_PVVs<-Phasing_results_PVVs[-which(duplicated(names(Phasing_results_PVVs)))] #remove any duplicates

## Add phasing information on the positive subclades ----
PVV_phasing_summary=sapply(Phasing_results_PVVs,extract_PVV_pos_clade_phasing_summary); names(PVV_phasing_summary)<-str_split(names(PVV_phasing_summary),pattern="\\.",simplify=T)[,1]
PVV_mutations<-left_join(PVV_mutations,data.frame(Chrom_pos=names(PVV_phasing_summary),PVV_pos_clade_phasing=PVV_phasing_summary),by="Chrom_pos")
#Simplify the output to 3 categories: same-phasing, non-matching phasing, or unable to confirm
PVV_mutations$PVV_pos_clade_phasing<-sapply(PVV_mutations$PVV_pos_clade_phasing,function(x) {
  if(grepl("Non-matching phasing confirmed",x)) {
    res<-"Non-matching phasing confirmed"
  } else if(grepl("Same phasing confirmed",x)){
    res<-"Same phasing confirmed"
  } else if(is.na(x)) {
    res<-"Unable to confirm phasing" #A few variants could not be run in the phasing analysis
  } else {
    res<-"Unable to confirm phasing"
  }
  return(res)
})
PVV_mutations$PVV_pos_clade_phasing<-factor(PVV_mutations$PVV_pos_clade_phasing,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing"))
PVV_mutations<-mutate(PVV_mutations,Class=ifelse(is.na(lesion_node)|ASCAT_result=="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected","FAIL","PASS"))

## Add in information on the negative subclades ----
# This annotates if reads from both alleles are confirmed present with the reference base, i.e. excludes local LOH
PVV_both_alleles_summary=sapply(Phasing_results_PVVs,extract_PVV_neg_clade_phasing_summary)
Neg_clades_phasing_df=data.frame(Chrom_pos=names(Phasing_results_PVVs),result=sapply(PVV_both_alleles_summary,function(x) return(x[1])))

## Simplify the output----
# Inclues simple categories of: "No loss of mutant allele", "Unable to confirm", "Likely LOH" or "Positive clade have non-matching phasing"
Neg_clades_phasing_df$basic_result=sapply(Neg_clades_phasing_df$result,function(result) {
  if(result%in%c("Alt allele reads present","Both alleles confirmed with reference allele","Both alleles confirmed with reference allele in at least one subclade")){
    return("No LOH")
  } else if(grepl("with maximum",result)){
    n_reads<-as.numeric(sapply(str_extract_all(result,"\\d"),paste,collapse=""))
    if(n_reads<=4) {return("Unable to confirm")} else if(n_reads>4) {return("LOH")}
  } else if(result=="Positive clades have non-matching phasing"){return(NA)} else {return("Unable to confirm")}
})
PVV_mutations=left_join(PVV_mutations,Neg_clades_phasing_df[,c("Chrom_pos","result","basic_result")],by="Chrom_pos") #Join this to the mutations dataframe
PVV_mutations$basic_result[is.na(PVV_mutations$basic_result)]<-"Unable to confirm"
PVV_mutations%>%filter(Class!="FAIL"&ASCAT_result!="LOH"&basic_result!="LOH")%>%pull(PVV_pos_clade_phasing)%>%table()

## Introduce a more stringent depth threshold of ≥13X ----
#for the minimum depth of the negative subclade - needs to be at least ≥13 (not 10 as previously)
PVV_mutations$max_neg_clade_depth=sapply(PVV_mutations$depth_per_negative_subclade,function(x) max(as.numeric(unlist(strsplit(x,split=",")))))
PVV_mutations<-PVV_mutations%>%mutate(Class=ifelse(max_neg_clade_depth<13|is.na(max_neg_clade_depth)|is.na(lesion_node),"FAIL",Class))

## Remove likely artefactual mutations found by manual review ----
#Several mutations were found to be likely artefactual by manual review during the peer-review process, these are manually annotated here
artefactual_muts=c("3-57906539-C-T",
                   "3-182038851-AT-A",
                   "1-32201632-C-T",
                   "10-63986636-C-T",
                   "14-82666754-G-A",
                   "6-104547705-C-A")

#Note that some "non-genuine" PVVs (by a variety of mechanisms) will still remain as not all individual mutations are fully assessable
PVV_mutations<-PVV_mutations%>%mutate(Class=ifelse(ASCAT_result!="LOH"& basic_result!="LOH" & PVV_pos_clade_phasing!="Non-matching phasing confirmed"&no_of_cell_divisions<=6 & !mut_ref1%in%artefactual_muts,Class,"FAIL"))

#========================================#
# Combine the annotated mutation lists and save ####
#========================================#
mutations_filt<-dplyr::bind_rows(MAV_mutations,PVV_mutations)
View(mutations_filt)

#========================================#
# Save this into the data directory ####
#========================================#
mutations_filt_file=paste0(data_dir,"mutations_filtered.tsv")
write.table(mutations_filt,file=mutations_filt_file,quote = F,sep="\t",row.names=F)


#========================================#
# Get into format for Table S2 for the manuscript ####
#========================================#
individual_metadata=read_csv(paste0(root_dir,"Tables/Individual_metadata.csv"))

mutations_filt%>%
  left_join(individual_metadata,by=c("Sample_ID"="Phylogeny ID","cat"="Category"))%>%
  mutate(phasing_summary=ifelse(is.na(phasing_summary),as.character(PVV_pos_clade_phasing),as.character(phasing_summary)))%>%
  dplyr::select("Phylogeny ID"=Sample_ID,Age, Sex,`Study ID`,"Category"=cat, `Sample type`,Type, Class,Chrom,Pos,"Mutation reference 1"=mut_ref1,
                "Mutation reference 2"=mut_ref2, "Reference allele"=Ref,"Alternate allele 1"=Alt1,"Alternate allele 2"=Alt2,
                "Mutation type 1"=Mut_type1,"Mutation type 2"=Mut_type2,"Assigned node for mutation 1"=Node1, "Assigned node for mutation 2"=Node2,
                "Lesion node"=lesion_node,"Lesion timing"=lesion_timing,"Lesion repair node"=lesion_repair_node, "Lesion repair timing"=lesion_repair_timing,
                "Lesion duration"=lesion_duration,"No of cell divisions with lesion"=no_of_cell_divisions,"Base order"=base_order,"Is MAV also a PVV?"=MAV_is_PVV,
                "Number of colonies in each wild type clade"=colonies_per_negative_subclade,
                "Depth in each wild type clade"=depth_per_negative_subclade,"Trinucleotide reference"=trinuc_ref,
                "Trinucleotide reference (pyrimidine based)"=trinuc_ref_py,"Substitution 1 (pyrimidine based)"=Sub1,
                "Substitution 2 (pyrimidine based)"=Sub2,"Mutation profile 1"=mut_profile_1,"Mutation profile 2"=mut_profile_2,
                "Number of wild type clades between mutation clades"=n_neg,"Phasing summary"=phasing_summary,
                "Reads-based LOH assessment of wild-type clades"=basic_result,"ASCAT-based LOH assessment of wild-type clades"=ASCAT_result,
                "Maximum depth of wild0type subclades (PVVs only)"=max_neg_clade_depth)%>%
  write_csv(file=paste0(root_dir,"Tables/Mutation_table_updated.csv"))

#========================================#
# Get numbers passing filters at each point ####
#========================================#

#Original number of PVVs
nrow(PVV_mutations)

#Number of PVVs remaining after removimng those not in a consistent orientation for persistent lesion
PVV_mutations%>%
  filter(!is.na(lesion_node)&ASCAT_result!="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected")%>%
  nrow()

PVV_mutations%>%
  filter(!is.na(lesion_node)&ASCAT_result!="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected")%>%
  filter(no_of_cell_divisions<=6)%>%
  nrow()

PVV_mutations%>%
  filter(!is.na(lesion_node)&ASCAT_result!="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected")%>%
  filter(no_of_cell_divisions<=6)%>%
  filter(PVV_pos_clade_phasing!="Non-matching phasing confirmed")%>%
  nrow()

PVV_mutations%>%
  filter(!is.na(lesion_node)&ASCAT_result!="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected")%>%
  filter(no_of_cell_divisions<=6)%>%
  filter(PVV_pos_clade_phasing!="Non-matching phasing confirmed")%>%
  filter(ASCAT_result!="LOH"& basic_result!="LOH")%>%
  nrow()

PVV_mutations%>%
  filter(!is.na(lesion_node)&ASCAT_result!="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected")%>%
  filter(no_of_cell_divisions<=6)%>%
  filter(PVV_pos_clade_phasing!="Non-matching phasing confirmed")%>%
  filter(ASCAT_result!="LOH"& basic_result!="LOH")%>%
  filter(max_neg_clade_depth>=13)%>%
  nrow()

PVV_mutations%>%
  filter(!is.na(lesion_node)&ASCAT_result!="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected")%>%
  filter(no_of_cell_divisions<=6)%>%
  filter(PVV_pos_clade_phasing!="Non-matching phasing confirmed")%>%
  filter(ASCAT_result!="LOH"& basic_result!="LOH")%>%
  filter(max_neg_clade_depth>=13)%>%
  filter(!mut_ref1%in%artefactual_muts)%>%
  nrow()
