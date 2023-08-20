#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("ggplot2","dplyr","tidyr","gridExtra","ggrepel","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
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
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Set the ggplot2 theme and palettes for plotting
#----------------------------------

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


#----------------------------------
# Set the root directory and read in the necessary files
#----------------------------------

options(stringsAsFactors = FALSE)
root_dir="~/R_work/Prolonged_persistence_of_DNA_lesions/"
data_dir=paste0(root_dir,"Data/")
plots_dir=paste0(root_dir,"Plots/")
source(paste0(data_dir,"Prolonged_persistence_functions.R")) #Source functions needed for the script
resave_plots=F

Phasing_MAV_file_path=paste0(data_dir,"Phasing_results_MAVs_all")
Phasing_PVV_file_path=paste0(data_dir,"Phasing_results_PVVs_all")
ASCAT_PVV_file_path=paste0(data_dir,"ASCAT_LOH_analysis_PVVs_all")
ASCAT_MAV_file_path=paste0(data_dir,"ASCAT_LOH_analysis_MAVs_all")
SN_phasing_file_path=paste0(data_dir,"Phasing_results_MAVs_SN.csv")

lesion_seg_input_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data"
lesion_seg_output_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/output2"
genome_file="~/R_work/reference_files/genome.fa"

mutations_file=paste0(data_dir,"mutations.tsv")
summary_df_file=paste0(data_dir,"summary_df.tsv")

if(file.exists(mutations_file)&file.exists(summary_df_file)) {
  
  mutations<-read.delim(mutations_file)
  summary_table_df<-read.delim(summary_df_file)
  
} else {
  #Set data file paths
  data_sets=c("NW","MF","EM","KY","MSC_BMT","MSC_fetal","SN")
  out_list=lapply(data_sets,function(data_set) {
    print(data_set)
    files=list.files(paste0(lesion_seg_output_dir,"/",data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
    print(files)
    mut_tables_list=lapply(files,read.delim)
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
  mutations$Ref[mutations$Ref=="TRUE"]<-"T"
  mutations$Alt1[mutations$Alt1=="TRUE"]<-"T"
  mutations$Alt2[mutations$Alt2=="TRUE"]<-"T"
  
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
  sample_ref=read.csv("Individual_ref.csv")
  mutations$cat=sapply(mutations$Sample_ID,function(sample) return(sample_ref$Category[sample_ref$Sample_ID==sample]))
  mutations<-mutations%>%filter(data_set!="PR")
  
  write.table(mutations,file=mutations_file,quote=F,sep="\t",row.names = F)
  
  
  #Summary table
  summary_table_list=lapply(data_sets,function(data_set) {
    print(data_set)
    data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
    data_set_df_list=lapply(data_set_samples,function(sample){
      sample_info=get_file_paths_and_project(data_set,Sample_ID=sample)
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

#----------------------------------
# START ANALYSIS - MAVs
#----------------------------------

#Analysis of MAV phasing
MAV_mutations=mutations%>%filter(Type=="MAV")

#Fail any "removed" mutations with more than 1 negative subclade (likely to be independently-acquired mutations)
MAV_mutations$n_neg=sapply(MAV_mutations$colonies_per_negative_subclade,function(x) length(unlist(strsplit(x[!is.na(x)],","))))

load(Phasing_MAV_file_path)
Phasing_results_MAVs<-Phasing_results_MAVs[-which(duplicated(names(Phasing_results_MAVs)))] #remove any duplicates

#Now extract the basic phasing info to add to the master dataframe
MAV_phasing_summary=sapply(Phasing_results_MAVs,extract_MAV_pos_clade_phasing_summary)
names(MAV_phasing_summary)<-str_split(names(MAV_phasing_summary),pattern="\\.",simplify = T)[,1]
MAV_phasing_summary<-MAV_phasing_summary[!duplicated(names(MAV_phasing_summary))]
MAV_mutations<-left_join(MAV_mutations,data.frame(Chrom_pos=names(MAV_phasing_summary),phasing_summary=MAV_phasing_summary),by="Chrom_pos")

#Now add the phasing information from the liver samples (this is analyzed separately as they are non-clonal samples so have a different framework)
SN_phasing_summary<-read.csv(SN_phasing_file_path,stringsAsFactors=F)
MAV_mutations$phasing_summary<-sapply(1:nrow(MAV_mutations),function(i) {
  if(MAV_mutations$Sample_ID[i]%in%SN_phasing_summary$exp_ID) {
    return(SN_phasing_summary%>%filter(Chrom_pos==MAV_mutations$Chrom_pos[i])%>%pull(res))
  } else {
    return(MAV_mutations$phasing_summary[i])
  }
})

#Create simplified phasing metric, with only 3 categories
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

myPhasingColors <- c("#33A02C","#FB9A99","grey")
names(myPhasingColors) <- levels(MAV_mutations$phasing_summary)

p.MAV.0.1<-MAV_mutations%>%
  filter(Class=="removed")%>%
  mutate(phasing_summary=factor(phasing_summary))%>%#filter(Class!="FAIL")%>%
  mutate(n_neg_cat=ifelse(n_neg>1,"n_neg > 1","n_neg <= 1"))%>%
  group_by(phasing_summary)%>%
  dplyr::count(phasing_summary,n_neg_cat)%>%
  ggplot(aes(x=n_neg_cat,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="stack",col="black",linewidth=0.2)+
  scale_fill_manual(name = "MAV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="Number of negative\nsubclades",y="Number of mutations")

#Now set these mutations with > 1 negative subclade to "FAIL"
MAV_mutations%>%group_by(phasing_summary)%>%summarise(n_exclude=sum(n_neg>1 & Class=="removed"))
MAV_mutations<-MAV_mutations%>%mutate(Class=ifelse(Class=="removed"&n_neg>1,"FAIL",Class))

#Review the non-matching phasing mutations
MAV_mutations%>%
  mutate(phasing_summary=factor(phasing_summary))%>%filter(phasing_summary%in%c("Non-matching phasing confirmed")&Class!="FAIL")

#By category, the phasing outcome of all mutations (FAIL mutations included)
p.MAV.1<-MAV_mutations%>%
  mutate(phasing_summary=factor(phasing_summary))%>%#filter(Class!="FAIL")%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  group_by(phasing_summary)%>%
  dplyr::count(phasing_summary,cat)%>%
  ggplot(aes(x=cat,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="stack",col="black",linewidth=0.2)+
  scale_fill_manual(name = "MAV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  #scale_colour_discrete()+
  theme_classic()+
  my_theme+
  labs(x="Data type",
       y="Number of mutations")

#Calculate if the difference between matching and non-matching is significant
cat_comp<-MAV_mutations%>%
  filter(phasing_summary%in%c("Non-matching phasing confirmed","Same phasing confirmed"))%>%
  mutate(phasing_summary=factor(phasing_summary),cat=factor(cat,level=c("Adult_HSPC","Foetal_HSPC","Chemo_HSPC","Bronchial","Colonic","Liver")))%>%
  mutate(Class=factor(Class,levels=c("FAIL","PASS","removed","simple")))%>%
  group_by(cat,Class,phasing_summary)%>%
  dplyr::count(phasing_summary,Class,cat)%>%
  tidyr::pivot_wider(names_from = phasing_summary,values_from=n)%>%
  tidyr::complete(fill=list(`Non-matching phasing confirmed`=0,`Same phasing confirmed`=0))

cat_comp$p.value<-mapply(FUN=function(x,y){binom.test(x=x,n=(x+y),p =0.5,alternative="greater")$p.value},y=cat_comp%>%pull(`Non-matching phasing confirmed`),x=cat_comp%>%pull(`Same phasing confirmed`))
cat_comp$significance=ifelse(cat_comp$p.value<0.005,"**",ifelse(cat_comp$p.value<0.05,"*","n.s."))

#Overall (across datasets), the number of phasing confirmed vs non-matching phasing confirmed
p.MAV.2<-MAV_mutations%>%
  mutate(phasing_summary=factor(phasing_summary))%>%filter(Class!="FAIL")%>%
  #group_by(Class,phasing_summary)%>%
  dplyr::count(Class,phasing_summary)%>%
  complete(Class,phasing_summary,fill=list(n=0))%>%
  filter(phasing_summary%in%c("Non-matching phasing confirmed","Same phasing confirmed"))%>%
  mutate(Class=ifelse(Class=="PASS","PASS (Liver)",Class))%>%
  ggplot(aes(x=Class,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="dodge",col="black",linewidth=0.2)+
  scale_fill_manual(name = "MAV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="MAV Class",
       y="Number of mutations")

#By category, the number of phasing confirmed vs non-matching phasing confirmed (including FAIL mutations)
p.MAV.3.1<-MAV_mutations%>%
  #filter(data_set!="SN")%>%
  mutate(Class=ifelse(Class=="PASS","simple",Class))%>%
  mutate(phasing_summary=factor(phasing_summary))%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  group_by(Class,phasing_summary)%>%
  dplyr::count(phasing_summary,cat)%>%
  tidyr::complete(phasing_summary,cat,fill=list(n=0))%>%
  filter(phasing_summary%in%c("Non-matching phasing confirmed","Same phasing confirmed"))%>%
  ggplot(aes(x=cat,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="dodge",col="black",linewidth=0.2)+
  facet_wrap(~Class,nrow=1)+
  scale_fill_manual(name = "MAV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="MAV Class",
       y="Number of mutations")

p.MAV.3.2<-MAV_mutations%>%
  filter(data_set!="SN")%>%
  mutate(phasing_summary=factor(phasing_summary))%>%
  group_by(Class,phasing_summary)%>%
  dplyr::count(phasing_summary)%>%
  complete(phasing_summary,fill=list(n=0))%>%
  filter(phasing_summary%in%c("Non-matching phasing confirmed","Same phasing confirmed"))%>%
  ggplot(aes(x=Class,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="dodge",col="black",size=0.2)+
  scale_fill_manual(name = "MAV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="MAV Class",
       y="Number of mutations")

order=MAV_mutations%>%
  filter(cat=="Bronchial")%>%
  group_by(Sample_ID)%>%
  summarise(n=n())%>%
  mutate(smoking_status=sapply(Sample_ID,function(sample) {smoking_status=sample_ref$Smoking_status[sample_ref$Sample_ID==sample]}))%>%
  mutate(smoking_status=factor(smoking_status,levels = c("Child","Never_smoker","Ex_smoker","Current_smoker")))%>%
  arrange(smoking_status)%>%
  pull(Sample_ID)

p.MAV.4<-MAV_mutations%>%
  mutate(smoking_status=sapply(Sample_ID,function(sample) {smoking_status=sample_ref$Smoking_status[sample_ref$Sample_ID==sample]}))%>%
  mutate(smoking_status=factor(smoking_status,levels = c("Child","Never_smoker","Ex_smoker","Current_smoker")),Sample_ID=factor(Sample_ID,levels=order))%>%
  filter(phasing_summary!="Non-matching phasing confirmed"&Class!="FAIL")%>%
  filter(cat=="Bronchial")%>%
  dplyr::count(Sample_ID)%>%
  complete(Sample_ID,fill=list(n=0))%>%
  mutate(smoking_status=sapply(Sample_ID,function(sample) {smoking_status=sample_ref$Smoking_status[sample_ref$Sample_ID==sample]}))%>%
  ggplot(aes(x=Sample_ID,y=n,fill=smoking_status))+
  geom_bar(stat="identity")+
  theme_classic()+
  my_theme

#Number of bronchial MAVs that include an MNV
p.MAV.5<-MAV_mutations%>%
  filter(phasing_summary!="Non-matching phasing confirmed"&Class!="FAIL")%>%
  filter(cat=="Adult_HSPC")%>%
  mutate(Mut_types=paste(Mut_type1,Mut_type2,sep="-"))%>%
  mutate(Contains_MNV=ifelse(grepl("MNV",Mut_types),"Contains MNV","SNVs only"))%>%
  group_by(Contains_MNV)%>%
  summarise(n=n())%>%
  ggplot(aes(x=Contains_MNV,y=n,fill=Contains_MNV))+
  geom_bar(stat="identity")+
  theme_classic()+
  labs(title="Number of bronchial MAVs that include an MNV")

#MAV signature analysis
multi_SNVs<-MAV_mutations%>%
  filter(phasing_summary!="Non-matching phasing confirmed"&Class!="FAIL")%>%
  filter(cat!="Bronchial"|Class!="removed")%>% #Remove the "removed" bronchial mutations
  filter(Mut_type1=="SNV"&Mut_type2=="SNV")

#Create consistent order for the multi-allelic mutations: transitions first (if either is), or else transversion to A first
multi_SNVs_ordered=dplyr::bind_rows(mapply(function(mut1,mut2) {
  pyr=c("T","C")
  if(substr(mut1,5,5)%in%pyr & !substr(mut2,5,5)%in%pyr) {
    return(data.frame(sub1=substr(mut1,3,5),sub2=substr(mut2,3,5)))
  } else if(!substr(mut1,5,5)%in%pyr & substr(mut2,5,5)%in%pyr) {
    return(data.frame(sub1=substr(mut2,3,5),sub2=substr(mut1,3,5)))
  } else if(!substr(mut1,5,5)%in%pyr & !substr(mut2,5,5)%in%pyr & substr(mut1,5,5)=="A") {
    return(data.frame(sub1=substr(mut1,3,5),sub2=substr(mut2,3,5)))
  } else if(!substr(mut1,5,5)%in%pyr & !substr(mut2,5,5)%in%pyr & substr(mut1,5,5)!="A") {
    return(data.frame(sub1=substr(mut2,3,5),sub2=substr(mut1,3,5))) 
  }
},mut1=multi_SNVs$mut_profile_1,mut2=multi_SNVs$mut_profile_2,SIMPLIFY = F))

#Create a "mut_combination" variable
multi_SNVs=cbind(multi_SNVs,multi_SNVs_ordered)
multi_SNVs$mut_combination=paste(multi_SNVs$sub1,multi_SNVs$sub2,sep="_")

#Divide by all individual samples
p.MAV.6<-multi_SNVs %>%
  mutate(mut_combination=as.factor(mut_combination)) %>%
  dplyr::count(cat,mut_combination)%>%
  complete(cat,mut_combination,fill=list(n=0))%>%
  ggplot(aes(x=mut_combination,y=n,fill=mut_combination)) +
  geom_bar(stat="identity",col="black",size=0.4) +
  my_theme+
  facet_wrap(~cat,nrow = 2) +
  theme(axis.text.x = element_text(size=8,angle = 90))

#Breakdown of MAV signature by sample (bronchial samples only)
p.MAV.7<-multi_SNVs %>%
  filter(cat=="Bronchial")%>%
  mutate(mut_combination=as.factor(mut_combination)) %>%
  dplyr::count(Sample_ID,mut_combination)%>%
  complete(Sample_ID,mut_combination,fill=list(n=0))%>%
  ggplot(aes(x=mut_combination,y=n,fill=mut_combination)) +
  geom_bar(stat="identity",col="black",size=0.4) +
  my_theme+
  facet_wrap(~Sample_ID,nrow = 2) +
  theme(axis.text.x = element_text(size=8,angle = 90))

#Now plot MAV signatures as 96-profile plots
#Design function to plot a bespoke 96-profile for 2 x SNV combinations at a single site
plot_MAV_96profile=function(multi_SNV_table,cat_types=NULL,relative=F,return_vec=F,compressed=F,ymax=0.1) {
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  #COLORS6=c("#2EBAED" ,"#000000" ,"#DE1C14" ,"#D4D2D2" ,"#ADCC54" ,"#F0D0CE")
  COLORS6=c("#66C2A5" ,"#FC8D62" ,"#8DA0CB" ,"#E78AC3" ,"#A6D854" ,"#FFD92F") #Different colours to the normal 96-profile sig to clarify different meaning
  if(!is.null(cat_types)) {
    mut_table=multi_SNV_table%>%
      mutate(mut_combination=as.factor(mut_combination))%>%
      mutate(trinuc_ref_py=as.factor(trinuc_ref_py))%>%
      filter(cat%in%cat_types)%>%
      mutate(cat=factor(cat,levels=cat_types))%>%
      dplyr::count(trinuc_ref_py,mut_combination,cat)%>%
      tidyr::complete(trinuc_ref_py,mut_combination,cat,fill=list(n=0))%>%
      filter(substr(trinuc_ref_py,2,2)==substr(mut_combination,1,1))
  } else {
    mut_table=multi_SNV_table%>%
      mutate(mut_combination=as.factor(mut_combination))%>%
      mutate(trinuc_ref_py=as.factor(trinuc_ref_py))%>%
      dplyr::count(trinuc_ref_py,mut_combination)%>%
      tidyr::complete(trinuc_ref_py,mut_combination,fill=list(n=0))%>%
      filter(substr(trinuc_ref_py,2,2)==substr(mut_combination,1,1))
  }
  if(relative) {
    mut_freqs=pivot_wider(mut_table,names_from = cat,values_from = n)
    mut_freqs[,3:ncol(mut_freqs)]<-apply(mut_freqs[,3:ncol(mut_freqs)],2,function(x) return(x/sum(x)))
    mut_table=gather(mut_freqs,-trinuc_ref_py,-mut_combination,key="cat",value="n")%>%
      mutate(cat=factor(cat,levels=cat_types))
  }
  mut_vec=mut_table$n; names(mut_vec)=paste(mut_table$mut_combination,mut_table$trinuc_ref_py,sep=",")
  if(return_vec) {stop(return(mut_vec))}
  ymax=max(mut_table$n)
  
  #THE PLOT
  plot<-mut_table%>%
    mutate(context=paste0(substr(trinuc_ref_py,1,1),".",substr(trinuc_ref_py,3,3)))%>%
    ggplot(aes(x = context, y = n, fill = mut_combination))+
    geom_bar(stat = "identity",width=0.8)+
    scale_fill_manual(values = COLORS6)+
    facet_grid(cat~mut_combination,scales = "free_y")+
    ylab(ifelse(relative,"Relative contribution","Frequency"))+
    coord_cartesian(ylim = c(0, ymax))+
    scale_y_continuous(breaks = NULL) + 
    guides(fill = FALSE)+
    theme_bw()+
    theme(axis.title.y = element_text(size = 8,vjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0.5,"lines"))
    
  if(compressed){
    plot+theme(rect=element_rect(size=0.3),
               axis.title.y = element_text(size = 7,vjust = 1),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(), 
               axis.title.x = element_text(size = 7),
               axis.text.x = element_blank(),
               strip.text.x = element_text(size = 5),
               strip.text.y = element_text(size = 5),
               panel.grid.major.x = element_blank(),
               panel.spacing.x = unit(0.2,"lines"))
  }
  return(plot)
}

p.MAV.8.1<-plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_types = c("Bronchial","Liver","Adult_HSPC","Chemo_HSPC"))
p.MAV.8.2<-plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_types="Bronchial")
p.MAV.8.3<-plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_types="Chemo_HSPC",compressed=T)
p.MAV.8.4<-plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_types="Adult_HSPC",compressed=T)
p.MAV.8.4.2<-plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_types="Liver",compressed=T)

#'FAIL' MAV signature analysis
multi_SNVs_fail<-MAV_mutations%>%
  filter(phasing_summary=="Non-matching phasing confirmed"|Class=="FAIL")%>%
  #filter(cat!="Bronchial"|Class!="removed")%>% #Remove the "removed" bronchial mutations
  filter(Mut_type1=="SNV"&Mut_type2=="SNV")

#Create consistent order for the multi-allelic mutations: transitions first (if either is), or else transversion to A first
multi_SNVs_fail_ordered=dplyr::bind_rows(mapply(function(mut1,mut2) {
  pyr=c("T","C")
  if(substr(mut1,5,5)%in%pyr & !substr(mut2,5,5)%in%pyr) {
    return(data.frame(sub1=substr(mut1,3,5),sub2=substr(mut2,3,5)))
  } else if(!substr(mut1,5,5)%in%pyr & substr(mut2,5,5)%in%pyr) {
    return(data.frame(sub1=substr(mut2,3,5),sub2=substr(mut1,3,5)))
  } else if(!substr(mut1,5,5)%in%pyr & !substr(mut2,5,5)%in%pyr & substr(mut1,5,5)=="A") {
    return(data.frame(sub1=substr(mut1,3,5),sub2=substr(mut2,3,5)))
  } else if(!substr(mut1,5,5)%in%pyr & !substr(mut2,5,5)%in%pyr & substr(mut1,5,5)!="A") {
    return(data.frame(sub1=substr(mut2,3,5),sub2=substr(mut1,3,5))) 
  }
},mut1=multi_SNVs_fail$mut_profile_1,mut2=multi_SNVs_fail$mut_profile_2,SIMPLIFY = F))

#Create a "mut_combination" variable
multi_SNVs_fail=cbind(multi_SNVs_fail,multi_SNVs_fail_ordered)
multi_SNVs_fail$mut_combination=paste(multi_SNVs_fail$sub1,multi_SNVs_fail$sub2,sep="_")

p.MAV.8.5<-plot_MAV_96profile(multi_SNV_table = multi_SNVs_fail,cat_types=c("Bronchial","Liver","Adult_HSPC","Chemo_HSPC"),relative=T)


#Get the MAV signatures in a savable data frame for comparison to anticipated MAV signatures
MAV_vec_list=lapply(c("Bronchial","Chemo_HSPC","Adult_HSPC","Liver"),function(cat_type) {
  mut_vec=plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_type=cat_type,relative = T,return_vec = T)
  return(mut_vec)
})
names(MAV_vec_list)<-c("Bronchial","Chemo_HSPC","Adult_HSPC","Liver")
MAV_sig_df=as.data.frame(dplyr::bind_cols(MAV_vec_list))
rownames(MAV_sig_df)<-names(MAV_vec_list[[1]])
write.csv(x = MAV_sig_df,file = "MAV_sig_df.csv",row.names = T)

#Analyse the NUMBERS of MAVs in different samples
#Generate summary table of numbers of MAV mutations
n_MAV_summary_table=MAV_mutations%>%
  filter(Class!="FAIL")%>%
  mutate(Sample_ID=factor(Sample_ID,levels=summary_table_df$Sample_ID))%>%
  group_by(Sample_ID)%>%
  summarise(n_MAV=n())

#Join this to the sample summary table
summary_table_df_MAV<-summary_table_df%>%
  left_join(n_MAV_summary_table)%>%
  mutate(n_MAV=replace_na(n_MAV,0))%>%
  mutate(smoking_status=factor(smoking_status,levels=c("unknown","Current_smoker","Ex_smoker","Never_smoker","Child")))

#Plot the number of MAVs vs the number of post-developmental nodes for all samples
p.MAV.9.1<-summary_table_df_MAV%>%
  filter(cat!="Chemo_HSPC")%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=post_development_nodes,y=n_MAV,col=cat,label=Sample_ID,shape=smoking_status))+
  geom_point(size=2,alpha=0.5)+
  geom_smooth(method="lm",data=summary_table_df_MAV%>%filter(cat=="Adult_HSPC"),col=RColorBrewer::brewer.pal(3,"Set1")[1])+
  geom_text_repel(show.legend=F,fontface="italic",size=2)+
  theme_classic()+
  my_theme+
  scale_color_brewer(palette="Set1")+
  theme(legend.key.size = unit(0.3, "cm"))+
  labs(x="Number of post-development internal nodes in phylogeny",
       y="Number of MAVs detected",
       col="Data type",
       shape="Smoking status")
  

p.MAV.9.2<-summary_table_df_MAV%>%filter(data_set%in%c("EM","MF","NW","MSC_BMT")&!Sample_ID%in%c("CB001_3_01","CB002_2_01","PX001_2_01"))%>%
  ggplot(aes(x=post_development_nodes,y=n_MAV,label=Sample_ID))+
  geom_point(col=RColorBrewer::brewer.pal(3,"Set1")[1],alpha=0.5)+
  #geom_text_repel(check_overlap = T,nudge_y = 1.5)+
  geom_smooth(method="lm",col="darkgrey")+
  theme_classic()+
  my_theme+
  labs(x="Number of post-development\ninternal nodes in phylogeny",
       y="Number of MAVs detected")

p.MAV.9.3<-summary_table_df_MAV%>%filter(data_set=="SN")%>%
  ggplot(aes(x=post_development_nodes,y=n_MAV,label=Sample_ID))+
  geom_jitter(col=RColorBrewer::brewer.pal(4,"Set1")[4],alpha=0.5,width=0.1,height=0.1)+
  geom_smooth(method="lm",col="darkgrey")+
  theme_classic()+
  my_theme+
  labs(x="Number of post-development\ninternal nodes in phylogeny",
       y="Number of MAVs detected")

p.MAV.9.4<-summary_table_df_MAV%>%filter(data_set=="KY")%>%
  ggplot(aes(x=post_development_nodes,y=n_MAV,label=Sample_ID,shape=smoking_status))+
  geom_jitter(col=RColorBrewer::brewer.pal(4,"Set1")[4],alpha=0.5,width=0.1,height=0.1)+
  geom_smooth(aes(x=post_development_nodes,y=n_MAV),method="lm",col="darkgrey",inherit.aes=F)+
  theme_classic()+
  scale_y_continuous(limits=function(x) c(0,x[2]))+
  my_theme+
  labs(x="Number of post-development\ninternal nodes in phylogeny",
       y="Number of MAVs detected",
       shape="Smoking status")

#Do a linear regression on the adult HSPC dataset for number of MAVs vs number of post-development nodes
mylm=lm(n_MAV~post_development_nodes,data=summary_table_df_MAV%>%filter(data_set%in%c("EM","MF","NW","MSC_BMT")&!Sample_ID%in%c("CB001_3_01","CB002_2_01","PX001_2_01")))
summary(mylm)

mylm.liver=lm(n_MAV~post_development_nodes,data=summary_table_df_MAV%>%filter(data_set=="SN"))
summary(mylm.liver)

mylm.lung=lm(n_MAV~post_development_nodes,data=summary_table_df_MAV%>%filter(data_set=="KY"))
summary(mylm.lung)

p.MAV.9.5<-summary_table_df_MAV%>%
  filter(cat!="Chemo_HSPC" & cat !="Foetal_HSPC")%>%
  mutate(cat=factor(cat,levels=c("Adult_HSPC","Liver","Bronchial","Chemo_HSPC")))%>%
  ggplot(aes(x=post_development_nodes,y=n_MAV,col=cat,label=Sample_ID,shape=smoking_status))+
  geom_jitter(alpha=0.5,width=0.1,height=0.1)+
  geom_smooth(aes(x=post_development_nodes,y=n_MAV),method="lm",col="darkgrey",inherit.aes=F)+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_y_continuous(limits=function(x) c(0,x[2]))+
  #scale_y_log10()+
  #scale_x_log10()+
  facet_grid(cols=vars(cat),scales="free")+
  #facet_grid(cols=vars(cat))+
  my_theme+
  labs(x="Number of post-development\ninternal nodes in phylogeny",
       y="Number of MAVs detected",
       col="Sample type",
       shape="Smoking status")

p.MAV.9.6<-summary_table_df_MAV%>%
  filter(cat!="Chemo_HSPC" & cat !="Foetal_HSPC")%>%
  mutate(cat=factor(cat,levels=c("Adult_HSPC","Liver","Bronchial","Chemo_HSPC")))%>%
  ggplot(aes(x=post_development_nodes,y=n_MAV,col=cat,label=Sample_ID,shape=smoking_status))+
  geom_jitter(alpha=0.5,width=0.1,height=0.1)+
  geom_smooth(aes(x=post_development_nodes,y=n_MAV,col=cat),method="lm",inherit.aes=F)+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_y_log10()+
  scale_x_log10()+
  my_theme+
  labs(x="Number of post-development\ninternal nodes in phylogeny",
       y="Number of MAVs detected",
       col="Sample type",
       shape="Smoking status")

p.MAV.9.7<-summary_table_df_MAV%>%
  mutate(MAVs_per_PDN=n_MAV/post_development_nodes)%>%
  filter(cat!="Foetal_HSPC")%>%
  mutate(cat=factor(cat,levels=c("Adult_HSPC","Liver","Bronchial","Chemo_HSPC")))%>%
  ggplot(aes(x=cat,y=MAVs_per_PDN,col=cat,shape=smoking_status))+
  geom_boxplot(aes(x=cat,y=MAVs_per_PDN),inherit.aes = F,outlier.shape = NA)+
  geom_jitter(width = 0.2,height=0,alpha=0.4)+
  scale_color_brewer(palette = "Set1")+
  scale_x_discrete(labels=function(x) {stringr::str_replace(x,"_","\n")})+
  theme_classic()+
  my_theme+
  labs(x="Tissue type",
       y="Number of MAVs per\npost-development node",
       col="Sample type",
       shape="Smoking status")

p.MAV.10.1<-MAV_mutations%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  filter(!is.na(lesion_duration)&Class=="removed"&lesion_duration!=0)%>%
  ggplot(aes(x=cat,y=lesion_duration,fill=cat)) +
  #geom_boxplot(outlier.size = 0.5)+
  geom_violin(width=0.7)+
  geom_jitter(col="grey",width=0.2,size=1,alpha=0.5)+
  #facet_grid(~cat,scales="free_x",space="free")+
  my_theme +
  scale_color_brewer(palette = "Set2")+
  theme(legend.position = "none",axis.text.x=element_text(size=7)) +
  scale_y_log10()+
  labs(y="Minimum molecular\nlesion duration (MMLD)",x="Phylogeny category")

p.MAV.11.1<-MAV_mutations%>%
  filter(Class!="FAIL")%>%
  group_by(no_of_cell_divisions)%>%
  summarise(n=n())%>%
  ggplot(aes(x=no_of_cell_divisions,y=n))+
  geom_bar(stat="identity")+
  my_theme+
  labs(x="Minimum cell divisions",y="Count")

p.MAV.11.1<-MAV_mutations%>%
  filter(Class!="FAIL")%>%
  mutate(no_of_cell_divisions=no_of_cell_divisions-1)%>%
  ggplot(aes(x=no_of_cell_divisions,fill=phasing_summary))+
  geom_bar(stat="count",position="stack",width=0.7)+
  scale_x_continuous(breaks=seq(1,5,1))+
  scale_fill_manual(name = "MAV phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="Minimum cell divisions",y="Count")

if(resave_plots){
  ggsave(p.MAV.0.1,filename=paste0(plots_dir,"MAV_outcomes_removed_by_nneg.pdf"),width=2,height=2)
  ggsave(p.MAV.1+theme(legend.position = "none"),filename=paste0(plots_dir,"MAV_outcomes_80321.pdf"),width=2,height=2)
  ggsave(p.MAV.3.1+theme(legend.position = "none"),filename=paste0(plots_dir,"MAV_outcomes_by_class_80321.pdf"),width=3,height = 2)
  ggsave(p.MAV.3.2,filename=paste0(plots_dir,"MAV_outcomes_by_class_only_80321.pdf"),width=2.5,height=2)
  ggsave(p.MAV.4,filename=paste0(plots_dir,"MAV_Bronchial_per_sample.pdf"),width=5,height=2)
  ggsave(p.MAV.5,filename=paste0(plots_dir,"MAV_MNVs.pdf"),width=8,height=4)
  ggsave(p.MAV.6,filename=paste0(plots_dir,"MAV_6profile_MutSigs_by_cat.pdf"),width=6,height=4)
  ggsave(p.MAV.8.1,filename=paste0(plots_dir,"MAV_signature_all.pdf"),width=7,height=4)
  ggsave(p.MAV.8.2,filename=paste0(plots_dir,"MAV_signature_bronchial.pdf"),width=7,height=2)
  ggsave(p.MAV.8.3,filename=paste0(plots_dir,"MAV_signature_chemo.pdf"),width=3,height=1.2)
  ggsave(p.MAV.8.4,filename=paste0(plots_dir,"MAV_signature_adult.pdf"),width=3,height=1.2)
  ggsave(p.MAV.8.5,filename=paste0(plots_dir,"MAV_fail_signature_all.pdf"),width=7,height=4)
  ggsave(p.MAV.9.1,filename=paste0(plots_dir,"MAV_numbers_by_PDN.pdf"),width=4,height=2)
  ggsave(p.MAV.9.2,filename=paste0(plots_dir,"MAV_numbers_by_PDN_AdultHSPC.pdf"),width=2,height=2)
  ggsave(p.MAV.9.3,filename=paste0(plots_dir,"MAV_numbers_by_PDN_Liver.pdf"),width=2,height=2)
  ggsave(p.MAV.9.4,filename=paste0(plots_dir,"MAV_numbers_by_PDN_Bronchial.pdf"),width=2,height=2)
  ggsave(p.MAV.9.5,filename=paste0(plots_dir,"MAV_numbers_by_PDN_faceted.pdf"),width=5,height=2)
  ggsave(p.MAV.9.6,filename=paste0(plots_dir,"MAV_numbers_by_PDN_logscale.pdf"),width=2.5,height=2)
  ggsave(p.MAV.9.7,filename=paste0(plots_dir,"MAVs_per_PDN.pdf"),width=3,height=2)
  ggsave(p.MAV.10.1,filename=paste0(plots_dir,"MAV_lesion_durations_by_cat.pdf"),width=2.5,height=2)
  ggsave(p.MAV.11.1,filename = paste0(plots_dir,"MAV_no_of_cell_divisions.pdf"),width = 2.5,height=2)
}


#----------------------------------
# START ANALYSIS - PVVs
#----------------------------------

PVV_mutations<-mutations%>%filter(Type=="PVV")

#Load the data to confirm no LOH
load(ASCAT_PVV_file_path)
ASCAT_LOH_analysis_PVVs<-ASCAT_LOH_analysis_PVVs[-which(duplicated(names(ASCAT_LOH_analysis_PVVs)))] #remove any duplicates
PVV_mutations<-left_join(PVV_mutations,data.frame(Chrom_pos=names(ASCAT_LOH_analysis_PVVs),ASCAT_result=unlist(ASCAT_LOH_analysis_PVVs)),by="Chrom_pos")

#Extract the phasing info regarding the positive subclades
load(Phasing_PVV_file_path)
Phasing_results_PVVs<-Phasing_results_PVVs[-which(duplicated(names(Phasing_results_PVVs)))] #remove any duplicates

#Add phasing information on the positive subclades
PVV_phasing_summary=sapply(Phasing_results_PVVs,extract_PVV_pos_clade_phasing_summary); names(PVV_phasing_summary)<-str_split(names(PVV_phasing_summary),pattern="\\.",simplify=T)[,1]
PVV_mutations<-left_join(PVV_mutations,data.frame(Chrom_pos=names(PVV_phasing_summary),PVV_pos_clade_phasing=PVV_phasing_summary),by="Chrom_pos")
#Simplify the output to 3 categories: same-phasing, non-matching phasing, or unable to confirm
PVV_mutations$PVV_pos_clade_phasing<-sapply(PVV_mutations$PVV_pos_clade_phasing,function(x) {
  if(grepl("Non-matching phasing confirmed",x)) {
    res<-"Non-matching phasing confirmed"
  } else if(grepl("Same phasing confirmed",x)){
    res<-"Same phasing confirmed"
  } else if(is.na(x)) {
    res<-NA
  } else {
    res<-"Unable to confirm phasing"
  }
  return(res)
})
PVV_mutations$PVV_pos_clade_phasing<-factor(PVV_mutations$PVV_pos_clade_phasing,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing"))
PVV_mutations<-mutate(PVV_mutations,Class=ifelse(is.na(lesion_node)|ASCAT_result=="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected","FAIL","PASS"))

#Add in information on the negative subclades
PVV_both_alleles_summary=sapply(Phasing_results_PVVs,extract_PVV_neg_clade_phasing_summary)
Neg_clades_phasing_df=data.frame(Chrom_pos=names(Phasing_results_PVVs),result=sapply(PVV_both_alleles_summary,function(x) return(x[1])))

#Simplify the output to only "No loss of mutant allele", "Unable to confirm", "Likely LOH" or "Positive clade have non-matching phasing"
Neg_clades_phasing_df$basic_result=sapply(Neg_clades_phasing_df$result,function(result) {
  if(result%in%c("Alt allele reads present","Both alleles confirmed with reference allele","Both alleles confirmed with reference allele in at least one subclade")){
    return("No LOH")
  } else if(grepl("with maximum",result)){
    n_reads<-as.numeric(sapply(str_extract_all(result,"\\d"),paste,collapse=""))
    if(n_reads<=4) {return("Unable to confirm")} else if(n_reads>4) {return("LOH")}
  } else if(result=="Positive clades have non-matching phasing"){return(NA)} else {return("Unable to confirm")}
})
PVV_mutations=left_join(PVV_mutations,Neg_clades_phasing_df[,c("Chrom_pos","result","basic_result")],by="Chrom_pos") #Join this to the mutations dataframe
PVV_mutations%>%filter(Class!="FAIL"&ASCAT_result!="LOH"&basic_result!="LOH")%>%pull(PVV_pos_clade_phasing)%>%table()

myLOHColors<-c(rev(RColorBrewer::brewer.pal(3,"Accent")),"grey")
names(myLOHColors)<-c("LOH","Unable to confirm","No LOH","Homozygous in all clades")

#Extract the adult HSPC fail mutations for comparatory signature analysis
FAIL_PVVs<-PVV_mutations%>%
  filter(Class=="FAIL"&cat=="Adult_HSPC")%>%
  group_by(cat)%>%
  pull(mut_ref1)
write.vcf(FAIL_PVVs,vcf_path = "Blood_fail_PVVs.vcf",vcf_header_path = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt")

#Introduce a more stringent threshold for the minimum depth of the negative subclade - needs to be at least ≥13 (not 10 as previously). This loses
PVV_mutations$max_neg_clade_depth=sapply(PVV_mutations$depth_per_negative_subclade,function(x) max(as.numeric(unlist(strsplit(x,split=",")))))
PVV_mutations<-PVV_mutations%>%filter(max_neg_clade_depth>12|is.na(max_neg_clade_depth))

p.PVV.1.1<-PVV_mutations%>%
  mutate(ASCAT_result=factor(ASCAT_result))%>%filter(Class!="FAIL")%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  group_by(ASCAT_result)%>%
  dplyr::count(ASCAT_result,cat)%>%
  ggplot(aes(x=cat,y=n,fill=ASCAT_result))+
  geom_bar(stat="identity",position="stack",col=NA,size=0.2)+
  scale_fill_manual(name = "ASCAT LOH result",values=myLOHColors,labels = function(x) str_wrap(x, width = 9))+
  theme_classic()+
  my_theme+
  labs(x="Data type",
       y="Number of mutations")

p.PVV.1.2<-PVV_mutations%>%
  mutate(PVV_pos_clade_phasing=factor(PVV_pos_clade_phasing))%>% #filter(Class!="FAIL")%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  group_by(PVV_pos_clade_phasing)%>%
  dplyr::count(PVV_pos_clade_phasing,cat)%>%
  ggplot(aes(x=cat,y=n,fill=PVV_pos_clade_phasing))+
  geom_bar(stat="identity",position="stack",col=NA,size=0.2)+
  scale_fill_manual(name = "PVV phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="Data type",
       y="Number of mutations")

p.PVV.1.3<-PVV_mutations%>%
  mutate(ASCAT_result=factor(ASCAT_result))%>%filter(Class!="FAIL")%>%
  #filter(ASCAT_result=="LOH")%>%
  mutate(basic_result=factor(basic_result))%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  #filter(basic_result=="Likely LOH")%>%
  dplyr::count(basic_result,cat)%>%
  ggplot(aes(x=cat,y=n,fill=basic_result))+
  geom_bar(stat="identity",position="stack",col=NA,size=0.2)+
  scale_fill_manual(name = "LOH result",values=myLOHColors,labels = function(x) str_wrap(x, width = 9))+
  theme_classic()+
  my_theme+
  labs(x="Data type",
       y="Number of mutations")

p.PVV.1.4<-PVV_mutations%>%
  mutate(ASCAT_result=factor(ASCAT_result))%>%filter(Class!="FAIL")%>%
  filter(ASCAT_result=="LOH")%>%
  mutate(basic_result=factor(basic_result))%>%
  #filter(basic_result=="Likely LOH")%>%
  dplyr::count(basic_result)%>%
  ggplot(aes(x="ASCAT LOH",y=n,fill=basic_result))+
  geom_bar(stat="identity",position="stack",col=NA,size=0.2)+
  scale_fill_manual(name = "LOH result",values=myLOHColors,labels = function(x) str_wrap(x, width = 9))+
  theme_classic()+
  my_theme+
  labs(x="",
       y="Number of mutations")

p.PVV.1.5<-PVV_mutations%>%
  mutate(ASCAT_result=factor(ASCAT_result))%>%filter(Class!="FAIL")%>%
  #filter(ASCAT_result=="LOH")%>%
  mutate(basic_result=factor(basic_result))%>%
  filter(basic_result=="LOH")%>%
  dplyr::count(ASCAT_result)%>%
  ggplot(aes(x="Read count LOH",y=n,fill=ASCAT_result))+
  geom_bar(stat="identity",position="stack",col=NA,size=0.2)+
  scale_fill_manual(name = "ASCAT result",values=myLOHColors,labels = function(x) str_wrap(x, width = 9))+
  theme_classic()+
  my_theme+
  labs(x="",
       y="Number of mutations")

p.PVV.2<-PVV_mutations%>%
  mutate(PVV_pos_clade_phasing=factor(PVV_pos_clade_phasing))%>%filter(Class!="FAIL")%>%
  group_by(Class,PVV_pos_clade_phasing)%>%
  dplyr::count(PVV_pos_clade_phasing)%>%
  complete(PVV_pos_clade_phasing,fill=list(n=0))%>%
  filter(PVV_pos_clade_phasing%in%c("Non-matching phasing confirmed","Same phasing confirmed"))%>%
  ggplot(aes(x=Class,y=n,fill=PVV_pos_clade_phasing))+
  geom_bar(stat="identity",position="dodge",col=NA,size=0.2)+
  scale_fill_manual(name = "PVV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="PVV Class",
       y="Number of mutations")

p.PVV.2.1<-PVV_mutations%>%
  filter(Class=="PASS")%>%
  select(ASCAT_result,basic_result)%>%
  mutate(ASCAT_result=ifelse(ASCAT_result=="Homozygous in all clades","No LOH",ASCAT_result))%>%
  table()%>%
  as.data.frame()%>%
  mutate(ASCAT_result=factor(ASCAT_result,levels=c("LOH","No LOH","Unable to confirm")),Read_based_result=factor(basic_result,levels=c("LOH","No LOH","Unable to confirm")))%>%
  ggplot(aes(x=ASCAT_result,y=Read_based_result,fill=Freq))+
  geom_tile()+
  scale_fill_gradient(low="lightgrey",high="orange")+
  geom_text(aes(label = Freq), vjust = 1)+
  coord_flip()+
  my_theme

ggsave(p.PVV.2.1,filename = paste0(plots_dir,"ASCATvsReads_confusionmat.pdf"),width = 2.5,height=2)

#By category, the number of phasing confirmed vs non-matching phasing confirmed (including FAIL mutations)
p.PVV.3<-PVV_mutations%>%
  mutate(PVV_pos_clade_phasing=factor(PVV_pos_clade_phasing))%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  group_by(Class,PVV_pos_clade_phasing)%>%
  dplyr::count(PVV_pos_clade_phasing,cat)%>%
  complete(PVV_pos_clade_phasing,cat,fill=list(n=0))%>%
  filter(PVV_pos_clade_phasing%in%c("Non-matching phasing confirmed","Same phasing confirmed"))%>%
  ggplot(aes(x=cat,y=n,fill=PVV_pos_clade_phasing))+
  geom_bar(stat="identity",position="dodge",col=NA,size=0.2)+
  #facet_wrap(~Class,nrow=1)+
  scale_fill_manual(name = "PVV_phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="Data type",
       y="Number of mutations")


p.PVV.3.1<-PVV_mutations%>%
  filter(lesion_duration>0)%>%
  mutate(PVV_pos_clade_phasing=factor(PVV_pos_clade_phasing))%>%
  mutate(no_of_cell_divisions=no_of_cell_divisions-1)%>%
  ggplot(aes(x=no_of_cell_divisions,fill=PVV_pos_clade_phasing))+
  geom_bar(stat="count",position="stack",width=0.7)+
  scale_x_continuous(breaks=seq(2,10,1))+
  geom_vline(xintercept = 6.5,linetype=2)+
  scale_fill_manual(name = "PVV phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="Minimum cell divisions",y="Count") +
  annotate("rect", xmin = 6.5, xmax = 11.5, ymin = 0, ymax = 400, 
           alpha = .3)

ggsave(p.PVV.3.1,filename = paste0(plots_dir,"PVV_no_of_divisions.pdf"),width = 2.5,height=2)

#REMOVE THE MUTATIONS FROM OTHER MECHANISMS FOR DOWNSTREAM ANALYSIS (i.e. signatures etc.)
#Note that some "non-genuine" PVVs (by a variety of mechanisms) will still remain as not all individual mutations were assessable
PVV_mutations<-PVV_mutations%>%filter(Class!="FAIL" & ASCAT_result!="LOH"& basic_result!="LOH" & PVV_pos_clade_phasing!="Non-matching phasing confirmed"&no_of_cell_divisions<=6)

#Or select only those that most stringently are not excluded
#PVV_mutations<-PVV_mutations%>%filter(Class!="FAIL" & ASCAT_result=="No LOH"& basic_result=="No LOH" & PVV_pos_clade_phasing=="Same phasing confirmed")
subs=c("C>A","C>G","C>T","T>A","T>C","T>G")
p.PVV.4<-PVV_mutations%>%
  mutate(Sub1=mapply(FUN=function(x,trinuc_ref_py) {if(nchar(x)!=3) {
    return("INDEL")
  } else if(x=="C>T" & substr(trinuc_ref_py,2,3)=="CG"){
    return("C>T at CpG")
  } else if(x=="C>T" & substr(trinuc_ref_py,2,3)!="CG"){
    return("C>T other")
  } else {
    return(x)
  }
  }, x=PVV_mutations$Sub1,trinuc_ref_py=PVV_mutations$trinuc_ref_py,SIMPLIFY = T))%>%
  filter(Type=="PVV")%>%
  #mutate(Cat=ifelse(Sample_ID=="PX001","PX001","Other"))%>%
  mutate(Sub1=factor(Sub1,levels = c("C>A","C>G","C>T at CpG","C>T other","T>A","T>C","T>G","INDEL")))%>%
  dplyr::count(cat,Sub1)%>%
  complete(cat,Sub1,fill=list(n=0))%>%
  ggplot(aes(x=Sub1,y=n,fill=Sub1)) +
  geom_bar(stat="identity",col="black",size=0.4) +
  my_theme+
  #geom_text(aes(label=n),vjust=-0.5,size=3)+
  theme(legend.position = "none")+
  #scale_y_continuous(breaks=seq(0,150,10),limits=c(0,150)) +
  scale_x_discrete(drop=FALSE)+
  scale_fill_manual(values=palette8,drop=FALSE)+
  facet_wrap(~cat,scales="free",ncol=4) +
  labs(title = "Phylogeny-violating variant (PVV) signatures",
       fill="Mutation type",
       x=NULL,
       y="Count")

#Display PVV MLDs by cat, using jittered scatter plot and logarithmic scale
p.PVV.5.1<-PVV_mutations%>%
  #filter(lesion_duration<400)%>%
  ggplot(aes(x=Sample_ID,y=lesion_duration,col=Sample_ID)) +
  geom_jitter(width = 0.2,height=0)+#,alpha=0.5) +
  facet_wrap(~cat,scales="free",ncol=2)+
  my_theme +
  theme(legend.position = "none") +
  labs(title="Minimum molecular lesion duration")+
  scale_y_log10()

#Display PVV MLDs by cat, using boxplot and logarithmic scale
p.PVV.5.2<-PVV_mutations%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  #filter(lesion_duration<400)%>%
  ggplot(aes(x=Sample_ID,y=lesion_duration,col=cat)) +
  geom_boxplot(outlier.size = 0.5,outlier.shape = NA)+
  geom_jitter(col="grey",width=0.1,alpha=0.3,size=0.5)+
  facet_grid(~cat,scales="free_x",space="free")+
  my_theme +
  scale_color_brewer(palette = "Set2")+
  theme(legend.position = "none",axis.text.x=element_text(angle=90,size=5)) +
  labs(y="Minimum molecular lesion duration (MMLD)",x="Sample")+
  scale_y_log10()

p.PVV.5.3<-PVV_mutations%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  #filter(lesion_duration<400)%>%
  ggplot(aes(x=cat,y=lesion_duration,fill=cat)) +
  #geom_boxplot(outlier.size = 0.5)+
  geom_violin()+
  geom_jitter(col="grey",width=0.3,alpha=0.4,size=0.6)+
  #facet_grid(~cat,scales="free_x",space="free")+
  my_theme +
  scale_color_brewer(palette = "Set2")+
  theme(legend.position = "none",axis.text.x=element_text(size=7)) +
  scale_y_log10()+
  labs(y="Minimum molecular\nlesion duration (MMLD)",x="Phylogeny category")

#Display the difference between C>T lesion durations (shorter) & non C>T lesion durations (longer)
p.PVV.6.1<-PVV_mutations%>%
  # mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%
  # mutate(sub_cat=factor(sub_cat,levels=c("C>T","T>C","other")))%>%
  mutate(sub_cat=ifelse(Sub1=="C>T","C>T","other"))%>%
  filter(cat=="Adult_HSPC")%>%
  filter(lesion_duration<400)%>%
  ggplot(aes(x=lesion_duration,fill=sub_cat)) +
#  geom_histogram()+
  geom_density()+
  facet_grid(rows=vars(sub_cat),scales="free")+
  #scale_x_log10()+
  scale_fill_brewer(palette = "Set2")+
  my_theme +
  theme(legend.position = "none") +
  labs(title="Minimum molecular lesion duration")

p.PVV.6.2<-PVV_mutations%>%
  mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%
  mutate(sub_cat=factor(sub_cat,levels=c("C>T","T>C","other")))%>%
  filter(cat=="Adult_HSPC")%>%
  ggplot(aes(x=sub_cat,y=lesion_duration,col=sub_cat)) +
  geom_boxplot(size=0.5,outlier.shape = NA)+
  geom_jitter(size=0.5,width=0.2,alpha=0.3)+
  scale_y_log10()+
  scale_color_brewer(palette = "Set2")+
  my_theme +
  theme(legend.position = "none") +
  labs(x="PVV mutation type",y="Minimum lesion duration")

PVV_mutations%>%
  mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%
  mutate(sub_cat=factor(sub_cat,levels=c("C>T","T>C","other")))%>%
  filter(cat=="Adult_HSPC")%>%
  group_by(sub_cat)%>%
  summarise(n=n(),mean=mean(lesion_duration),median=median(lesion_duration))

#Test the distribution of the C>T MMLDs with the T>C and 'other' categories  
wilcox.test(x = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="C>T"&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="T>C"&cat=="Adult_HSPC")%>%pull(lesion_duration))

wilcox.test(x = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="C>T"&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="other"&cat=="Adult_HSPC")%>%pull(lesion_duration))

#Now test the distribution of adult HSPC removed MAVs with the PVVs as a whole, and the 'other' categories
wilcox.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%filter(cat=="Adult_HSPC")%>%pull(lesion_duration))

wilcox.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="C>T"&cat=="Adult_HSPC")%>%pull(lesion_duration))

ks.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="T>C"&cat=="Adult_HSPC")%>%pull(lesion_duration))

ks.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="other"&cat=="Adult_HSPC")%>%pull(lesion_duration))


#Display minimum lesion onset & earliest lesion repair bars by sample
p.PVV.7.1<-PVV_mutations%>%
  mutate(Sub1=mapply(FUN=function(x,trinuc_ref_py) {if(nchar(x)!=3) {
    return("INDEL")
  } else if(x=="C>T" & substr(trinuc_ref_py,2,3)=="CG"){
    return("C>T at CpG")
  } else if(x=="C>T" & substr(trinuc_ref_py,2,3)!="CG"){
    return("C>T other")
  } else {
    return(x)
  }
  }, x=PVV_mutations$Sub1,trinuc_ref_py=PVV_mutations$trinuc_ref_py,SIMPLIFY = T))%>%
  filter(Sample_ID%in%c("KX003_5_01","KX004_5_01","KX008_2_01"),!is.na(lesion_timing))%>%
  arrange(Sample_ID,lesion_timing)%>%
  mutate(order=1:sum(Type=="PVV"))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(xmin=(order-0.9),xmax=order,ymin=lesion_timing-2,ymax=lesion_repair_timing+2))+
  geom_rect(aes(fill=factor(Sub1))) +
  scale_y_continuous(limits=function(x){c(0,max(x[2],1200))})+
  scale_fill_brewer(palette="Set2")+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_grid(~Sample_ID,scales = "free",space="free")+
  labs(y="Molecular time",
       x="Individual PVVs, by minimum lesion onset",
       fill="Mutation class")

p.PVV.7.1a<-lapply(c("KX003_5_01","KX004_5_01","KX008_2_01"),function(Sample_ID) {
  tree<-read.tree(get_file_paths_and_project("EM",Sample_ID)$tree_file_path)
  internal_node_heights=nodeHeights(tree)[which(!tree$edge[,2]%in%1:length(tree$tip.label)),2]
  internal_node_heights<-internal_node_heights[!internal_node_heights<=50] #exclude developmental nodes
  return(data.frame(Sample_ID=Sample_ID,node_heights=internal_node_heights))
})%>%dplyr::bind_rows()%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=Sample_ID,y=node_heights))+
  geom_violin(fill="grey")+
  theme_classic()+
  my_theme


p.PVV.7.2<-PVV_mutations%>%
  mutate(Sub1=mapply(FUN=function(x,trinuc_ref_py) {if(nchar(x)!=3) {
    return("INDEL")
  } else if(x=="C>T" & substr(trinuc_ref_py,2,3)=="CG"){
    return("C>T at CpG")
  } else if(x=="C>T" & substr(trinuc_ref_py,2,3)!="CG"){
    return("C>T other")
  } else {
    return(x)
  }
  }, x=PVV_mutations$Sub1,trinuc_ref_py=PVV_mutations$trinuc_ref_py,SIMPLIFY = T))%>%
  filter(Sample_ID%in%c("PX001_2_01"),!is.na(lesion_timing))%>%
  arrange(Sample_ID,lesion_timing)%>%
  mutate(order=1:sum(Type=="PVV"))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(xmin=(order-0.9),xmax=order,ymin=lesion_timing-2,ymax=lesion_repair_timing+2))+
  geom_rect(aes(fill=factor(Sub1))) +
  scale_y_continuous(limits=function(x){c(0,max(x[2],3500))})+
  scale_fill_brewer(palette="Set2")+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_grid(~Sample_ID,scales = "free",space="free")+
  labs(y="Molecular time",
       x="Individual PVVs, by minimum lesion onset",
       fill="Mutation class")

p.PVV.7.2a<-lapply(c("PX001_2_01"),function(Sample_ID) {
  tree<-read.tree(get_file_paths_and_project("EM",Sample_ID)$tree_file_path)
  internal_node_heights=nodeHeights(tree)[which(!tree$edge[,2]%in%1:length(tree$tip.label)),2]
  internal_node_heights<-internal_node_heights[!internal_node_heights<=50] #exclude developmental nodes
  return(data.frame(Sample_ID=Sample_ID,node_heights=internal_node_heights))
})%>%dplyr::bind_rows()%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=Sample_ID,y=node_heights))+
  scale_y_continuous(limits=function(x){c(0,max(x[2],3500))})+
  geom_violin(fill="grey")+
  theme_classic()+
  my_theme
#ASSESS NUMBERS of PVVs per sample
n_PVV_summary<-PVV_mutations%>%
  filter(!is.na(lesion_node))%>%
  mutate(Sample_ID=factor(Sample_ID,levels=summary_table_df$Sample_ID))%>%
  group_by(Sample_ID)%>%
  summarise(n_PVV=n())

#Join this to the sample summary table
summary_table_df_PVV<-summary_table_df%>%
  left_join(n_PVV_summary)%>%
  mutate(n_PVV=replace_na(n_PVV,0))%>%
  mutate(smoking_status=factor(smoking_status,levels=c("unknown","Current_smoker","Ex_smoker","Never_smoker","Child")))

#PVV analysis by sample - comparisons with the simulated capture rate
capture_rate=read.delim("capture_rate_sim_poisson.tsv")
#Summary of the capture rates per sample - divided by data set
p.PVV.8.0<-left_join(summary_table_df_PVV,capture_rate,by="Sample_ID")%>%
  filter(cat!="Foetal_HSPC" & !is.na(capture_rate))%>%
  ggplot(aes(x=data_set,col=data_set,y=capture_rate,label=Sample_ID))+
  geom_jitter(width=0.1,height=0.1)+
  geom_text_repel(show.legend=F,fontface="italic",size=3)+
  scale_y_log10()+
  theme_classic()

#Capture rate vs n PVV
p.PVV.8.1<-left_join(summary_table_df_PVV,capture_rate,by="Sample_ID")%>%
  filter(cat!="Foetal_HSPC" & !is.na(capture_rate))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(data_set=ifelse(data_set=="MF"|Sample_ID%in%c("Pair28"),"Clonal\nhaemopoiesis",ifelse(data_set=="NW","MPN",ifelse(Sample_ID=="PX001","Chemo\nHSPC",ifelse(data_set=="KY","Bronchial\nepithelium","Adult\nHSPC")))))%>%
  ggplot(aes(x=capture_rate,col=data_set,y=n_PVV,label=Sample_ID))+
  geom_point(size=0.6,alpha=0.5)+
  geom_text_repel(show.legend=F,fontface="italic",size=2)+
  theme_classic()+
  my_theme+
  labs(x="Simulated PVV capture rate",y="Number of PVVs",col="Phylogeny type")
  

#Capture rate vs n PVV - excluding the chemo samples & those with large clonal expansions
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
p.PVV.8.2<-left_join(summary_table_df_PVV,capture_rate,by="Sample_ID")%>%
  filter(cat!="Foetal_HSPC" & data_set!="NW" & data_set!="MF"&!is.na(capture_rate) & !Sample_ID%in%c("Pair28","PX001_2_01"))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=capture_rate,col=cat,y=n_PVV,label=Sample_ID))+
  geom_point(size=0.6,alpha=0.5)+
  geom_smooth(method="lm",aes(x=capture_rate,y=n_PVV),col="black",inherit.aes=F)+
  geom_text_repel(show.legend=F,fontface="italic",size=2)+
  scale_color_brewer(palette = "Set2")+
  scale_x_continuous(labels = fancy_scientific)+
  theme_classic()+
  my_theme+
  labs(x="Simulated PVV capture rate",y="Number of PVVs",col="Data category")

lm.data<-left_join(summary_table_df_PVV,capture_rate,by="Sample_ID")%>%
  filter(cat!="Foetal_HSPC" & data_set!="NW" & data_set!="MF"&!is.na(capture_rate) & !Sample_ID%in%c("Pair28","PX001_2_01"))
capturerate.lm<-lm(n_PVV~capture_rate,data=lm.data)
summary(capturerate.lm)



if(resave_plots){
  #Save these plots
  ggsave(p.PVV.1.1,filename=paste0(plots_dir,"PVV_ASCAT_LOH.pdf"),width=2.7,height=2)
  ggsave(p.PVV.1.2,filename=paste0(plots_dir,"PVV_Phasable_mutations.pdf"),width=2.7,height=2)
  ggsave(p.PVV.1.3,filename=paste0(plots_dir,"PVV_reads_LOH.pdf"),width=2.7,height=2)
  ggsave(p.PVV.1.4,filename=paste0(plots_dir,"PVV_reads_LOH_ASCAT_pos_only.pdf"),width=1.5,height=2)
  ggsave(p.PVV.1.5,filename=paste0(plots_dir,"PVV_ASCAT_reads_LOH_only.pdf"),width=1.5,height=2)
  ggsave(p.PVV.3,filename=paste0(plots_dir,"PVV_phasing_outcomes_by_data.pdf"),width=2.7,height=2)
  ggsave(p.PVV.4,filename = paste0(plots_dir,"PVV_7profile_mutsig.pdf"),width = 12,height=5)
  ggsave(p.PVV.5.1,filename = paste0(plots_dir,"PVV_lesion_durations_by_sample.pdf"),width = 5,height=5)
  ggsave(p.PVV.5.2,filename = paste0(plots_dir,"PVV_lesion_durations_by_sample_boxplot.pdf"),width = 7,height=2)
  ggsave(p.PVV.5.3,filename = paste0(plots_dir,"PVV_lesion_durations_by_cat_violinplot.pdf"),width = 2.3,height=2)
  ggsave(p.PVV.6.1,filename = paste0(plots_dir,"PVV_TC_vs_CT_durations.pdf"),width = 2.3,height=2)
  ggsave(p.PVV.6.2,filename = paste0(plots_dir,"PVV_durations_by_type_boxplot.pdf"),width = 2,height=2)
  ggsave(p.PVV.7.1,filename = paste0(plots_dir,"PVV_timings_through_life.pdf"),width = 5,height=2)
  ggsave(p.PVV.7.1a,filename = paste0(plots_dir,"PVV_timings_through_life_nodedensity.pdf"),width = 2,height=2)
  ggsave(p.PVV.7.2,filename = paste0(plots_dir,"PVV_timings_through_life_PX001.pdf"),width = 4.5,height=2)
  ggsave(p.PVV.7.2a,filename = paste0(plots_dir,"PVV_timings_through_life_PX001_nodedensity.pdf"),width = 1,height=2)
  ggsave(p.PVV.8.0,filename = paste0(plots_dir,"PVV_capture_rate_by_sample_dataset.pdf"),width = 6,height=4)
  ggsave(p.PVV.8.1,filename = paste0(plots_dir,"PVV_nPVV_vs_capturerate_all.pdf"),width = 3.5,height=2)
  ggsave(p.PVV.8.2,filename = paste0(plots_dir,"PVV_nPVV_vs_capturerate_CH_excluded.pdf"),width = 3.5,height=2)
}

summary_table_df_comb<-summary_table_df_MAV%>%
  dplyr::select(Sample_ID,n_MAV)%>%
  right_join(summary_table_df_PVV,by="Sample_ID")%>%
  dplyr::select(Sample_ID,data_set,cat,n_sample,n_MAV,n_PVV)%>%
  gather(-Sample_ID,-data_set,-n_sample,-cat,key="Mutation_type",value="n")

summary_table_df_comb$cat<-factor(summary_table_df_comb$cat,levels=c("Foetal_HSPC","Adult_HSPC","Chemo_HSPC","Bronchial","Liver"))
sample_order<-summary_table_df_comb%>%mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%arrange(cat,n_sample)%>%pull(Sample_ID)%>%unique()

p.summary<-summary_table_df_comb%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Sample_ID=factor(Sample_ID,levels=sample_order))%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  mutate(cat=factor(cat,levels=c("Foetal\nHSPC","Adult\nHSPC","Chemo\nHSPC","Bronchial","Liver")))%>%
  mutate(Mutation_type=stringr::str_remove(Mutation_type,pattern = "n_"))%>%
  ggplot(aes(x=Sample_ID,y=n,fill=Mutation_type))+
  geom_bar(stat="identity",position="dodge")+
  facet_grid(~cat,drop = F,scales="free",space="free")+
  scale_y_continuous(limits=c(0,120))+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,size=4))+
  labs(y="Number detected",fill="Mutation\ntype")

ggsave(p.summary,filename=paste0(plots_dir,"Summary_nPVV_nMAV.pdf"),width=7,height=2)

#WRITE VCFS OF MUTATION SETS FOR SIGNATURE/ STRAND BIAS ANALYSIS
details=PVV_mutations%>%filter(cat=="Chemo_HSPC")%>%dplyr::select(Chrom,Pos,Ref,Alt1)%>%dplyr::rename("Alt"=Alt1)
PVV_vcf_path="Chemo_PVVs.vcf"
write.vcf(details,vcf_path = PVV_vcf_path,vcf_header_path = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt")

details=PVV_mutations%>%filter(cat=="Adult_HSPC" & Class=="PASS")%>%dplyr::select(Chrom,Pos,Ref,Alt1)%>%dplyr::rename("Alt"=Alt1)
PVV_vcf_path="Blood_PVVs.vcf"
write.vcf(details,vcf_path = PVV_vcf_path,vcf_header_path = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt")

#Write combined vcf of lung and liver MAVs
details=MAV_mutations%>%filter(Class!="FAIL"&(cat=="Liver"|cat=="Bronchial"))%>%dplyr::select(Chrom,Pos,Ref,Alt1)%>%dplyr::rename("Alt"=Alt1)%>%
  bind_rows(MAV_mutations%>%filter(Class!="FAIL"&(cat=="Liver"|cat=="Bronchial"))%>%dplyr::select(Chrom,Pos,Ref,Alt2)%>%dplyr::rename("Alt"=Alt2))
MAV_vcf_path="Lung_and_liver_MAVs.vcf"
write.vcf(details,vcf_path = MAV_vcf_path,vcf_header_path = "/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt")

##SUMMARISE MUTATION NUMBERS IN EACH SAMPLE##
#How many of each type of mutation do we have?
summary<-mutations%>%
  group_by(Type) %>%
  summarise(count=n())

summary_by_sample<-mutations%>%
  group_by(Type,Sample_ID) %>%
  summarise(count=n())

#How many MAVs are also PVVs
overlap<-mutations%>%
  filter(Type=="MAV")%>%
  dplyr::count(MAV_is_PVV)

#Draw Venn diagram of this overlap
font="Helvetica"
dev.off();plot.new()
draw.pairwise.venn(area1 =summary$count[summary$Type=="MAV"],
                   area2=summary$count[summary$Type=="PVV"],
                   cross.area = overlap$n[overlap$MAV_is_PVV=="Yes"],
                   category=c("MAVs","PVVs"),
                   fill=palette2,
                   alpha=0.5,
                   cat.pos = c(0,0),
                   cat.cex=1.5,
                   cat.fontfamily=c(rep(font,2)),
                   fontfamily = font
)
#Now visualize
#AS SEPARATE FACETS FOR MAVs/ PVVs
sample_order=mutations%>%dplyr::select(Sample_ID,cat)%>%filter(!duplicated(Sample_ID))%>%arrange(cat)%>%dplyr::pull(Sample_ID)

mutations%>% #and then also broken down by sample
  mutate(Sample_ID=factor(Sample_ID,levels=sample_order),Variant_type=as.factor(Type))%>%
  dplyr::count(Sample_ID,Variant_type)%>%
  tidyr::complete(Sample_ID,Variant_type,fill=list(n=0))%>%
  ggplot(aes(x=Sample_ID,y=n,fill=Variant_type)) +
  geom_bar(stat="identity",col="black",size=0.4) +
  facet_wrap(~Variant_type) +
  my_theme+
  #geom_text(aes(label=n),vjust=-0.5,size=2)+
  scale_y_continuous(breaks=seq(0,90,20),limits=c(0,90)) +
  scale_fill_manual(values=palette2)+
  labs(title="MAVs and PVVs by sample",x=NULL,y="Count") +
  theme(axis.text.x=element_text(size=7,angle = 90),legend.position = "none")

#AS DODGED BAR PLOT FOR MAVs/ PVVs
p1<-mutations%>% #and then also broken down by sample
  mutate(Sample_ID=factor(Sample_ID,levels=sample_order),Variant_type=as.factor(Type)) %>%
  dplyr::count(Sample_ID,Variant_type)%>%
  complete(Sample_ID,Variant_type,fill=list(n=0))%>%
  ggplot(aes(x=Sample_ID,y=n,fill=Variant_type)) +
  geom_bar(width=0.8,stat="identity",position = position_dodge(preserve = "single"),col="black",size=0.4) +
  my_theme+
  scale_y_continuous(breaks=seq(0,90,20),limits=c(0,90)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values=palette2)+
  theme(legend.position = "right")+
  labs(title="MAVs and PVVs by sample",
       x=NULL,
       y="Count",
       fill="Variant type")

mutations%>% #and then also broken down by sample
  mutate(Sample_ID=as.factor(Sample_ID),Variant_type=as.factor(Type)) %>%
  dplyr::count(cat,Variant_type)%>%
  complete(cat,Variant_type,fill=list(n=0))%>%
  ggplot(aes(x=cat,y=n,fill=Variant_type)) +
  geom_bar(width=0.8,stat="identity",position = position_dodge(preserve = "single"),col="black",size=0.4) +
  my_theme+
  scale_y_continuous(breaks=seq(0,150,20),limits=c(0,150)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values=palette2)+
  theme(legend.position = "right")+
  labs(title="MAVs and PVVs by sample",
       x=NULL,
       y="Count",
       fill="Variant type")

#2. Now view the poor fit 96-profile (with trinucleotide context)
mutations$mut_profile_1=paste0(substr(mutations$trinuc_ref_py,1,1),"[",mutations$Sub1,"]",substr(mutations$trinuc_ref_py,3,3))
mutations$mut_profile_2=paste0(substr(mutations$trinuc_ref_py,1,1),"[",mutations$Sub2,"]",substr(mutations$trinuc_ref_py,3,3))
freqs=table(mutations$mut_profile_1[mutations$Type=="PVV" & mutations$Sample_ID!="PX001"])
freqs_PX001=table(mutations$mut_profile_1[mutations$Type=="PVV" & mutations$Sample_ID=="PX001"])
freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
freqs_full_PX001 = freqs_PX001[full_vec]; freqs_full_PX001[is.na(freqs_full_PX001)] = 0; names(freqs_full_PX001) = full_vec
plot_96_profile(as.matrix(cbind(freqs_full,freqs_full_PX001)),ymax=0.3)

#3. Looks at MAVs
#Plot the 96-profile of MAVs with both variants (1 and 2) including separately
freqs=table(c(mutations$mut_profile_1[mutations$Type=="MAV"&mutations$Sample_ID!="PX001"],mutations$mut_profile_2[mutations$Type=="MAV"&mutations$Sample_ID!="PX001"]))
freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
plot_96_profile(as.matrix(freqs_full),ymax = 0.1)


##INCLUDE LESION DURATION/ TIMING
#Review the poor fit lesion durations

#Review the MAVs with lesion durations
mutations%>%
  filter(Type=="MAV" & lesion_duration>0)%>%
  ggplot(aes(x=Sample_ID,y=lesion_duration,col=Sample_ID)) +
  geom_jitter(width = 0.1,height=0,alpha=0.5) +
  #geom_violin()+
  facet_wrap(~cat,scales="free",ncol=4)+
  my_theme+
  theme(legend.position = "none")

#Review all MAVs or PVVs with lesion durations
p4.1<-MAV_mutations%>%
  filter(lesion_duration>0)%>%
  ggplot(aes(x=cat,y=lesion_duration,col=Type)) +
  geom_jitter(width = 0.1,height=1) +
  my_theme+
  scale_color_manual(values=palette2)+
  theme(axis.text.x=element_text(angle=0))+
  labs(title="Minimum lesion durations (molecular time)",
       x="",
       y="Molecular lesion duration")

p4.2<-PVV_mutations%>%
  filter(cat=="Adult_HSPC")%>%
  filter(lesion_duration>0)%>%
  ggplot(aes(x=Sample_ID,y=lesion_duration,col=Type)) +
  geom_jitter(width = 0.1,height=1) +
  my_theme+
  scale_color_manual(values=palette2)+
  theme(axis.text.x=element_text(angle=0))+
  labs(title="Minimum lesion durations (molecular time) - Adult blood",
       x="",
       y="Molecular lesion duration")

mutations%>%
  filter(Sample_ID=="PX001")%>%
  filter(lesion_duration>0)%>%
  ggplot(aes(x=Type,y=lesion_duration,col=Type)) +
  geom_jitter(width = 0.2,height=10,alpha=0.5) +
  #geom_violin()+
  my_theme+
  theme()+
  labs(title = "Lesion duration of MAVs and PVVs in individual PX001")

#Review numbers of cell divisions of each MAV/ PVV

p5<-mutations%>%
  #filter(Type=="MAV")%>%
  ggplot(aes(x=factor(no_of_cell_divisions),fill=Type))+
  geom_bar(col="black")+
  scale_fill_manual(values=palette2)+
  my_theme + theme(axis.text.x = element_text(angle=0)) +
  facet_wrap(~Type)+
  labs(title = "Minimum number of cell divisions of persistent DNA lesions (MAVs)",
       x="Minimum cell divisions")


#Show all multi-allelic variants in a single plot
multi_SNVs%>%
  filter(!is.na(lesion_timing))%>%
  arrange(Sample_ID,lesion_timing,Chrom_pos)%>%
  mutate(order=1:(length(Sample_ID)))%>%
  ggplot(aes(xmin=(order-0.9),xmax=order,ymin=lesion_timing,ymax=lesion_repair_timing+10))+
  geom_rect(aes(fill=factor(mut_combination))) +
  my_theme+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-10,NA))+
  facet_wrap(~Sample_ID,scales = "free")+
  labs(title = "Multi-allelic variants across all samples",
       x="Mutations",
       y="Molecular time",
       fill="Mutation type")

comb_part1=arrangeGrob(p2,p3,ncol=1)
comb_plot=arrangeGrob(p1,comb_part1,ncol=2,widths = c(5,3))
plot(comb_plot)
figures_dir="~/Documents/Lesion_segregation/Figures/"
ggsave(comb_plot,filename = paste0(figures_dir,"Figure_2.pdf"),device=cairo_pdf,width = 183,height = 100,units = "mm")
ggsave(p1,filename = paste0(figures_dir,"MAVs_and_PVVs_by_sample.pdf"),device=cairo_pdf,width = 183,height = 100,units = "mm")
ggsave(p2,filename = paste0(figures_dir,"PVVs_simple_signature.pdf"),device=cairo_pdf,width = 100,height = 50,units = "mm")
ggsave(p3,filename = paste0(figures_dir,"MAV_signatures_by_category.pdf"),device=cairo_pdf,width = 140,height = 50,units = "mm")
ggsave(arrangeGrob(p4.1,p4.2,ncol=1),filename = paste0(figures_dir,"Molecular_lesion_durations.pdf"),device=cairo_pdf,width = 140,height = 100,units = "mm")
ggsave(p5,filename = paste0(figures_dir,"Minimum_cell_divisions.pdf"),device=cairo_pdf,width = 140,height = 70,units = "mm")

#Review of lesion nodes on the tree
pdf("Affected_MAV_nodes_on_tree.pdf",width=15,height=10)
sapply(MAV_mutations%>%filter(Class!="FAIL")%>%pull(Sample_ID)%>%unique(), function(sampleID) {
  data_set=MAV_mutations%>%filter(Sample_ID==sampleID)%>%pull(data_set)%>%unique()
  sample_info=get_file_paths_and_project(dataset = data_set,Sample_ID=sampleID)
  tree=read.tree(sample_info$tree_file_path)
  node_vec=MAV_mutations%>%filter(Sample_ID==sampleID & Class!="FAIL")%>%pull(lesion_node)

  tree=plot_tree(tree,cex.label = 0) #Here need to make sure that have the correct
  add_annotation(tree,
                 details,
                 matrices = NULL,
                 annot_function = function(tree=tree,details=details,matrices=matries,node,node_vec=node_vec){
                   if(sum(node_vec==node)>0){
                     n_node=base::table(node_vec)[as.character(node)]
                     radius=sqrt(n_node)
                     info=get_edge_info(tree,details,node)
                     plotrix::draw.circle(info$x,info$yb,radius = radius,col="red",density=50)
                   }
                 },
                 node_vec=node_vec)
})
dev.off()


#CONFIRM THE PHYLOGENIES BY LOOKING AT THE MUTATIONS CONFIRMING THE STRUCTURE
confirm_phylogeny=lapply(data_sets, function(dataset) {
  dataset="EM"
  phasing_output_dir=paste0("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/phasing_output/",dataset)
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
    pdf(paste0(output_dir,"/",sample,"_confirm_PVVs.pdf"),width=20,height=6)
    sample_out=lapply(1:nrow(data_set_sample_PVVs),function(i) {
      print(i)
      Chrom=str_split(data_set_sample_PVVs$Chrom_pos[i],pattern = "-",simplify=T)[,1]
      mut=data_set_sample_PVVs$mut_ref1[i]
      print(mut)
      
      Pos=as.numeric(str_split(data_set_sample_PVVs$Chrom_pos[i],pattern = "-",simplify=T)[,2])
      Ref<-data_set_sample_PVVs$Ref[i]
      Alt<-data_set_sample_PVVs$Alt1[i]
      
      if(is.na(data_set_sample_MAVs$lesion_node[i])) {stop(return("No lesion node"))}
      
      #Here to identify the mixed subclade & otherwise confirm the mutations that cause this phylogeny
      mixed_subclade=get_mixed_subclades(mut1 = mut,lesion_node = data_set_sample_PVVs$lesion_node[i],tree=tree_curr,matrices=matrices)
      print(mixed_subclade)
      if(any(mixed_subclade=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")) {stop(return("More than one mixed subclade"))}
      
      clade_muts=details$mut_ref[details$node==mixed_subclade]
      print(paste("There are",length(clade_muts),"mutations confirming this clade"))
      
      if(length(clade_muts)>5) {
        to_visualize=sample(clade_muts,5)
      } else {
        to_visualize=clade_muts
      }
      
      sapply(to_visualize,function(confirm_mut) {
        tree_curr=plot_tree(tree_curr,cex.label = 0)
        add_annotation(tree_curr,
                       details=details,
                       matrices=matrices,
                       annot_function = confirm_PVV_phylogeny,
                       mut=confirm_mut,
                       PVV_mut=mut,
                       lesion_node=data_set_sample_PVVs$lesion_node[i],
                       lwd=2,
                       cex=0.4)
      })
      return(length(clade_muts))
    })
    dev.off()
    return(sample_out)
  })
  return(data_set_out)
})

