#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr")
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

palette2=c("#E30613","#1D71B8")
palette6=c("#72e5ef", "#074d65", "#3d99ce", "#c257d3", "#d1add5", "#61356e")
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

#========================================#
# GENERATE SUMMARY OF MUTATION NUMBERS ####
#========================================#

## Assess the NUMBERS of MAVs in different samples----
summary_table_df_MAV=mutations%>%
  filter(Type=="MAV" & Class!="FAIL")%>%
  mutate(Sample_ID=factor(Sample_ID,levels=summary_table_df$Sample_ID))%>%
  group_by(Sample_ID)%>%
  summarise(n_MAV=n())%>%
  right_join(summary_table_df)%>%
  mutate(n_MAV=replace_na(n_MAV,0))%>%
  mutate(smoking_status=factor(smoking_status,levels=c("unknown","Current_smoker","Ex_smoker","Never_smoker","Child")))


## Assess NUMBERS of PVVs per sample ----
summary_table_df_PVV<-mutations%>%
  filter(Type=="PVV" & !is.na(lesion_node) & Class!="FAIL")%>%
  mutate(Sample_ID=factor(Sample_ID,levels=summary_table_df$Sample_ID))%>%
  group_by(Sample_ID)%>%
  summarise(n_PVV=n())%>%
  right_join(summary_table_df)%>%
  mutate(n_PVV=replace_na(n_PVV,0))%>%
  mutate(smoking_status=factor(smoking_status,levels=c("unknown","Current_smoker","Ex_smoker","Never_smoker","Child")))

## Combine summary tables ----
summary_table_df_comb<-summary_table_df_MAV%>%
  dplyr::select(Sample_ID,n_MAV)%>%
  right_join(summary_table_df_PVV,by="Sample_ID")%>%
  dplyr::select(Sample_ID,data_set,cat,n_sample,n_MAV,n_PVV)%>%View()
  gather(-Sample_ID,-data_set,-n_sample,-cat,key="Mutation_type",value="n")

sample_order<-summary_table_df_comb%>%mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%arrange(cat,n_sample)%>%pull(Sample_ID)%>%unique()

rename_cat_vec=c("Adult HSPC","C","F","Bronchial","Liver (MAVs only)")
names(rename_cat_vec)=c("Adult_HSPC","Chemo_HSPC","Foetal_HSPC","Bronchial","Liver")

## Generate Extended Data Fig. 9a ------
p.summary<-summary_table_df_comb%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Sample_ID=factor(Sample_ID,levels=sample_order))%>%
  mutate(cat=factor(rename_cat_vec[cat],levels=rename_cat_vec))%>%
  mutate(Mutation_type=stringr::str_remove(Mutation_type,pattern = "n_"))%>%
  ggplot(aes(x=Sample_ID,y=n,fill=Mutation_type))+
  geom_bar(stat="identity",position="dodge")+
  facet_grid(~cat,drop = F,scales="free",space="free")+
  scale_y_continuous(limits=c(0,120))+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,size=4),axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "inside",legend.position.inside=c(0.1,0.7),strip.background = element_rect(fill=NA))+
  labs(y="Number detected")

ggsave(p.summary,filename=paste0(plots_dir,"ExtDatFig9a.pdf"),width=7,height=2)


#========================================#
# SIMULATED CAPTURE RATE ####
#========================================#

#This is calculated by simulations in the script 'Tree_structure_PVV_sensitivity_2.R'
capture_rate=read.delim(paste0(root_dir,"Data/simulation_results/capture_rate_sim_poisson.tsv"))

## Generate Extended Data Fig. 9b ------
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
  labs(x="Simulated PVV capture rate",y="Number of PVVs",col="")+
  theme(legend.position = "inside",legend.position.inside = c(0.5,0.9))

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

## Generate Extended Data Fig. 9c ------
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
  labs(x="Simulated PVV capture rate",y="Number of PVVs",col="")+
  theme(legend.position="inside",legend.position.inside=c(0.2,0.8))

lm.data<-left_join(summary_table_df_PVV,capture_rate,by="Sample_ID")%>%
  filter(cat!="Foetal_HSPC" & data_set!="NW" & data_set!="MF"&!is.na(capture_rate) & !Sample_ID%in%c("Pair28","PX001_2_01"))
capturerate.lm<-lm(n_PVV~capture_rate,data=lm.data)
summary(capturerate.lm)


#========================================#
# LESION SIMULATOR ####
#========================================#

#This has a custom function to simulate 'types' of detected lesions resulting from a persistent lesion, depending on
#1. the specific pairing probabilities of the lesion
#2. the number of lesions persisting across 1, or more than 1 node

calculated_PVV_to_MAV_ratio=function(Pairing_probabilities,
                                     Two_to_three_replication_ratio) {
  
  #Define the different variables
  Pairing_possibilities=c("G*/C","G*/T","G*/A","G*/G")
  Correct_pairing="G*/C"
  Possible_outcomes=c("No mutation","Mutation","Simple MAV","Separated MAV","Separated MAV (triallelic)","PVV")
  Recognized_persistent_lesion_outcomes=c("Simple MAV","Separated MAV","Separated MAV (triallelic)","PVV")
  
  ## CREATE THE 'TWO REPLICATION' DATA FRAME
  Two_replications_df=expand.grid(Pairing_possibilities,Pairing_possibilities)
  colnames(Two_replications_df)=c("Replication_1","Replication_2")
  Two_replications_df$Result<-apply(Two_replications_df,1,function(vec) {
    if(vec[1]==vec[2]&vec[1]==Correct_pairing) {
      return("No mutation")
    } else if(vec[1]==vec[2]&vec[1]!=Correct_pairing) {
      return("Mutation")
    } else if(vec[1]==Correct_pairing|vec[2]==Correct_pairing) {
      return("Mutation")
    } else if(vec[1]!=Correct_pairing & vec[2]!=Correct_pairing & vec[1]!=vec[2]) {
      return("Simple MAV")
    }
  })
  
  Two_replications_df$Probability<-apply(Two_replications_df,1,function(vec) {
    prod(Pairing_probabilities[vec[1:2]])
  })
  
  ## CREATE THE 'THREE REPLICATION' DATA FRAME
  Three_replications_df=data.frame("Replication_1"=rep(Pairing_possibilities,each=16),
                                   "Replication_2"=rep(Two_replications_df$Replication_1,times=4),
                                   "Replication_3"=rep(Two_replications_df$Replication_2,times=4))
  
  Three_replications_df$Result<-apply(Three_replications_df,1,function(vec) {
    if(sum(vec==Correct_pairing)==3) {
      return("No mutation")
    } else if(vec[1]==Correct_pairing & vec[2]==vec[3] & vec[2]!=Correct_pairing) {
      return("Mutation")
    } else if(sum(vec==Correct_pairing)==2) {
      return("Mutation")
    } else if(vec[1]==Correct_pairing & vec[2]!=vec[3]) {
      return("Simple MAV")
    } else if(sum(vec==Correct_pairing)==1 & length(unique(vec))==2 & vec[2]!=vec[3]){
      return("PVV")
    } else if(sum(vec==Correct_pairing)==1 & length(unique(vec))==3 & vec[2]!=vec[3]){
      return("Separated MAV")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==3){
      return("Separated MAV (triallelic)")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==1){
      return("Mutation")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==2 & vec[2]!=vec[3]){
      return("Separated MAV")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==2 & vec[2]==vec[3]){
      return("Simple MAV")
    } else {
      return("Undefined")
    }
  })
  
  Three_replications_df$Probability<-apply(Three_replications_df,1,function(vec) {
    prod(Pairing_probabilities[vec[1:3]])
  })
  
  ##Create outcome summaries
  Two_replications_summary<-sapply(Possible_outcomes,function(outcome) {Two_replications_df%>%filter(Result==outcome)%>%pull(Probability)%>%sum()})
  Three_replications_summary<-sapply(Possible_outcomes,function(outcome) {Three_replications_df%>%filter(Result==outcome)%>%pull(Probability)%>%sum()})
  
  outcomes_summary=data.frame(Outcome=Possible_outcomes,
                              Two_replications=Two_replications_summary,
                              Three_replications=Three_replications_summary,
                              Totals=sapply(Possible_outcomes,function(outcome) {Two_to_three_replication_ratio*Two_replications_summary[outcome] + Three_replications_summary[outcome]})
  )
  
  outcomes_summary$Totals_normalized=outcomes_summary$Totals/sum(outcomes_summary$Totals)
  
  outcomes_summary$Recognized_as_persistent_lesions=sapply(Possible_outcomes,function(outcome) {
    if(outcome%in%Recognized_persistent_lesion_outcomes) {
      return(outcomes_summary%>%filter(Outcome==outcome)%>%pull(Totals))
    } else {
      return(0)
    }
  })
  
  outcomes_summary$Recognized_as_persistent_lesions_normalized=outcomes_summary$Recognized_as_persistent_lesions/sum(outcomes_summary$Recognized_as_persistent_lesions)
  
  MAV_total=sum(outcomes_summary$Recognized_as_persistent_lesions_normalized[grepl("MAV",outcomes_summary$Outcome)])
  PVV_total=sum(outcomes_summary$Recognized_as_persistent_lesions_normalized[grepl("PVV",outcomes_summary$Outcome)])
  
  PVV_to_MAV_ratio=PVV_total/MAV_total
  return(PVV_to_MAV_ratio)
}

## Try the function out with this toy data ----
Pairing_probabilities=c("G*/C"=0.49,
                        "G*/T"=0.49,
                        "G*/A"=0.005,
                        "G*/G"=0.005)


calculated_PVV_to_MAV_ratio(Pairing_probabilities=Pairing_probabilities,
                            Two_to_three_replication_ratio=3)


## Can express this more simply as a 'transition to transversion' ratio ----
#Assumes that the two tranversions & two transitions equally likely (a simplification)
Transition_to_transversion_ratio=100
Pairing_probabilities_unnormalized=c("G*/C"=Transition_to_transversion_ratio,
                                     "G*/T"=Transition_to_transversion_ratio,
                                     "G*/A"=1,
                                     "G*/G"=1)

Pairing_probabilities=Pairing_probabilities_unnormalized/sum(Pairing_probabilities_unnormalized)
calculated_PVV_to_MAV_ratio(Pairing_probabilities=Pairing_probabilities,
                            Two_to_three_replication_ratio=3)

## Apply function across a range of 'transition to tranversion' ratios ----
all_outcomes_df<-lapply(c(1,10,50),function(Two_to_three_replication_ratio) {
  outcomes_df=data.frame(Transition_to_transversion_ratio=10^seq(-3,3,0.1))
  outcomes_df$PVV_to_MAV_ratio=sapply(outcomes_df$Transition_to_transversion_ratio,function(Transition_to_transversion_ratio) {
    Pairing_probabilities_unnormalized=c("G*/C"=Transition_to_transversion_ratio,
                                         "G*/T"=Transition_to_transversion_ratio,
                                         "G*/A"=1,
                                         "G*/G"=1)
    
    Pairing_probabilities=Pairing_probabilities_unnormalized/sum(Pairing_probabilities_unnormalized)
    
    calculated_PVV_to_MAV_ratio(Pairing_probabilities=Pairing_probabilities,
                                Two_to_three_replication_ratio=Two_to_three_replication_ratio)
  })
  return(outcomes_df%>%mutate(Two_to_three_replication_ratio=Two_to_three_replication_ratio))
})%>%dplyr::bind_rows()

## Generate Extended Data Fig. 9d ------
breaks=c(0.01,0.1,1,10,100,1000)
PVV_to_MAV_ratio<-all_outcomes_df%>%
  ggplot(aes(x=Transition_to_transversion_ratio,y=PVV_to_MAV_ratio,col=factor(Two_to_three_replication_ratio)))+
  geom_line()+
  theme_classic()+
  scale_x_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  scale_y_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  geom_hline(yintercept=5,linetype=2)+
  geom_vline(xintercept=240,linetype=2)+
  labs(x="Transition to transversion ratio of TLS",
       y="PVV:MAV ratio of resulting\n(recognized) persistent lesions",
       col="Ratio of lesions persisting across only 1 node,\ncompared to those persisting across >=2")+
  my_theme

ggsave(filename=paste0(plots_dir,"ExtDatFig9d.pdf"),plot=PVV_to_MAV_ratio,width = 5,height=3)
