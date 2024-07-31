#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr","ggrepel")
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
data_dir=paste0(root_dir,"Data/")
plots_dir=paste0(root_dir,"plots/")
output_dir=paste0(root_dir,"/output/")
source(paste0(root_dir,"/Data/Prolonged_persistence_functions.R"))
ref_table=read.csv(paste0(root_dir,"/Data/metadata/Individual_ref.csv"))
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa" ##Set to the location of GRCh37 genome file
vcf_header_path=paste0(root_dir,"Data/reference_files/vcfHeader.txt")

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
# Visualize the embryonic mutation in KX001 ####
#========================================#

#Plot the early embryonic MAV from KX001
sample_info=get_file_paths_and_project("EM","KX001_4_01",input_data_dir=paste0(root_dir,"Data/input_data/"))
tree=read.tree(sample_info$tree_file_path)
load(sample_info$filtered_muts_path)
details=filtered_muts$COMB_mats.tree.build$mat
NV=filtered_muts$COMB_mats.tree.build$NV
NR=filtered_muts$COMB_mats.tree.build$NR

## Generate Figure 4a ------
mut1="X-147677255-A-G"
mut2="X-147677255-A-T"
tree.squash=squash_tree(tree,cut_off=15)
pdf(paste0(plots_dir,"Fig4a.pdf"),width=7,height=2.5)
tree.squash=plot_tree(tree.squash,cex.label=0,scale=5)
temp=add_annotation(tree.squash,
               details=details,
               matrices=list(NV=NV,NR=NR),
               annot_function = plot_MAV_mut,
               mut1=mut1,
               mut2=mut2,
               lwd=2,
               cex=0.4)
dev.off()

#========================================#
# Analyse the PVV mutations ####
#========================================#
PVV_mutations<-mutations%>%filter(Type=="PVV")%>%
  mutate(Class=ifelse(is.na(lesion_node),"FAIL",Class))

PVV_mutations<-PVV_mutations%>%filter(Class!="FAIL")

## Generate Figure 4b ------
#Display minimum lesion onset & earliest lesion repair bars by sample
#This is the main plot, but doesn't include the node density information - need to manually add this in via illustrator

sub_categories=PVV_mutations%>%
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
  arrange(desc(Mut_type1),Sub1)%>%
  pull(Sub1)%>%unique()
sub_cols=RColorBrewer::brewer.pal(7,"Set2")[c(1,8,2,3,5,6,7,4)]
names(sub_cols)<-sub_categories

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
  mutate(Sub1=factor(Sub1,levels=sub_categories))%>%
  filter(Sample_ID%in%c("KX003_5_01","KX004_5_01","KX008_2_01"),!is.na(lesion_timing))%>%
  arrange(Sample_ID,lesion_timing)%>%
  mutate(order=1:sum(Type=="PVV"))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(xmin=(order-0.9),xmax=order,ymin=lesion_timing-2,ymax=lesion_repair_timing+2))+
  geom_rect(aes(fill=Sub1),show.legend = T) +
  scale_y_continuous(breaks=seq(0,1250,250),limits=c(0,1500))+
  scale_fill_manual(values=sub_cols,drop=F)+
  guides(fill=guide_legend(title=element_blank(),ncol=3,position="inside",))+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position.inside = c(0.5,0.9),strip.background = element_blank(),strip.text.x = element_blank())+
  facet_grid(~Sample_ID,scales = "free",space="free")+
  labs(y="Molecular time (muts)",
       x="Individual PVVs")

p.PVV.7.1a<-lapply(c("KX003_5_01","KX004_5_01","KX008_2_01"),function(Sample_ID) {
  tree<-read.tree(get_file_paths_and_project("EM",Sample_ID,input_data_dir=paste0(root_dir,"Data/input_data/"))$tree_file_path)
  internal_node_heights=nodeHeights(tree)[which(!tree$edge[,2]%in%1:length(tree$tip.label)),2]
  internal_node_heights<-internal_node_heights[!internal_node_heights<=50] #exclude developmental nodes
  return(data.frame(Sample_ID=Sample_ID,node_heights=internal_node_heights))
})%>%dplyr::bind_rows()%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=Sample_ID,y=node_heights))+
  geom_violin(fill="grey")+
  scale_y_continuous(breaks=seq(0,1250,250),limits=c(0,1500))+
  theme_classic()+
  my_theme+
  theme(axis.title = element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  labs(y="Molecular time (muts)")
  
PVV_timing_aged_blood<-gridExtra::arrangeGrob(grobs=list(p.PVV.7.1,p.PVV.7.1a),widths = c(4,1))
ggsave(PVV_timing_aged_blood,filename = paste0(plots_dir,"Fig4b.pdf"),width = 7,height=2)
  
## Generate Figure 4c ------
#This is the same but for the chemotherapy tree
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
  mutate(Sub1=factor(Sub1,levels=sub_categories))%>%
  filter(Sample_ID%in%c("PX001_2_01"),!is.na(lesion_timing))%>%
  arrange(Sample_ID,lesion_timing)%>%
  mutate(order=1:sum(Type=="PVV"))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(xmin=(order-0.9),xmax=order,ymin=lesion_timing-2,ymax=lesion_repair_timing+2))+
  geom_rect(aes(fill=Sub1),show.legend = T) +
  scale_y_continuous(limits=c(0,3500))+
  scale_fill_manual(values=sub_cols,drop=F)+
  guides(fill=guide_legend(title=element_blank(),ncol=3,position="right",drop=F))+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())+
  facet_grid(~Sample_ID,scales = "free",space="free")+
  labs(y="Molecular time (muts)",
       x="Individual PVVs",
       fill="Mutation class")

p.PVV.7.2a<-lapply(c("PX001_2_01"),function(Sample_ID) {
  tree<-read.tree(get_file_paths_and_project("EM",Sample_ID,input_data_dir=paste0(root_dir,"Data/input_data/"))$tree_file_path)
  internal_node_heights=nodeHeights(tree)[which(!tree$edge[,2]%in%1:length(tree$tip.label)),2]
  internal_node_heights<-internal_node_heights[!internal_node_heights<=50] #exclude developmental nodes
  return(data.frame(Sample_ID=Sample_ID,node_heights=internal_node_heights))
})%>%dplyr::bind_rows()%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  ggplot(aes(x=Sample_ID,y=node_heights))+
  scale_y_continuous(limits=c(0,3500))+
  geom_violin(fill="grey")+
  theme_classic()+
  my_theme+
  labs(y="Molecular time (muts)")+
  theme(axis.title = element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())

PVV_timing_chemo_blood<-gridExtra::arrangeGrob(grobs=list(p.PVV.7.2,p.PVV.7.2a),widths = c(5,1))
ggsave(PVV_timing_chemo_blood,filename = paste0(plots_dir,"Fig4c.pdf"),width = 7,height=2)

## Generate Figure 4d ------
p.PVV.5.3<-PVV_mutations%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  ggplot(aes(x=cat,y=lesion_duration,fill=cat)) +
  geom_violin()+
  geom_jitter(col="grey",width=0.3,alpha=0.4,size=0.6)+
  my_theme +
  scale_fill_manual(values=c( "#E41A1C","#4DAF4A" ,"#984EA3" ,"#FF7F00"))+
  theme(legend.position = "none",axis.text.x=element_text(size=7)) +
  scale_y_log10()+
  labs(y="Minimum molecular\nlesion duration (MMLD)",x="Phylogeny category")

## Generate Figure 4e ------
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
  labs(x="PVV mutation type",y="Minimum molecular\nlesion duration (MMLD)")

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


#========================================#
# Analyse the MAV mutations ####
#========================================#
MAV_mutations<-mutations%>%
  filter(Type=="MAV")%>%
  mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

## Generate Figure 4f ------
p.MAV.10.1<-MAV_mutations%>%
  filter(Class!="FAIL")%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  filter(!is.na(lesion_duration)&Class=="removed"&lesion_duration!=0)%>%
  ggplot(aes(x=cat,y=lesion_duration,fill=cat)) +
  geom_violin(width=0.7)+
  geom_jitter(col="grey",width=0.2,size=1,alpha=0.5)+
  my_theme +
  scale_fill_manual(values=c( "#E41A1C","#4DAF4A" ,"#984EA3" ,"#FF7F00"))+
  theme(legend.position = "none",axis.text.x=element_text(size=7),axis.title.x = element_blank(),title = element_text(size=7,face="bold")) +
  scale_y_log10()+
  labs(y="Minimum molecular\nlesion duration (MMLD)",title = "Separated MAV durations")

ggsave(p.MAV.10.1,filename=paste0(plots_dir,"Fig4f.pdf"),width=2.5,height=2)


#Now test the distribution of adult HSPC removed MAVs with the PVVs as a whole, and the 'other' categories
wilcox.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%filter(cat=="Adult_HSPC")%>%pull(lesion_duration))

wilcox.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="C>T"&cat=="Adult_HSPC")%>%pull(lesion_duration))

ks.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="T>C"&cat=="Adult_HSPC")%>%pull(lesion_duration))

ks.test(x = MAV_mutations%>%filter(Class=="removed"&lesion_duration!=0&!is.na(lesion_duration)&cat=="Adult_HSPC")%>%pull(lesion_duration),
        y = PVV_mutations%>%mutate(sub_cat=ifelse(Sub1=="C>T","C>T",ifelse(Sub1=="T>C","T>C","other")))%>%filter(sub_cat=="other"&cat=="Adult_HSPC")%>%pull(lesion_duration))



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

