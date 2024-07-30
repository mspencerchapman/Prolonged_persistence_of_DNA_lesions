#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr","ggrepel")
bioconductor_packages=c("GenomicRanges","IRanges","Rsamtools","MutationalPatterns","BSgenome","TxDb.Hsapiens.UCSC.hg19.knownGene",ref_genome)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
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

sig_theme=theme(rect=element_rect(linewidth =0.3),
                axis.title.y = element_text(size = 7,vjust = 1),
                axis.text.y = element_text(size = 6),
                axis.title.x = element_text(size = 7),
                axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
                strip.text.x = element_text(size = 6),
                strip.text.y = element_text(size = 6),
                legend.text = element_text(size=5),
                legend.key.size = unit(0.25,"cm"),
                panel.grid.major.x = element_blank(),
                panel.spacing.x = unit(0.5,"lines"))

condensed_sig_theme= theme(rect=element_rect(linewidth=0.3),
                           axis.title.y = element_text(size = 7,vjust = 1),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(), 
                           axis.title.x = element_text(size = 7),
                           axis.text.x = element_blank(),
                           strip.text.x = element_text(size = 5),
                           strip.text.y = element_text(size = 5),
                           legend.text = element_text(size=5),
                           panel.grid.major.x = element_blank(),
                           panel.spacing.x = unit(0.2,"lines"))

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
# Analyse the MAV mutations ####
#========================================#
MAV_mutations<-mutations%>%
  filter(Type=="MAV")%>%
  mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

#Generate summary table of numbers of MAV mutations ----
n_MAV_summary_table=MAV_mutations%>%
  filter(Class!="FAIL")%>%
  mutate(Sample_ID=factor(Sample_ID,levels=summary_table_df$Sample_ID))%>%
  group_by(Sample_ID)%>%
  summarise(n_MAV=n())

#Join this to the sample summary table ----
summary_table_df_MAV<-summary_table_df%>%
  left_join(n_MAV_summary_table)%>%
  mutate(n_MAV=replace_na(n_MAV,0))%>%
  mutate(smoking_status=factor(smoking_status,levels=c("unknown","Current_smoker","Ex_smoker","Never_smoker","Child")))%>%
  mutate(MAVs_per_PDN=n_MAV/post_development_nodes)

## Generate Figure 3a ------
MAVs_per_PDN<-summary_table_df_MAV%>%
  filter(cat!="Foetal_HSPC")%>%
  mutate(cat=factor(cat,levels=c("Adult_HSPC","Liver","Bronchial","Chemo_HSPC")))%>%
  ggplot(aes(x=cat,y=MAVs_per_PDN,col=cat,shape=smoking_status))+
  geom_boxplot(aes(x=cat,y=MAVs_per_PDN),inherit.aes = F,outlier.shape = NA,data=summary_table_df_MAV%>%
                 filter(!cat%in%c("Foetal_HSPC","Chemo_HSPC"))%>%
                 mutate(cat=factor(cat,levels=c("Adult_HSPC","Liver","Bronchial"))))+
  geom_jitter(width = 0.15,height=0,alpha=0.8)+
  scale_color_brewer(palette = "Set1",guide="none")+
  guides(shape=guide_legend(position="inside",keywidth=0.1,keyheight=0.3,default.unit="cm"))+
  scale_x_discrete(labels=function(x) {stringr::str_replace(x,"_","\n")})+
  theme_classic()+
  my_theme+
  theme(legend.position.inside =c(0.2,0.9),plot.margin = margin(t=8,r = 3,l=3,b=3,unit="mm"),legend.background=element_rect(fill=NA),axis.title.x = element_blank())+
  labs(y="Number of MAVs per\npost-development node",
       col="Sample type",
       shape="Smoking status")

ggsave(MAVs_per_PDN,filename=paste0(plots_dir,"Fig3a.pdf"),width=3.5,height=2.5)

#Do significance testing for differences in 'MAVs_per_node' statistic in different groups
df_for_kruskal_testing<-summary_table_df_MAV%>%
  filter(cat!="Foetal_HSPC" & post_development_nodes>0)%>%
  mutate(cat=factor(cat,levels=c("Adult_HSPC","Liver","Bronchial","Chemo_HSPC")))

kruskal.test(x = df_for_kruskal_testing$MAVs_per_PDN,
             g=df_for_kruskal_testing$cat)


#Do significance testing for differences in 'MAVs_per_node' statistic for SMOKERS vs NON-SMOKERS
smoking_df_for_kruskal_testing=df_for_kruskal_testing%>%
  filter(cat=="Bronchial")%>%
  mutate(smoking_status=ifelse(smoking_status%in%c("Never_smoker","Child"),"Never_smoker","Smoking_history"))%>%
  left_join(sample_ref%>%dplyr::select(Sample_ID,Pack_years))%>%
  mutate(Pack_years=replace_na(Pack_years,0))%>%
  filter(Pack_years==0|Pack_years>=30) #Include only 'never smokers' or individuals with at least 30 pack years

kruskal.test(x=smoking_df_for_kruskal_testing$MAVs_per_PDN,
             g=smoking_df_for_kruskal_testing$smoking_status)


#========================================#
# Analyse MAV signatures ####
#========================================#

## Create table of MAVs to plot signatures - these should all be SNV combinations (no indels) ----
multi_SNVs<-MAV_mutations%>%
  filter(Class!="FAIL")%>%
  filter(cat!="Bronchial"|Class!="removed")%>% #Remove the "removed" bronchial mutations
  filter(Mut_type1=="SNV"&Mut_type2=="SNV")

## Create consistent order for the multi-allelic mutations ----
#Transitions first (if either is), or else transversion to A first
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

## Create a "mut_combination" variable ----
multi_SNVs=cbind(multi_SNVs,multi_SNVs_ordered)
multi_SNVs$mut_combination=paste(multi_SNVs$sub1,multi_SNVs$sub2,sep=" + ")

#Now plot MAV signatures as 96-profile plots
#Design function to plot a bespoke 96-profile for 2 x SNV combinations at a single site
plot_MAV_96profile=function(multi_SNV_table,cat_types=NULL,relative=F,return_vec=F,compressed=F,ymax=0.1,breaks=NULL) {
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
    scale_y_continuous(breaks = breaks) + 
    guides(fill = "none")+
    theme_bw()+
    theme(axis.title.y = element_text(size = 8,vjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          panel.grid.major = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.x = unit(0.5,"lines"))
    
  if(compressed){
    plot+theme(rect=element_rect(linewidth=0.3),
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

## Generate Figure 3b ------
MAV_sig_all<-plot_MAV_96profile(multi_SNV_table = multi_SNVs,cat_types = c("Bronchial","Liver","Adult_HSPC","Chemo_HSPC"),breaks=seq(0,10,2))

ggsave(MAV_sig_all,filename=paste0(plots_dir,"Fig3b.pdf"),width=7,height=4)

#========================================#
# Analyse the PVV mutations ####
#========================================#

PVV_mutations<-mutations%>%filter(Type=="PVV")%>%
  mutate(Class=ifelse(is.na(lesion_node),"FAIL",Class))

## Extract the adult HSPC fail mutations for comparatory signature analysis----
FAIL_PVVs<-PVV_mutations%>%
  filter(is.na(lesion_node) & cat=="Adult_HSPC")%>%
  group_by(cat)%>%
  pull(mut_ref1)
write.vcf(FAIL_PVVs,vcf_path = paste0(root_dir,"Data/VCFs/Blood_fail_PVVs.vcf"),vcf_header_path = vcf_header_path)


PVV_mutations<-PVV_mutations%>%filter(Class!="FAIL")

## Write VCFs of the chemotherapy & blood PVVs ----
details=PVV_mutations%>%filter(cat=="Chemo_HSPC" & Class=="PASS")%>%dplyr::select(Chrom,Pos,Ref,Alt1)%>%dplyr::rename("Alt"=Alt1)
ChemoPVV_vcf_path=paste0(root_dir,"Data/VCFs/Chemo_PVVs.vcf")
write.vcf(details,vcf_path = ChemoPVV_vcf_path,vcf_header_path = vcf_header_path)

details=PVV_mutations%>%filter(cat=="Adult_HSPC" & Class=="PASS")%>%dplyr::select(Chrom,Pos,Ref,Alt1)%>%dplyr::rename("Alt"=Alt1)
BloodPVV_vcf_path=paste0(root_dir,"Data/VCFs/Blood_PVVs.vcf")
write.vcf(details,vcf_path = BloodPVV_vcf_path,vcf_header_path = vcf_header_path)

## Read back in SNVs from VCFs and plot the standard 96 profile matrix----
vcfs=read_vcfs_as_granges(vcf_files = c(BloodPVV_vcf_path,ChemoPVV_vcf_path),sample_names = c("Blood PVVs","Chemo PVVs"),genome=ref_genome,type="snv",predefined_dbs_mbs=T)

## Generate Figure 3c ------
mut_mat=mut_matrix(vcfs,ref_genome = ref_genome)
  p.PVV.1.1<-plot_96_profile(mut_mat[,c("Blood PVVs","Chemo PVVs")],ymax=0.32)+sig_theme
ggsave(p.PVV.1.1,filename = paste0(plots_dir,"Fig3c.pdf"),width=7,height=2.5)

## Generate Figure 3d ------
### Transcription strand bias ----
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand_trans <- mut_strand(vcfs[[1]], genes_hg19, mode = "transcription")
mut_mat_s_trans <- mut_matrix_stranded(vcfs[[1]], ref_genome, genes_hg19,mode = "transcription")
strand_counts_trans <- strand_occurrences(mut_mat_s_trans)
strand_bias_trans <- strand_bias_test(strand_counts_trans)
print(strand_bias_trans)
Blood_PVV_trans_bias<-plot_192_profile_single_context(mut_mat_s_trans[,1,drop=F],select_context="C>T",col="#DE1C14",ymax=0.3)+sig_theme+theme(strip.text.x=element_text(size=9))
ggsave(Blood_PVV_trans_bias,filename=paste0(plots_dir,"Fig3d.pdf"),width=2.5,height=2)

## Generate Figure 3e ------
#Write the complete set of mutations for PX001 (then used as a comparison)
load(get_file_paths_and_project(dataset="EM",Sample_ID = "PX001_2_01",input_data_dir = paste0(root_dir,"Data/input_data/"))$filtered_muts_path)
details<-filtered_muts$COMB_mats.tree.build$mat
PX001_all_VCF_path=paste0(root_dir,"Data/VCFs/PX001.vcf")
write.vcf(details,vcf_path = PX001_all_VCF_path,vcf_header_path = vcf_header_path)

#Read in indels from VCFs and plot the indel profile
indel_vcfs=read_vcfs_as_granges(vcf_files = c(ChemoPVV_vcf_path,PX001_all_VCF_path),sample_names = c("Chemo PVVs","Overall chemo sig."),genome=ref_genome,type="indel")
indel_vcfs<-get_indel_context(indel_vcfs,ref_genome=ref_genome)
indel_counts <- count_indel_contexts(indel_vcfs)
p.indelPVV.1.1<-plot_indel_contexts(indel_counts[,c("Chemo PVVs","Overall chemo sig.")])+sig_theme+theme(panel.spacing.x = unit(0.1,"lines"),legend.title = element_blank())+labs(x="Sequence context (count of repeat units)",y="Number of indels")
ggsave(p.indelPVV.1.1,filename = paste0(plots_dir,"Fig3e.pdf"),width=7,height=2.5)



