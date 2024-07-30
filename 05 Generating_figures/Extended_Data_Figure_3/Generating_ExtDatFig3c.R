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

Phasing_MAV_file_path=paste0(data_dir,"Phasing_results_MAVs_all")
Phasing_PVV_file_path=paste0(data_dir,"Phasing_results_PVVs_all")
ASCAT_PVV_file_path=paste0(data_dir,"ASCAT_LOH_analysis_PVVs_all")
ASCAT_MAV_file_path=paste0(data_dir,"ASCAT_LOH_analysis_MAVs_all")
SN_phasing_file_path=paste0(data_dir,"Phasing_results_MAVs_SN.csv")

lesion_seg_input_dir=paste0(data_dir,"input_data/")
lesion_seg_output_dir=output_dir

mutations_filt_file=paste0(data_dir,"mutations_filtered.tsv")
summary_df_file=paste0(data_dir,"summary_df.tsv")


mutations<-readr::read_delim(mutations_filt_file,col_types = "cccccccciicccciiiiiiccccccnccccccciccccci")
summary_table_df<-read.delim(summary_df_file)
sample_ref=read.csv(paste0(data_dir,"metadata/Individual_ref.csv"))


#========================================#
# Set the root directory and read in the necessary files ####
#========================================#
MAV_mutations<-mutations%>%
  filter(Type=="MAV")%>%
  mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

rename_phasing=c("Matching phasing","Opposite phasing","Unable to confirm")
names(rename_phasing)=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")

myPhasingColors <- c("#478EB0","#A35563","grey")
names(myPhasingColors) <- rename_phasing

# Generate Extended Data Fig. 3c ----
p.MAV.0.1<-MAV_mutations%>%
  filter(Original_Class=="removed")%>%
  mutate(phasing_summary=factor(phasing_summary))%>%
  mutate(n_neg_cat=ifelse(n_neg>1,">=2","<=1"))%>%
  group_by(phasing_summary)%>%
  dplyr::count(phasing_summary,n_neg_cat)%>%
  mutate(phasing_summary=factor(rename_phasing[phasing_summary],levels=rename_phasing))%>%
  ggplot(aes(x=n_neg_cat,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="stack",col="black",linewidth=NA)+
  scale_fill_manual(name = "MAV phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  guides(fill=guide_legend(position="right",ncol=1,override.aes = list(linewidth=0)))+
  theme_classic()+
  my_theme+
  theme(legend.title=element_blank(),legend.key.size = unit(3,"mm"))+
  labs(x="Number of negative\nsubclades",y="Number of mutations")

ggsave(p.MAV.0.1,filename = paste0(plots_dir,"ExtDatFig3c.pdf"),width=2.5,height=2)
