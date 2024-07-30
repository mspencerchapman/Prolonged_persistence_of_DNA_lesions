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
                axis.text = element_text(size = 8),
                axis.title = element_text(size=10),
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
# Review the read-based method LOH results of the PVVs ####
#========================================#
PVV_mutations<-mutations%>%
  filter(Type=="PVV")%>%
  mutate(PVV_pos_clade_phasing=factor(PVV_pos_clade_phasing,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

myLOHColors<-c(rev(RColorBrewer::brewer.pal(3,"Accent")),"grey")
names(myLOHColors)<-c("LOH","Unable to confirm","No LOH","Homozygous in all clades")

## Generate Extended Data Figure 5a ------
p.LOH.readbased<-PVV_mutations%>%
  mutate(ASCAT_result=factor(ASCAT_result))%>%
  filter(!is.na(lesion_node))%>%
  mutate(basic_result=factor(basic_result))%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  dplyr::count(basic_result,cat)%>%
  ggplot(aes(x=cat,y=n,fill=basic_result))+
  geom_bar(stat="identity",position="stack",col=NA,linewidth=0.2)+
  scale_fill_manual(name = "LOH result",values=myLOHColors)+
  theme_classic()+
  my_theme+
  labs(x="Data type",
       y="Number of mutations")+
  theme(legend.position="inside",
        legend.position.inside=c(0.6,0.7),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust=+0.5,vjust = +0.5),
        legend.key.size = unit(3,"mm"))

ggsave(p.LOH.readbased,filename=paste0(plots_dir,"ExtDatFig5a.pdf"),width=2,height=2)

## Generate Extended Data Figure 5b ------
p.LOHvsASCAT<-PVV_mutations%>%
  filter(!is.na(lesion_node) & !is.na(basic_result) &!ASCAT_result=="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected" &
           max_neg_clade_depth>=13)%>%
  dplyr::select(ASCAT_result,basic_result)%>%
  mutate(ASCAT_result=ifelse(ASCAT_result=="Homozygous in all clades","No LOH",ASCAT_result))%>%
  table()%>%
  as.data.frame()%>%
  mutate(ASCAT_result=factor(ASCAT_result,levels=c("LOH","No LOH","Unable to confirm")),Read_based_result=factor(basic_result,levels=c("LOH","No LOH","Unable to confirm")))%>%
  ggplot(aes(x=ASCAT_result,y=Read_based_result,fill=Freq))+
  geom_tile()+
  scale_fill_gradient(low="lightgrey",high="orange")+
  geom_text(aes(label = Freq), vjust = 1)+
  coord_flip()+
  my_theme+
  labs(x="Read-based result",y="ASCAT result")+
  theme(legend.key.size = unit(5,"mm"))

ggsave(p.LOHvsASCAT,filename=paste0(plots_dir,"ExtDatFig5b.pdf"),width=3,height=2)

## Generate Extended Data Figure 5c ------
rename_phasing=c("Matching phasing","Opposite phasing","Unable to confirm")
names(rename_phasing)=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")

myPhasingColors <- c("#478EB0","#A35563","grey")
names(myPhasingColors) <- rename_phasing

p.PVV.3.1<-PVV_mutations%>%
  filter(!is.na(lesion_node) & !is.na(basic_result) &!ASCAT_result=="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected" &
           max_neg_clade_depth>=13)%>%
  filter(lesion_duration>0)%>%
  mutate(PVV_pos_clade_phasing=factor(rename_phasing[PVV_pos_clade_phasing],levels=rename_phasing))%>%
  mutate(no_of_cell_divisions=no_of_cell_divisions-1)%>%
  ggplot(aes(x=no_of_cell_divisions,fill=PVV_pos_clade_phasing))+
  geom_bar(stat="count",position="stack",width=0.7)+
  scale_x_continuous(breaks=seq(2,10,1),limits=c(1.5,12))+
  geom_vline(xintercept = 6.5,linetype=2)+
  scale_fill_manual(values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  theme_classic()+
  my_theme+
  labs(x="Minimum cell divisions",y="Number of mutations")+
  theme(legend.position="inside",legend.position.inside = c(0.75,0.7),legend.key.size = unit(3,"mm"),legend.title = element_blank())

ggsave(p.PVV.3.1,filename=paste0(plots_dir,"ExtDatFig5c.pdf"),width=2.5,height=2)
