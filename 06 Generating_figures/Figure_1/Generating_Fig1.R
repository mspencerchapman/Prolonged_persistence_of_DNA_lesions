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
# PLOT THE THREE TREES ####
#========================================#
## Generate Figure 1g ----
info=get_file_paths_and_project(dataset="MSC_fetal",Sample_ID="18pcw",input_data_dir = paste0(root_dir,"Data/input_data"))
tree=read.tree(info$tree_file)
load(info$filtered_muts_path)
pdf(file=paste0(plots_dir,"Fig1g.pdf"),width = 5,height = 2)
tree=plot_tree(tree=tree,cex.label=0,lwd=0.5)
temp=add_annotation(tree,
               details=filtered_muts$COMB_mats.tree.build$mat,
               matrices=filtered_muts$COMB_mats.tree.build,
               annot_function = plot_MAV_mut,
               mut1=mutations%>%filter(Sample_ID=="18pcw")%>%pull(mut_ref1)%>%.[1],
               mut2=mutations%>%filter(Sample_ID=="18pcw")%>%pull(mut_ref2)%>%.[1],
               lwd=1,
               cex=0)
dev.off()

## Generate Figure 1h ----
info=get_file_paths_and_project(dataset="MSC_BMT",Sample_ID="Pair28",input_data_dir = paste0(root_dir,"Data/input_data"))
tree=read.tree(info$tree_file)
load(info$filtered_muts_path)
pdf(file=paste0(plots_dir,"Fig1h.pdf"),width = 5,height = 2)
tree=plot_tree(tree=tree,cex.label=0,lwd=0.5)
temp=add_annotation(tree,
                    details=filtered_muts$COMB_mats.tree.build$mat,
                    matrices=filtered_muts$COMB_mats.tree.build,
                    annot_function = plot_MAV_mut,
                    mut1=mutations%>%filter(Sample_ID=="Pair28")%>%pull(mut_ref1)%>%.[5],
                    mut2=mutations%>%filter(Sample_ID=="Pair28")%>%pull(mut_ref2)%>%.[5],
                    lwd=1,
                    cex=0)
dev.off()

## Generate Figure 1i ----
info=get_file_paths_and_project(dataset="MSC_fetal",Sample_ID="8pcw",input_data_dir = paste0(root_dir,"Data/input_data"))
tree=read.tree(info$tree_file)
load(info$filtered_muts_path)
pdf(file=paste0(plots_dir,"Fig1i.pdf"),width = 5,height = 2)
tree=plot_tree(tree=tree,cex.label=0,lwd=0.5)
temp=add_annotation(tree,
                    details=filtered_muts$COMB_mats.tree.build$mat,
                    matrices=filtered_muts$COMB_mats.tree.build,
                    annot_function = plot_MAV_mut,
                    mut1=mutations%>%filter(Sample_ID=="8pcw")%>%pull(mut_ref1),
                    lwd=1,
                    cex=0)
dev.off()
