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

mut_examples=c("15-29890184-G-A","15-31430023-G-A","16-86458492-C-T")

PVV_examples<-mutations%>%filter(mut_ref1%in%mut_examples)

## Generate Figure 4a ------
pdf(paste0(plots_dir,"ExtDatFig1e.pdf"),width=7,height=5)
par(mfrow=c(1,3))
temp=sapply(1:nrow(PVV_examples),function(i) {
  cat(i,sep="\n")
  
  SampleID=PVV_examples$Sample_ID[i]
  data_set=PVV_examples$data_set[i]
  LN=PVV_examples$lesion_node[i]
  LRN=PVV_examples$lesion_repair_node[i]
  mut_ref=PVV_examples$mut_ref1[i]
  
  #Import the info
  sample_info=get_file_paths_and_project(dataset=data_set,SampleID,input_data_dir=paste0(root_dir,"Data/input_data/"))
  tree=read.tree(sample_info$tree_file_path)
  load(sample_info$filtered_muts_path)
  details=filtered_muts$COMB_mats.tree.build$mat
  NV=filtered_muts$COMB_mats.tree.build$NV
  NR=filtered_muts$COMB_mats.tree.build$NR
  
  #Extract the subtree of the PVV
  sub_tree<-extract.clade(phy=tree,node=LN)
  if(!"node"%in% colnames(details)) {stop(return(NULL))}
  LRN_ancestors<-get_ancestral_nodes(LRN,tree$edge)
  LN_ancestors<-get_ancestral_nodes(LN,tree$edge)
  LP_nodes<-LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
  confirmatory_mutations<-details%>%filter(node%in%LP_nodes)%>%pull(mut_ref)%>%sample(size=15)
  combined_muts=c(mut_ref,confirmatory_mutations)
  print(filtered_muts$COMB_mats.tree.build$NR[combined_muts,sub_tree$tip.label]) #confirm the coverage
  
  ##Plot the original PVV on the mpboot tre
  hm=create_mutation_heatmap(mut_refs1=mut_ref,mut_refs2=confirmatory_mutations,matrices=filtered_muts$COMB_mats.tree.build,tree=sub_tree)
  
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=1.5,title = mut_ref)
  add_mut_heatmap_PVV(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.5)
  
})
dev.off()