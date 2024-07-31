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
genome_file=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/reference_files/genome.fa","/nfs/cancer_ref02/human/GRCh37d5/genome.fa") ##Set to the location of GRCh37 genome file
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
# DEFINE CUSTOM FUNCTION FOR SCRIPT ####
#========================================#

plot_anticipated_MAV_sig_from_96profile=function(mut_mat_props,return_mat=F,relative=T,compressed=F,triplet_freqs=triplet_freqs_all) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  #COLORS6=c("#2EBAED" ,"#000000" ,"#DE1C14" ,"#D4D2D2" ,"#ADCC54" ,"#F0D0CE")
  COLORS6=c("#66C2A5" ,"#FC8D62" ,"#8DA0CB" ,"#E78AC3" ,"#A6D854" ,"#FFD92F")
  sig_names=colnames(mut_mat_props)
  
  mut_table=data.frame(mut1=rep(NA,96))
  mut_contexts=rownames(mut_mat_props)
  mut_table$mut1=c(mut_contexts[1:16],mut_contexts[33:48],mut_contexts[33:48],mut_contexts[49:64],mut_contexts[65:80],mut_contexts[65:80])
  mut_table$mut2=c(mut_contexts[17:32],mut_contexts[1:16],mut_contexts[17:32],mut_contexts[81:96],mut_contexts[49:64],mut_contexts[81:96])
  
  #Add columns of the trinuc_ref_py and mut_comb so that can plug into MAV_96_profile_sig function
  mut_table$trinuc_ref_py=paste0(substr(mut_table$mut1,1,1),substr(mut_table$mut1,3,3),substr(mut_table$mut1,7,7))
  mut_table$mut_combination=paste(substr(mut_table$mut1,3,5),substr(mut_table$mut2,3,5),sep="_")
  
  #Can input different signatures here
  res_list=lapply(1:ncol(mut_mat_props),function(i) {
    mut1_prop=mut_mat_props[mut_table$mut1,i]; mut2_prop=mut_mat_props[mut_table$mut2,i]
    n=mut1_prop*mut2_prop/triplet_freqs[mut_table$trinuc_ref_py] #Work out the product of the proportions, and divide by the trinucleotide abundances
    n=n/sum(n) #Normalize
    return(data.frame(Sample=colnames(mut_mat_props)[i],mut_combination=mut_table$mut_combination,trinuc_ref_py=mut_table$trinuc_ref_py,n=n))
  })
  res_df=dplyr::bind_rows(res_list)
  
  mut_freqs=pivot_wider(res_df,names_from = Sample,values_from = n)
  if(return_mat){stop(return(mut_freqs))}
  ymax=max(res_df$n)
  
  #THE PLOT
  plot<-res_df%>%
    mutate(Sample=factor(Sample,levels=sig_names),context=paste0(substr(trinuc_ref_py,1,1),".",substr(trinuc_ref_py,3,3)))%>%
    ggplot(aes(x = context, y = n, fill = mut_combination))+
    geom_bar(stat = "identity", colour = "black",linewidth = 0.2,width=0.6)+
    scale_fill_manual(values = COLORS6)+
    facet_grid(Sample~mut_combination)+
    ylab(ifelse(relative,"Relative contribution","Frequency"))+
    coord_cartesian(ylim = c(0, ymax))+
    scale_y_continuous(breaks = seq(0, ymax, ifelse(relative,0.1,1))) + 
    guides(fill = FALSE)+
    theme_bw()+
    theme(rect=element_rect(linewidth=0.3),
          axis.title.y = element_text(size = 8,vjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0.5,"lines"))
  
  if(compressed){
    plot<-plot+scale_y_continuous(breaks = NULL)+geom_bar(stat = "identity", width=0.8)+
      theme(rect=element_rect(linewidth=0.3),
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

#========================================#
# Get the trinucleotide frequencies across the genome ####
#========================================#

#Load the trinucleotide frequencies of the human genome - needed for correction of signatures
Hg_19_strings<-readDNAStringSet(genome_file, format="fasta",
                                nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
triplet_freqs=Biostrings::trinucleotideFrequency(Hg_19_strings)
triplet_freqs_all=colSums(triplet_freqs)
triplet_freqs_all=triplet_freqs_all/sum(triplet_freqs_all) #Normalize the frequencies

#========================================#
# Import blood mutations - represent 'normal blood' sig ####
#========================================#
Blood_sig_vcf_path=paste0(root_dir,"Data/VCFs/Pair11.vcf")

if(!file.exists(Blood_sig_vcf_path)) {
  load(get_file_paths_and_project(dataset="MSC_BMT",Sample_ID="Pair11",input_data_dir=lesion_seg_input_dir)$filtered_muts_path)
  write.vcf(filtered_muts$COMB_mats.tree.build$mat,vcf_path = Blood_sig_vcf_path,vcf_header_path=vcf_header_path)
}

vcfs=read_vcfs_as_granges(vcf_files = Blood_sig_vcf_path,sample_names = "Blood_signature",genome=ref_genome)
mut_mat=mut_matrix(vcfs,ref_genome = ref_genome)
plot_96_profile(mut_mat,ymax=0.07)
mut_mat_props=mut_mat/sum(mut_mat)

#========================================#
# Calculated MAV signature from LUNG signatures ####
#========================================#

#Anticipated MAV signature for the de novo extracted signatures in Kenichi's manuscript
lung_sigs<-read.delim(paste0(root_dir,"Data/reference_files/known_signatures/ExtendedFigure5_hdp_96probabilities.txt"))
lung_sigs<-lung_sigs[,-c(1,2)]
rownames(lung_sigs)<-rownames(mut_mat)

## Generate Extended Data Fig. 8b ------
#Create plot of the original extracted signatures
p.Lung.1.1<-plot_96_profile(lung_sigs[,c("SBS1","SBS4","SBS5","SBS16","SBS18","Sig.A","Sig.B"),drop=F],ymax=0.12)+geom_bar(stat = "identity", width=0.8)+condensed_sig_theme+labs(x="Sequence context")

#Plot the inferred MAV signatures from these
p.MAV.1.1<-plot_anticipated_MAV_sig_from_96profile(mut_mat_props = lung_sigs[,c("SBS1","SBS4","SBS5","SBS16","SBS18","Sig.A","Sig.B"),drop=F],compressed=T)+theme(axis.title.y = element_blank())+labs(x="Sequence context")

## Combine
lung_comb_plot<-gridExtra::arrangeGrob(grobs=list(p.Lung.1.1,p.MAV.1.1),nrow=1)
plot(lung_comb_plot)

ggsave(lung_comb_plot,filename = paste0(plots_dir,"ExtDatFig8b.pdf"),width = 7,height=3.5)


#========================================#
# Calculated MAV signature from LIVER signatures ####
#========================================#

#Anticipated MAV signature for the de novo extracted signatures in Stan's manuscript
liver_sigs<-read.delim(paste0(root_dir,"Data/reference_files/known_signatures/Liversigs_96profile.txt"))

## Generate Extended Data Fig. 8c ------
#Create plot of the original extracted signatures
p.Liver.1.1<-plot_96_profile(liver_sigs[,c(3,5)],ymax=0.25)+geom_bar(stat = "identity", width=0.8)+condensed_sig_theme+labs(x="Sequence context")

#Plot the inferred MAV signatures from these
p.MAV.1.1a<-plot_anticipated_MAV_sig_from_96profile(mut_mat_props = liver_sigs[,c(3,5),drop=F],compressed=T)+theme(axis.title.y = element_blank())+labs(x="Sequence context")

## Combine
liver_comb_plot<-gridExtra::arrangeGrob(grobs=list(p.Liver.1.1,p.MAV.1.1a),nrow=1)
plot(liver_comb_plot)

ggsave(liver_comb_plot,filename = paste0(plots_dir,"ExtDatFig8c.pdf"),width = 7,height=2.5)

#Can also review the anticipated signatures from all the other components if desired - though not in manuscript.
p.MAV.1.1b<-plot_anticipated_MAV_sig_from_96profile(mut_mat_props = liver_sigs[,11:20,drop=F],compressed=T)+theme(axis.title.y = element_blank())
p.MAV.1.1c<-plot_anticipated_MAV_sig_from_96profile(mut_mat_props = liver_sigs[,21:30,drop=F],compressed=T)+theme(axis.title.y = element_blank())

#========================================#
# Calculated MAV signature from BLOOD signatures ####
#========================================#

#Anticipated MAV signature for the signatures from Emily's chemo patient PX001
EM_sigs<-t(read.delim(paste0(root_dir,"Data/mutsig_extraction/EM_extracted_signatures.tsv")))
rownames(EM_sigs)<-rownames(mut_mat);colnames(EM_sigs)<-c("0","BM_signature","Chemo 1","Chemo 2")

## Generate Extended Data Fig. 8d ------
## Create plot of the original extracted signatures
p.Blood.1.1<-plot_96_profile(EM_sigs[,c("BM_signature","Chemo 1","Chemo 2")],ymax=0.1)+geom_bar(stat = "identity", width=0.8)+condensed_sig_theme+labs(x="Sequence context")

## Plot the inferred MAV signatures from these
p.MAV.1.2<-plot_anticipated_MAV_sig_from_96profile(mut_mat_props = EM_sigs[,c("BM_signature","Chemo 1","Chemo 2")],triplet_freqs=triplet_freqs_all,compressed=T)+theme(axis.title.y = element_blank())+labs(x="Sequence context")

## Combine
blood_comb_plot<-gridExtra::arrangeGrob(grobs=list(p.Blood.1.1,p.MAV.1.2),nrow=1)
plot(blood_comb_plot)

ggsave(blood_comb_plot,filename = paste0(plots_dir,"ExtDatFig8d.pdf"),width = 7,height=2.5)



