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
# DEFINE CUSTOM FUNCTIONS FOR SCRIPT ####
#========================================#

#The function to plot MAV profiles from the data
plot_MAV_96profile=function(multi_SNV_table,cat_types=NULL,relative=F,return_mat=F,compressed=F,ymax=0.1,breaks=NULL) {
  require(ggplot2)
  require(tidyr)
  require(dplyr)
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
    
    mut_freqs=pivot_wider(mut_table,names_from = cat,values_from = n)
    mut_mat=as.matrix(mut_freqs[,cat_types]); rownames(mut_mat)=paste(mut_freqs$mut_combination,mut_freqs$trinuc_ref_py,sep=",")
    
  } else {
    mut_table=multi_SNV_table%>%
      mutate(mut_combination=as.factor(mut_combination))%>%
      mutate(trinuc_ref_py=as.factor(trinuc_ref_py))%>%
      dplyr::count(trinuc_ref_py,mut_combination)%>%
      tidyr::complete(trinuc_ref_py,mut_combination,fill=list(n=0))%>%
      filter(substr(trinuc_ref_py,2,2)==substr(mut_combination,1,1))
    
    mut_vec=mut_table$n; names(mut_vec)=paste(mut_table$mut_combination,mut_table$trinuc_ref_py,sep=",")
    mut_mat=as.matrix(mut_vec,ncol=1); colnames(mut_mat)<-"MAVs"
  }
  if(relative) {
    mut_freqs=pivot_wider(mut_table,names_from = cat,values_from = n)
    mut_freqs[,3:ncol(mut_freqs)]<-apply(mut_freqs[,3:ncol(mut_freqs)],2,function(x) return(x/sum(x)))
    mut_table=gather(mut_freqs,-trinuc_ref_py,-mut_combination,key="cat",value="n")%>%
      mutate(cat=factor(cat,levels=cat_types))
  }
  
  if(return_mat) {stop(return(mut_mat))}
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

#The function to plot MAV profiles from the data
plot_MAV_from_mat=function(mut_mat,cat_types=NULL,compressed=F,ymax=0.1,breaks=NULL) {
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  
  COLORS6=c("#66C2A5" ,"#FC8D62" ,"#8DA0CB" ,"#E78AC3" ,"#A6D854" ,"#FFD92F") #Different colours to the normal 96-profile sig to clarify different meaning
  
  if(is.null(cat_types)){
    cat_types<-colnames(mut_mat)
  }
  
  #THE PLOT
  plot<-as.data.frame(mut_mat)%>%
    tibble::rownames_to_column(var = "var")%>%
    tidyr::separate(var,into=c("mut_combination","context"),sep = ",")%>%
    mutate(context=paste0(substr(context,1,1),".",substr(context,3,3)))%>%
    tidyr::gather(-mut_combination,-context,key="cat",value="n")%>%
    ggplot(aes(x = context, y = n, fill = mut_combination))+
    geom_bar(stat = "identity",width=0.8)+
    scale_fill_manual(values = COLORS6)+
    facet_grid(cat~mut_combination,scales = "free_y",space = "free")+
    ylab("Relative contribution")+
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

#The function to infer MAV profiles from signatures
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
  mut_table$mut_combination=paste(substr(mut_table$mut1,3,5),substr(mut_table$mut2,3,5),sep=" + ")
  
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
    geom_bar(stat = "identity", colour = "black",size = 0.2,width=0.6)+
    scale_fill_manual(values = COLORS6)+
    facet_grid(Sample~mut_combination)+
    ylab(ifelse(relative,"Relative contribution","Frequency"))+
    coord_cartesian(ylim = c(0, ymax))+
    scale_y_continuous(breaks = seq(0, ymax, ifelse(relative,0.1,1))) + 
    guides(fill = FALSE)+
    theme_bw()+
    theme(rect=element_rect(size=0.3),
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
      theme(rect=element_rect(size=0.3),
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
# DEFINE CUSTOM FUNCTION FOR SCRIPT ####
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
# CREATE VCFs OF REPRESENTATIVE SIGNATURES ####
#========================================#

# Write VCFs of the lung/ liver mutations ----
bronchial_samples<-mutations%>%filter(Type=="MAV" & data_set=="KY")%>%group_by(Sample_ID)%>%summarise(n=n())%>%filter(n>5)%>%pull(Sample_ID)
liver_samples<-mutations%>%filter(Type=="MAV" & data_set=="SN")%>%pull(Sample_ID)

temp=Map(Sample_IDs=list(bronchial_samples,liver_samples),data_set=c("KY","SN"),tissue=c("bronchial","liver"),function(Sample_IDs,data_set,tissue) {
  
  full_mutation_set<-lapply(Sample_IDs,function(Sample_ID) {
    info=get_file_paths_and_project(dataset=data_set,Sample_ID=Sample_ID,input_data_dir = lesion_seg_input_dir)
    load(info$filtered_muts_path)
    return(filtered_muts$COMB_mats.tree.build$mat%>%dplyr::select(mut_ref,Chrom,Pos,Ref,Alt))
  })%>%dplyr::bind_rows()
  
  #Reduce the size down to 20000 mutations - sufficient to have high confidence profiles
  if(nrow(full_mutation_set)>2e4) {
    set.seed(101)
    full_mutation_set<-full_mutation_set[sample(1:nrow(full_mutation_set),size = 2e4),]
  }
  
  #Write the VCF file
  write.vcf(details = full_mutation_set,vcf_path = paste0(root_dir,"Data/VCFs/",tissue,"_muts.vcf"),vcf_header_path = vcf_header_path)
})

#Reimport these VCFs and the Normal blood ('Pair11') and chemo blood ('PX001')----
VCF_dir=paste0(root_dir,"Data/VCFs/")
VCF_file_paths=paste0(VCF_dir,c("bronchial_muts","liver_muts","Pair11","PX001"),".vcf")
vcfs=read_vcfs_as_granges(vcf_files = VCF_file_paths,
                          sample_names = c("Bronchial","Liver","Adult_HSPC","Chemo_HSPC"),genome=ref_genome)
mut_mat=mut_matrix(vcfs,ref_genome = ref_genome)
mut_mat_props=apply(mut_mat,2,function(x) {x/sum(x)})
plot_96_profile(mut_mat_props,ymax = 0.07)

#Calculate anticipated MAV signatures from each tissue ----
mut_mat_anticipated_ind<-plot_anticipated_MAV_sig_from_96profile(mut_mat_props = mut_mat_props,return_mat = T)
mut_mat_anticipated_ind<-mut_mat_anticipated_ind%>%
  tidyr::unite(col="mut_type",mut_combination,trinuc_ref_py,sep = ",")%>%
  tibble::column_to_rownames(var="mut_type")

## Create table of MAVs to plot signatures - these should all be SNV combinations (no indels) ----
multi_SNVs<-mutations%>%
  filter(Type=="MAV"&Class=="FAIL")%>%
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

MAV_mut_mat=plot_MAV_96profile(multi_SNVs,cat_types = c("Bronchial","Liver","Adult_HSPC","Chemo_HSPC"),return_mat = T)
MAV_mut_mat<-apply(MAV_mut_mat,2, function(x) {x/sum(x)})
colnames(MAV_mut_mat)<-paste(colnames(MAV_mut_mat),"Unrelated_MAVs",sep="_")

## Combine the 'anticipated' with 'actual' unrelated MAV profiles
MAV_comb_mat<-cbind(mut_mat_anticipated_ind,MAV_mut_mat[rownames(mut_mat_anticipated_ind),])

## Generate Extended Data Fig. 4a-d ----
## Plot the anticipated MAV signatures besides the actual MAV signatures
temp=Map(tissue=c("Bronchial","Liver","Adult_HSPC","Chemo_HSPC"),Fig=c("a","b","c","d"),function(tissue,Fig) {
  cat(tissue,sep="\n")
  p.MAV_vs_expected<-plot_MAV_from_mat(MAV_comb_mat[,c(tissue,paste(tissue,"Unrelated_MAVs",sep = "_"))],ymax = 0.1,compressed = T)
  ggsave(p.MAV_vs_expected,filename = paste0(plots_dir,"ExtDatFig4",Fig,".pdf"),width = 7,height=3)
  
})
