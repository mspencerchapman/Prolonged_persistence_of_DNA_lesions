#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","phangorn","MASS","tidyr","deconstructSigs")
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
# Define custom functions for script ####
#========================================#

mutlist_to_96_contexts=function(mutlist,build=37,genomeFile=lustre_ref){
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  if(is.null(mutlist$SampleID)){
    samples="this_sample"
    mutlist$SampleID="this_sample"
  } else {
    samples=unique(mutlist$SampleID)
  }
  
  trinuc_mut_mat=matrix(0,ncol=96,nrow=length(samples))
  if(build==37) {
    chromosomes=c(1:22,"X","Y")
  } else if(build==38) {
    chromosomes=paste0("chr",c(1:22,"X","Y"))
  } else {
    stop(cat("Specified genome build must be either '37' or '38'"))
  }
  for (n in 1:length(samples)){
    s=samples[n]
    mutations=as.data.frame(mutlist[mutlist$SampleID==s,c("Chr","Pos","Ref","Alt")])
    colnames(mutations) = c("chr","pos","ref","mut")
    mutations$pos=as.numeric(mutations$pos)
    mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% chromosomes,]
    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                       as.numeric(mutations$pos)+1))))
    ntcomp = c(T="A",G="C",C="G",A="T")
    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
    mutations$trinuc_ref_py = mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A","G")) { # Purine base
        mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
        mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    trinuc_mut_mat[n,]=freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat)=full_vec
  rownames(trinuc_mut_mat)=samples
  return(trinuc_mut_mat)
}

plot_comp=function(decomp,signatures.ref=signatures.nature2013,sig_cols=NULL) {
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  require(stringr)
  require(deconstructSigs)
  weights<-as.numeric(decomp$weights)
  signif_sigs=names(decomp$weights)[weights>0]
  weights<-weights[weights>0]
  names(weights)<-signif_sigs
  
  
  #Get color palette
  if(is.null(sig_cols)){
    if(length(weights)<=9){
      sig_cols<-RColorBrewer::brewer.pal(length(weights),name = "Set1")
    } else {
      sig_cols<-colorRampPalette(colors=RColorBrewer::brewer.pal(9,name = "Set1"))(length(weights))
    }
  }
  
  
  weights_matrix=matrix(rep(weights,times=96),nrow=length(weights),dimnames = list(names(weights),colnames(signatures.ref)))
  sigs_matrix=as.matrix(signatures.ref[names(weights),])
  product_df<-as.data.frame(weights_matrix*sigs_matrix)
  p1<-product_df%>%
    tibble::rownames_to_column(var="Signature")%>%
    gather(-Signature,key = "Context",value="Contribution")%>%
    mutate(Change=str_extract(Context,pattern ="\\[.*\\]" ))%>%
    mutate(Change=str_remove_all(Change,"\\[|\\]"))%>%
    mutate(Context=str_replace(Context,pattern = "\\[.*\\]",replacement = "\\."))%>%
    ggplot(aes(x=Context,y=Contribution,fill=Signature))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=sig_cols)+
    facet_grid(cols=vars(Change))+
    guides(fill=guide_legend(position = "inside",nrow = 1))+
    theme_classic()+
    theme(text = element_text(family="Helvetica"),
          axis.text = element_text(size = 7),
          axis.title = element_text(size=8),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3),
          legend.text = element_text(size=6),
          legend.title = element_blank(),
          strip.text = element_text(size=8),
          strip.background = element_rect(fill=NA,linewidth = 0.4),
          legend.spacing = unit(1,"mm"),
          legend.key.size= unit(5,"mm"),
          axis.text.x=element_text(angle=90),
          legend.position.inside = c(0.8,0.8))
  plot(p1)
  Sys.sleep(time = 1)
  return(product_df)
}


#========================================#
# Deconvolute blood signature into SBS1,5 and 19 ####
#========================================#

#Import the mutational signatures extracted using HDP on the 'EM' dataset (normal blood and chemotherapy-exposed blood)
EM_sigs<-read.delim(paste0(root_dir,"Data/mutsig_extraction/EM_extracted_signatures.tsv"),stringsAsFactors = F)
COSMIC_sigs<-as.data.frame(t(MutationalPatterns::get_known_signatures()))
colnames(EM_sigs)<-colnames(COSMIC_sigs)<-colnames(signatures.cosmic)

MutationalPatterns::plot_96_profile(t(signatures.nature2013[c("Signature.1A","Signature.5","Signature.19"),]))
MutationalPatterns::plot_96_profile(t(signatures.cosmic[c("Signature.1","Signature.5","Signature.19"),]))
MutationalPatterns::plot_96_profile(t(COSMIC_sigs[c("SBS1","SBS5","SBS19"),]),ymax=0.4)

MutationalPatterns::plot_96_profile(t(EM_sigs[2,]),ymax=0.07)

signatures_to_use=c("SBS1","SBS5","SBS19")
decomp.blood<-whichSignatures(tumor.ref = EM_sigs,
                              sample.id = "1",
                              signatures.ref = COSMIC_sigs[signatures_to_use,],
                              signature.cutoff = 0.06)
decomp.blood$weights

#How close does the decomposed signature match the real signature?
MutationalPatterns::cos_sim(as.numeric(decomp.blood$tumor),as.numeric(decomp.blood$product))

## Generate Extended Data Fig. 10a ----
myCols=c("#f6568b", "#7dac22", "#0cc0aa")
names(myCols)=c("SBS19","SBS5","SBS1")
product_df<-plot_comp(decomp.blood,signatures.ref=COSMIC_sigs,sig_cols=myCols)
ggsave(filename = paste0(plots_dir,"ExtDatFig10a.pdf"),width=7,height = 2)
sum(abs(decomp$diff))

SBS19_prop=apply(product_df,2,function(x) x[3]/sum(x))


ggsave(filename = paste0(plots_dir,"ExtDatFig10a.pdf"),width=3.5,height=2)

#========================================#
# Now look at contributions in other tissues ####
#========================================#

#Files relating to normal tissue mutation profiling are in this directory
normal_tissue_sigs_dir=paste0(root_dir,"Data/reference_files/normal_tissue_mutation_profiles/")

## 1. PROSTATE----
#This is data from "Development, maturation, and maintenance of human prostate inferred from somatic mutations."
# Cell Stem Cell. 2021 Jul 1;28(7):1262-1274.e5. doi: 10.1016/j.stem.2021.02.005. Epub 2021 Mar 2. PMID: 33657416; PMCID: PMC8260206.

tissue="Prostate"
tissue_dir=paste0(normal_tissue_sigs_dir,tissue)

files=list.files(tissue_dir,pattern=".vcf")
individuals=unique(str_sub(files,1,8))

prostate_muts<-lapply(individuals,function(individual) {
  individual_files=grep(individual,files,value=T)
  individual_muts<-lapply(individual_files,function(file) {
    muts=data.table::fread(input=paste0(tissue_dir,"/",file),skip="#CHROM")
    muts<-muts[,c("#CHROM","POS","REF","ALT")]
    colnames(muts)=c("Chr","Pos","Ref","Alt")
    muts<-muts[!duplicated(muts),]
    muts$Chr=as.character(muts$Chr)
    return(muts)
  })%>%dplyr::bind_rows()
  return(individual_muts)
})%>%dplyr::bind_rows()

if(nrow(prostate_muts)>10000) {
  prostate_muts<-prostate_muts[sample(1:nrow(prostate_muts),size=10000),]
}

prostate_contexts<-mutlist_to_96_contexts(mutlist=prostate_muts,genomeFile = genome_file)

prostate_contexts<-prostate_contexts/sum(prostate_contexts)
prostate_contexts<-as.data.frame(prostate_contexts)
colnames(prostate_contexts)=colnames(signatures.cosmic)

plot_96_profile(t(prostate_contexts),ymax=0.1)

signatures_to_use=c("SBS1","SBS5","SBS19")
decomp.prostate<-whichSignatures(tumor.ref = prostate_contexts,
                                 sample.id = "this_sample",
                                 signatures.ref = COSMIC_sigs[signatures_to_use,],
                                 signature.cutoff = 0.01)


## 2. ENDOMETRIUM----
# This is data from: "The mutational landscape of normal human endometrial epithelium."
# Moore, L., Leongamornlert, D., Coorens, T.H.H. et al.  Nature 580, 640–646 (2020). https://doi.org/10.1038/s41586-020-2214-z

tissue="Endometrium"
tissue_dir=paste0(normal_tissue_sigs_dir,tissue)
files=list.files(tissue_dir,pattern=".txt")

endometrium_muts<-lapply(files,function(file) {
  muts<-read.delim(paste0(tissue_dir,"/",file),stringsAsFactors = F)
  muts<-muts[,c("Chr","Pos","Ref","Alt")]
  colnames(muts)=c("Chr","Pos","Ref","Alt")
  muts$Chr=as.character(muts$Chr)
  return(muts)
})%>%dplyr::bind_rows()


if(nrow(endometrium_muts)>10000) {
  endometrium_muts<-endometrium_muts[sample(1:nrow(endometrium_muts),size=10000),]
}

endometrium_contexts<-mutlist_to_96_contexts(mutlist=endometrium_muts,genomeFile = genome_file)
endometrium_contexts<-endometrium_contexts/sum(endometrium_contexts)
endometrium_contexts<-as.data.frame(endometrium_contexts)
colnames(endometrium_contexts)=colnames(signatures.cosmic)

plot_96_profile(t(endometrium_contexts),ymax=0.1)
decomp.endometrium<-whichSignatures(tumor.ref = endometrium_contexts,
                                    sample.id = "this_sample",
                                    signatures.ref = COSMIC_sigs[signatures_to_use,],
                                    signature.cutoff = 0.01)

## 3. COLON----
# This is data from "Somatic mutation landscapes at single-molecule resolution."
# Abascal, F., Harvey, L.M.R., Mitchell, E. et al. Nature 593, 405–410 (2021). https://doi.org/10.1038/s41586-021-03477-4

seqs<-readxl::read_excel(path = paste0(normal_tissue_sigs_dir,"41586_2021_3477_MOESM3_ESM.xlsx"),sheet=7)%>%
  filter(grepl("coloniccrypts",SampleID))%>%
  tibble::column_to_rownames(var="SampleID")%>%
  as.matrix()%>%
  colSums()
seqs<-seqs/sum(seqs)
seqs<-seqs%>%as.data.frame()%>%
  tibble::rownames_to_column(var="Sub")%>%
  pivot_wider(names_from = "Sub",values_from = ".")
colnames(seqs)<-colnames(signatures.cosmic)
rownames(seqs)<-"ColonicCrypts"
class(seqs)<-"data.frame"

decomp.colon<-whichSignatures(tumor.ref = seqs,
                              sample.id = "ColonicCrypts",
                              signatures.ref = COSMIC_sigs[signatures_to_use,],
                              signature.cutoff = 0.01)


## 4. SMOOTH MUSCLE----
# This is data from "Somatic mutation landscapes at single-molecule resolution."
# Abascal, F., Harvey, L.M.R., Mitchell, E. et al. Nature 593, 405–410 (2021). https://doi.org/10.1038/s41586-021-03477-4

seqs<-readxl::read_excel(path = paste0(normal_tissue_sigs_dir,"41586_2021_3477_MOESM3_ESM.xlsx"),sheet=7)%>%
  filter(grepl("smoothmuscle",SampleID))%>%
  tibble::column_to_rownames(var="SampleID")%>%
  as.matrix()%>%
  colSums()
seqs<-seqs/sum(seqs)
seqs<-seqs%>%as.data.frame()%>%
  tibble::rownames_to_column(var="Sub")%>%
  pivot_wider(names_from = "Sub",values_from = ".")
colnames(seqs)<-colnames(signatures.cosmic)
rownames(seqs)<-"smoothmuscle"
class(seqs)<-"data.frame"

decomp.smoothmuscle<-whichSignatures(tumor.ref = seqs,
                                     sample.id = "smoothmuscle",
                                     signatures.ref = COSMIC_sigs[signatures_to_use,],
                                     signature.cutoff = 0.01)

## 5. CORD BLOOD----
# This is data from "Somatic mutation landscapes at single-molecule resolution."
# Abascal, F., Harvey, L.M.R., Mitchell, E. et al. Nature 593, 405–410 (2021). https://doi.org/10.1038/s41586-021-03477-4

seqs<-readxl::read_excel(path = paste0(normal_tissue_sigs_dir,"41586_2021_3477_MOESM3_ESM.xlsx"),sheet=7)%>%
  filter(!grepl("botseq",SampleID) & !grepl("mung",SampleID) &grepl("cord",SampleID))%>%
  tibble::column_to_rownames(var="SampleID")%>%
  as.matrix()%>%
  colSums()
seqs<-seqs/sum(seqs)
seqs<-seqs%>%as.data.frame()%>%
  tibble::rownames_to_column(var="Sub")%>%
  pivot_wider(names_from = "Sub",values_from = ".")
colnames(seqs)<-colnames(signatures.cosmic)
rownames(seqs)<-"Cord_Blood"
class(seqs)<-"data.frame"

decomp.cordblood<-whichSignatures(tumor.ref = seqs,
                                  sample.id = "Cord_Blood",
                                  signatures.ref = COSMIC_sigs[signatures_to_use,],
                                  signature.cutoff = 0.01)

decomp.cordblood$weights

###Plot these differences ----
all.decomps=list(Adult_blood=decomp.blood,
                 Cord_blood=decomp.cordblood,
                 Prostate=decomp.prostate,
                 Endometrium=decomp.endometrium,
                 Colonic_crypts=decomp.colon,
                 Smooth_muscle=decomp.smoothmuscle)

decomp_weights_list=lapply(all.decomps,function(decomp) decomp$weights)
contributions_df<-Map(vec=decomp_weights_list,type=names(decomp_weights_list),function(vec,type) {
  return(cbind(data.frame(Tissue=type),vec))
})%>%dplyr::bind_rows()%>%
  gather(-Tissue,key="Signature",value="Contribution")

## Generate Extended Data Fig. 10b ----
tissue_order<-rev(c("Cord_blood","Adult_blood","Endometrium","Colonic_crypts","Smooth_muscle","Prostate"))
contributions_plot<-contributions_df%>%
  mutate(Signature=factor(Signature,levels=paste0("SBS",c(1,5,19))),Tissue=factor(Tissue,levels=tissue_order))%>%
  ggplot(aes(y=Tissue,x=Contribution,fill=forcats::fct_rev(Signature)))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=myCols)+
  scale_x_continuous(breaks=seq(0,1,0.2))+
  theme_classic()+
  labs(fill="Signature",y="Normal tissue")+
  my_theme+
  theme(legend.title = element_blank())

ggsave(filename = paste0(plots_dir,"ExtDatFig10b.pdf"),width=3.5,height=2)

#========================================#
# ASSESS CONTRIBUTIONS OF SBS19 TO BLOOD DRIVER LANDSCAPE ####
#========================================#
#Do this using mutation sets from in well-established CH drivers in COSMIC

GENES=c("GNAS","GNB2","KDM6A","CBL","PHIP","BRCC3","PPM1D","TP53","ASXL1","DNMT3A","TET2","SRSF2","U2AF1","EZH2","SF3B1")

out_df<-lapply(GENES,function(GENE) {
  GENE_muts<-read.delim(paste0(root_dir,"Data/reference_files/COSMIC_CH_drivers/COSMIC_",GENE,".tsv"))
  ntcomp = c(T="A",G="C",C="G",A="T")
  GENE_ns<-GENE_muts%>%filter(Primary.Tissue=="Haematopoietic and lymphoid"&
                                !grepl("delins|ins",CDS.Mutation))%>% #only include SNVs
    mutate(Pos=str_extract(CDS.Mutation,"[0-9]+"),NC=str_sub(CDS.Mutation,-3,-1))%>%
    mutate(chr=as.character(stringr::str_split(Genomic.Co.ordinates,pattern = ":",simplify = T)[,1]),pos=stringr::str_split(Genomic.Co.ordinates,pattern = "\\.\\.",simplify = T)[,2])%>%
    mutate(ref=str_sub(NC,1,1),mut=str_sub(NC,3,3))%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select(-Census.Tier.1,-Transcript,-Sample.Name,-Tissue.Subtype.1,-Tissue.Subtype.2,-Histology.Subtype.2,-CGP.Study,-Somatic.Status,-Zygosity,-Pos)
  #View(GENE_ns)
  
  mutations=as.data.frame(GENE_ns)
  chromosomes=c(1:22,"X","Y")
  mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")),]
  mutations$trinuc_ref = as.vector(scanFa(file = genome_file, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                     as.numeric(mutations$pos)+1))))
  
  if(str_sub(mutations$trinuc_ref[1],2,2)!=mutations$ref[1]) {
    mutations$ref<-ntcomp[mutations$ref]
    mutations$mut<-ntcomp[mutations$mut]
  }
  
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  mutations$sub_and_context=paste0(str_sub(mutations$trinuc_ref_py,1,1),"[",mutations$sub,"]",str_sub(mutations$trinuc_ref_py,3,3))
  
  #View(mutations)
  
  mutations%>%
    pull(sub_and_context)%>%table()%>%
    as.data.frame()%>%
    tibble::column_to_rownames(var=".")%>%
    MutationalPatterns::plot_96_profile(ymax=0.42)
  
  df<-mutations%>%
    #filter(!grepl("R882",AA.Mutation))%>%
    mutate(SBS19_prob=SBS19_prop[sub_and_context])%>%
    summarise(sum_SBS19_prob=sum(SBS19_prob),n=n())%>%
    mutate(SBS19_proportion=sum_SBS19_prob/n,Gene=GENE)
  return(df)
})%>%dplyr::bind_rows()

## Generate Extended Data Fig. 10d ----
gene_order<-out_df%>%arrange(desc(n))%>%pull(Gene)
muts_from_SBS19<-out_df%>%
  mutate(SBS19=SBS19_proportion*n)%>%
  mutate(Other=n-SBS19)%>%
  dplyr::select(Gene,SBS19,Other)%>%
  gather(-Gene,key="Signature",value="N_mutations_in_COSMIC")%>%
  mutate(Gene=factor(Gene,levels=gene_order))%>%
  ggplot(aes(y=Gene,x=N_mutations_in_COSMIC,fill=Signature))+
  geom_bar(stat="identity")+
  theme_classic()+
  my_theme+
  labs(x="Number of mutations in COSMIC")

ggsave(filename=paste0(plots_dir,"ExtDatFig10d.pdf"),muts_from_SBS19,width=3.5,height=2)

## Look at total proportion of mutations that are most likely from SBS19
out_df%>%
  mutate(SBS19=SBS19_proportion*n)%>%
  mutate(Other=n-SBS19)%>%
  dplyr::select(Gene,SBS19,Other)%>%
  summarise(n_SBS19=sum(SBS19),n_other=sum(Other))%>%
  mutate(prop_SBS19=n_SBS19/(n_SBS19+n_other))


