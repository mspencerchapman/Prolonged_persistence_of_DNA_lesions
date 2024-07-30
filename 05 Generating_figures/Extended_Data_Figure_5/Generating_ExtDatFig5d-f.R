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

condensed_sig_theme= theme(rect=element_rect(linewidth=0.3),
                           axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(), 
                           axis.title.x = element_text(size = 7),
                           axis.text.x = element_blank(),
                           strip.text.x = element_text(size = 8,face="bold"),
                           strip.text.y = element_blank(),
                           panel.border = element_blank(),
                           axis.line.y = element_blank(),
                           strip.background = element_blank(),
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
# Get the trinucleotide frequencies across the genome ####
#========================================#

#Load the trinucleotide frequencies of the human genome - needed for correction of signatures
Hg_19_strings<-readDNAStringSet(genome_file, format="fasta",
                                nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
triplet_freqs=Biostrings::trinucleotideFrequency(Hg_19_strings)
triplet_freqs_all=colSums(triplet_freqs)
triplet_freqs_all=triplet_freqs_all/sum(triplet_freqs_all) #Normalize the frequencies

#========================================#
# GET A 'normal blood' MUTATIONAL PROFILE ####
#========================================#
#Use Pair11 for this - arbitrary, but has sufficient (but not excessive) numbers of mutations
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
# EXPECTED MULTIPLE INDEPENDENT MUTATION SIGNATURE (PVVs) ####
#========================================#
#This would be expected to reflect:
#1. the likelihood of any mutation occuring TWICE
#2. the overall trinucleotide frequency in the human genome

#Plot the underlying likelihood for mutation of any individual base
#i.e. correct the observed mutational signature in blood for the trinucleotide frequencies in the human genome

EM_sigs<-t(read.delim(paste0(root_dir,"Data/mutsig_extraction/EM_extracted_signatures.tsv")))
rownames(EM_sigs)<-rownames(mut_mat);colnames(EM_sigs)<-c("0","BM_signature","Chemo.sig.1","Chemo.sig.2")
plot_96_profile(EM_sigs,ymax=0.07)

blood_sig<-as.data.frame(EM_sigs[,"BM_signature",drop=F],stringsAsFactors=F)
blood_sig=blood_sig/sum(blood_sig)
blood_sig$trinuc=paste0(substr(rownames(blood_sig),1,1),substr(rownames(blood_sig),3,3),substr(rownames(blood_sig),7,7))
blood_sig$trinuc_lik=sapply(1:nrow(blood_sig),function(i) {
  lik=blood_sig$BM_signature[i]/triplet_freqs_all[blood_sig$trinuc[i]]
  return(lik)
})
blood_sig$trinuc_lik=blood_sig$trinuc_lik/sum(blood_sig$trinuc_lik)
plot_96_profile(blood_sig[,c("trinuc_lik"),drop=F],ymax = 0.3)

## Generate Extended Data Fig. d (i) ------
#Now do the square of the likelihood of a particular mutation (at a particular context) x the trinucleotide frequencies
multi_ind_sig<-blood_sig[,c("trinuc_lik"),drop=F]*blood_sig[,c("trinuc_lik"),drop=F]*as.matrix(triplet_freqs_all[blood_sig$trinuc])
colnames(multi_ind_sig)<-"Independent\nmutations sig"
multi_ind_sig<-multi_ind_sig/sum(multi_ind_sig)
p.MultiSig.1<-plot_96_profile(multi_ind_sig,ymax = 0.5)+condensed_sig_theme
ggsave(p.MultiSig.1,filename = paste0(plots_dir,"ExtDatFig5di"),width=3.5,height=1)

#========================================#
# EXPECTED SPONETANOUS REVERSION MUTATION SIGNATURE (PVVs) ####
#========================================#
#This would be expected to reflect:
#1. the likelihood of an original mutation
#2. the likelihood of the reversion mutation
#3. the overall trinucleotide frequency in the human genome

#The product of the probability of a mutation at a given context x probability of the reversion mutation

#Work out what is the mutation & trinucleotide context of the reversion mutation
spont_rev=data.frame(mut1=rep(NA,96))
spont_rev$mut1=rownames(mut_mat_props)
spont_rev$mut_rev_Ref=substr(spont_rev$mut1,5,5)
spont_rev$mut_rev_Alt=substr(spont_rev$mut1,3,3)
ntcomp = c(T="A",G="C",C="G",A="T")
spont_rev$Sub = paste(spont_rev$mut_rev_Ref,spont_rev$mut_rev_Alt,sep=">")
spont_rev$trinuc_ref=paste0(substr(spont_rev$mut1,1,1),substr(spont_rev$mut1,5,5),substr(spont_rev$mut1,7,7))
spont_rev$trinuc_ref_py = spont_rev$trinuc_ref
for (j in 1:nrow(spont_rev)) {
  if (spont_rev$mut_rev_Ref[j] %in% c("A","G")) { # Purine base
    spont_rev$Sub[j] = paste(ntcomp[spont_rev$mut_rev_Ref[j]],ntcomp[spont_rev$mut_rev_Alt[j]],sep=">")
    spont_rev$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(spont_rev$trinuc_ref[j],split="")[[1]])],collapse="")
  }
}
spont_rev$mut_rev<-paste0(substr(spont_rev$trinuc_ref_py,1,1),"[",spont_rev$Sub,"]",substr(spont_rev$trinuc_ref_py,3,3))

## Generate Extended Data Fig. d (ii) ------
#Now multiply (1) the blood signature probabilities by (2) the blood signature probabilities of the reversion mutation by (3) the trinucleotide frequencies
blood_sig<-cbind(blood_sig[,"BM_signature",drop=F],blood_sig[,"trinuc_lik"]*blood_sig[spont_rev$mut_rev,"trinuc_lik"]*as.matrix(triplet_freqs_all[blood_sig$trinuc]))
colnames(blood_sig)=c("BM_signature","Spont_rev_signature")
blood_sig<-apply(blood_sig,2,function(x) x/sum(x)) #Normalize the signatures

p.SRsig<-plot_96_profile(blood_sig[,"Spont_rev_signature",drop=F],ymax=0.4)+condensed_sig_theme
ggsave(p.SRsig,filename = paste0(plots_dir,"ExtDatFig5dii.pdf"),width=5,height=1.5)

#========================================#
# EXPECTED LOH/ INCORRECT TREE-BUILDING SIGNATURE ####
#========================================#
#This will just reflect the complete BM signature
## Generate Extended Data Fig. d (iii) ------
p.LOH.sig<-plot_96_profile(blood_sig[,"BM_signature",drop=F],ymax=0.07)+condensed_sig_theme
ggsave(p.LOH.sig,filename = paste0(plots_dir,"ExtDatFig5diii.pdf"),width=5,height=1.5)

#========================================#
# ACTUAL SIGNATURE OF 'FAIL' PVVs ####
#========================================#

#Compare the true PVV signature with the spontanous reversion signature
vcfs=read_vcfs_as_granges(vcf_files = paste0(root_dir,"Data/VCFs/Blood_fail_PVVs.vcf"),sample_names = "Blood_fail_PVVs",genome=ref_genome,type="snv")
Blood_PVV_mut_mat=mut_matrix(vcfs,ref_genome = ref_genome)
Blood_PVV_mut_mat_props=as.vector(Blood_PVV_mut_mat/sum(Blood_PVV_mut_mat))

p.failPVV.sig<-plot_96_profile(Blood_PVV_mut_mat,ymax=0.25)+condensed_sig_theme+theme(axis.title.x = element_blank())
ggsave(p.failPVV.sig,filename = paste0(plots_dir,"ExtDatFig5e.pdf"),width=5,height=1.5)

#Compare the fail signature with different possible mechanisms
cos_sim(x=as.vector(Blood_PVV_mut_mat/sum(Blood_PVV_mut_mat)),y=mut_mat_props[,"Spont_rev_signature"])
cos_sim(x=as.vector(Blood_PVV_mut_mat/sum(Blood_PVV_mut_mat)),y=mut_mat_props[,"Blood_signature"])
cos_sim(x=as.vector(Blood_PVV_mut_mat/sum(Blood_PVV_mut_mat)),y=as.vector(multi_ind_sig[,1]))


#========================================#
# CO-OCCURRENCE OF MUTATIONS BETWEEN INDIVIDUALS ####
#========================================#
#Compare multiple occurrences of mutations *within* versus *between* individuals i.e. mutation hotspots vs lesion persistence
#Need to have lots of RAM (approximately 64GB) to import so much mutation info at the same time
data_sets=c("MSC_BMT","MF","EM")
lesion_seg_input_dir=paste0(root_dir,"/Data/input_data")
all_muts_df=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  all_data_set_mat<-lapply(data_set_samples,function(sample_ID) {
    cat(sample_ID,sep="\n")
    load(get_file_paths_and_project(Sample_ID =sample_ID,dataset=data_set,input_data_dir = lesion_seg_input_dir)$filtered_muts_path)
    mutation_mat<-filtered_muts$COMB_mats.tree.build$mat%>%dplyr::mutate(Sample_ID=sample_ID,.before=1)%>%dplyr::select(Sample_ID,mut_ref,Chrom,Pos,Ref,Alt)
    rm(filtered_muts)
    return(mutation_mat)
  })%>%dplyr::bind_rows()
  return(all_data_set_mat)
})%>%dplyr::bind_rows()

all_trees=lapply(data_sets,function(data_set) {
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  all_data_set_trees<-lapply(data_set_samples,function(sample_ID) {
    tree<-read.tree(get_file_paths_and_project(Sample_ID =sample_ID,dataset=data_set,input_data_dir = lesion_seg_input_dir)$tree_file_path)
    return(tree)
  })
  names(all_data_set_trees)<-data_set_samples
  return(all_data_set_trees)
})%>%unlist(recursive = F)

#Summarise the number of mutations per individual
summary_df=data.frame(Sample_ID=unique(all_muts_df$Sample_ID))
summary_df$n_muts<-sapply(summary_df$Sample_ID,function(sample_ID) {nrow(all_muts_df%>%filter(Sample_ID==sample_ID))})

#Find all mutations that are present more than once in the full dataset
dup_muts<-all_muts_df$mut_ref[duplicated(all_muts_df$mut_ref)]

#Filter the dataframe to include only these mutation
dup_muts_df<-all_muts_df%>%
  filter(mut_ref%in%dup_muts)%>%
  arrange(mut_ref)

#Annotate if mutations are on private or shared branches
dup_muts_df$branch_type<-sapply(1:nrow(dup_muts_df),function(i) {ifelse(dup_muts_df$node[i] <= length(all_trees[[dup_muts_df$Sample_ID[i]]]$tip.label),"private","shared")})

#Give a sense of how many doublets/ triplets etc there are
table(table(dup_muts_df$mut_ref))

#Create a matrix of these mutations and which samples they are present in
dup_mat<-matrix(0,nrow = length(dup_muts),ncol=length(unique(all_muts_df$Sample_ID)))
rownames(dup_mat)<-dup_muts
colnames(dup_mat)<-unique(all_muts_df$Sample_ID)

for(sample_ID in unique(all_muts_df$Sample_ID)) {
  Sample_dups<-dup_muts_df%>%filter(Sample_ID==sample_ID)%>%pull(mut_ref)
  dup_mat[Sample_dups,sample_ID]<-1
}
head(dup_mat[order(rowSums(dup_mat),decreasing = F),])

#Now do a similar matrix, but one which shows the number of overlaps, relative to the number of opportunities for overlaps (n_mut1 * nmut2)
sample_by_sample_dups<-matrix(NA,nrow=length(unique(all_muts_df$Sample_ID)),ncol=length(unique(all_muts_df$Sample_ID)),dimnames = list(unique(all_muts_df$Sample_ID),(unique(all_muts_df$Sample_ID))))

for(i in 1:nrow(sample_by_sample_dups)) {
  for(j in 1:ncol(sample_by_sample_dups)) {
    Sample1=rownames(sample_by_sample_dups)[i];Sample2=colnames(sample_by_sample_dups)[j]
    nmut_1=summary_df$n_muts[summary_df$Sample_ID==Sample1]; nmut_2=summary_df$n_muts[summary_df$Sample_ID==Sample2]
    sample_by_sample_dups[i,j] <- sum(dup_mat[,i]==1 & dup_mat[,j]==1)/prod(nmut_1,nmut_2)
  }
  #Set the same sample comparisons to 'NA'
  sample_by_sample_dups[i,i]<-NA
}

#Show this as a heatmap
library(pheatmap)
pheatmap(sample_by_sample_dups,cluster_cols=F,cluster_rows = F)

###Look at the mutation signatures of these duplicates
mutations<-tidyr::separate(data.frame(mut_ref=rownames(dup_mat)),col = mut_ref,into = c("chr","pos","ref","mut"),sep = "-",remove = F)%>%mutate(SampleID="BS_duplicates",.before=1)
mutations$pos=as.numeric(mutations$pos)
library("GenomicRanges")
library("Rsamtools")
library("MASS")
mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
mutations$trinuc_ref = as.vector(scanFa(genome_file, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
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

plot_96profile(mutations,colnames=c("SampleID","chr","pos","ref","mut"),genomeFile = genome_file,savefile=paste0(plots_dir,"ExtDatFig5f.pdf"))

