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
# Define custom function ####
#========================================#

#Updated version of this function that allows for a 'select vector'
create_vcf_files=function(mat, select_vector = NULL,from_mut_ref=F,mut_ref_col="mut_ref") {
  if(from_mut_ref){
    if(is.null(select_vector)) {vcf_file = stringr::str_split(mat[,mut_ref_col],pattern = "-",simplify=T)} else {vcf_file = stringr::str_split(mat[select_vector,mut_ref_col],pattern = "-",simplify=T)}
    vcf_file=as.data.frame(vcf_file,stringsAsFactors=F)
  } else {
    if(is.null(select_vector)) {vcf_file = mat[,c("Chrom","Pos","Ref","Alt")]} else {vcf_file = mat[select_vector,c("Chrom","Pos","Ref","Alt")]}
  }
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  vcf_file = vcf_file[,c(1,2,8,3,4,7,6,5)]
  return(vcf_file)
}

#Updated version of this function that allows for a 'select vector' - used by the 'write_branch_muts_as_vcfs' function
write.vcf=function(details,vcf_path,select_vector=NULL,from_mut_ref=F,mut_ref_col="mut_ref",vcf_header_path="~/Documents/vcfHeader.txt") {
  vcf=create_vcf_files(mat=details,select_vector=select_vector,from_mut_ref = from_mut_ref,mut_ref_col=mut_ref_col)
  write.table(vcf,sep = "\t", quote = FALSE,file=paste0(vcf_path,".temp"),row.names = F)
  system(paste0("cat ",vcf_header_path," ",vcf_path,".temp > ",vcf_path))
  system(paste0("rm ",vcf_path,".temp"))
}

write_branch_muts_as_vcfs=function(details,Sample_ID,output_dir,vcf_header_path="/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt") {
  out=sapply(unique(details$node),function(node) {
    write.vcf(details,
              vcf_path=paste0(output_dir,"/",Sample_ID,"_",node,".vcf"),
              select_vector = details$node==node,
              from_mut_ref=T,
              vcf_header_path = vcf_header_path)
  })
  return(NULL)
}

#========================================#
# Perform the lesion segregation analysis if not already done ####
#========================================#

# Define autosomal chromosomes & additional file paths
chromosomes <- seqnames(get(ref_genome))[1:22]
output_dir=paste0(root_dir,"Data/lesion_segregation_analysis")
system(paste0("mkdir -p ",output_dir))

data_sets=c("NW","MF","EM","KY","SN","MSC_BMT","MSC_fetal")

#Run the analysis
#If actually running the lesion segregation analysis from scratch, this takes a long time
summary_table_list=lapply(data_sets,function(data_set) {
  print(data_set)
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  lapply(data_set_samples,function(Sample_ID){
    print(Sample_ID)
    complete=file.exists(paste0(output_dir,"/",Sample_ID,"_binom",".tsv"))&
      file.exists(paste0(output_dir,"/",Sample_ID,"_wald",".tsv"))&
      file.exists(paste0(output_dir,"/",Sample_ID,"_rl20",".tsv"))
    if(!complete) {
      sample_info=get_file_paths_and_project(dataset=data_set,Sample_ID=Sample_ID,input_data_dir = lesion_seg_input_dir)
      load(sample_info$filtered_muts_path)
      details<-filtered_muts$COMB_mats.tree.build$mat
      system(paste0("mkdir -p ",output_dir,"/temp"))
      write_branch_muts_as_vcfs(details,Sample_ID,output_dir=paste0(output_dir,"/temp"),vcf_header_path = vcf_header_path)
      vcf_files=grep(Sample_ID,list.files(path=paste0(output_dir,"/temp"),pattern = ".vcf",full.names = T),value = T)
      sample_names = gsub(".vcf","",grep(Sample_ID,list.files(path=paste0(output_dir,"/temp"),pattern = ".vcf"),value = T))
      vcfs=read_vcfs_as_granges(vcf_files = vcf_files,sample_names = sample_names,genome=ref_genome)
      
      if(!file.exists(paste0(output_dir,"/",Sample_ID,"_binom",".tsv"))){
        lesion_segregation <- calculate_lesion_segregation(vcfs, sample_names)
        write.table(lesion_segregation,file = paste0(output_dir,"/",Sample_ID,"_binom",".tsv"),quote = F,row.names = F,sep = "\t")
      }
      
      if(!file.exists(paste0(output_dir,"/",Sample_ID,"_wald",".tsv"))){
        lesion_segregation_wald <- calculate_lesion_segregation(vcfs, sample_names,test = "wald-wolfowitz")
        write.table(lesion_segregation_wald,file=paste0(output_dir,"/",Sample_ID,"_wald",".tsv"),quote = F,row.names = F,sep = "\t")
      }
      
      if(!file.exists(paste0(output_dir,"/",Sample_ID,"_rl20",".tsv"))){
        lesion_segregation_rl20 <- calculate_lesion_segregation(vcfs,sample_names,test = "rl20",ref_genome = ref_genome,chromosomes = chromosomes)
        write.table(lesion_segregation_rl20,file = paste0(output_dir,"/",Sample_ID,"_rl20",".tsv"),quote = F,row.names = F,sep = "\t")
      }

    } else {
      print(paste("Analysis for",Sample_ID,"already completed."))
    }
  })
})

#Re-import analysis
lesion_seg_results_list=lapply(data_sets,function(data_set) {
  print(data_set)
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  data_set_out=lapply(data_set_samples,function(Sample_ID){
    lesion_segregation<-read.delim(file = paste0(output_dir,"/",Sample_ID,"_binom",".tsv"))
    lesion_segregation_wald<-read.delim(file=paste0(output_dir,"/",Sample_ID,"_wald",".tsv"))
    lesion_segregation_rl20<-read.delim(file = paste0(output_dir,"/",Sample_ID,"_rl20",".tsv"))
    return(list(data_set=data_set,lesion_segregation=lesion_segregation,lesion_segregation_wald=lesion_segregation_wald,lesion_segregation_rl20=lesion_segregation_rl20))
  })
  return(data_set_out)
})
lesion_seg_results_list=unlist(lesion_seg_results_list,recursive=F)
names(lesion_seg_results_list)<-sapply(lesion_seg_results_list,function(list) {
  #Re-extract the sample name
  sample_name_vec=head(stringr::str_split(list$lesion_segregation$sample_name[1],pattern="_")[[1]],-1)
  sample=paste(sample_name_vec,collapse = "_")
  return(sample)
})

#Create a summary table of how many branches show lesion segregation by the different tests
out_list=lapply(lesion_seg_results_list,function(list) {
  #Re-extract the sample name
  sample_name_vec=head(stringr::str_split(list$lesion_segregation$sample_name[1],pattern="_")[[1]],-1)
  sample=paste(sample_name_vec,collapse = "_")
  n_binom=sum(list$lesion_segregation$fdr<0.05,na.rm = T)
  n_wald=sum(list$lesion_segregation_wald$fdr<0.05,na.rm = T)
  n_rl20_6_or_over=sum(list$lesion_segregation_rl20$rl20>=6,na.rm = T)
  return(data.frame(data_set=list$data_set,Sample_ID=sample,n_binom=n_binom,n_wald=n_wald,n_rl20=n_rl20_6_or_over))
})
lesion_seg_df=dplyr::bind_rows(out_list)
lesion_seg_df$any_pos<-apply(lesion_seg_df[,3:5],1,function(x) {return(any(as.logical(x)))})
lesion_seg_df$all_pos<-apply(lesion_seg_df[,3:5],1,function(x) {return(all(as.logical(x)))})
head(lesion_seg_df)

## Generate Figure 5b ------
### For PX001 correlation of lesion segregation with chemo sig 1 proportion ----
data_set="EM"
HDP_dir=paste0(root_dir,"Data/HDP/")
exposures_file_path=paste0(HDP_dir,data_set,"/exposures.csv")

if(file.exists(exposures_file)) {
  exposures_df<-read_csv(exposures_file)
} else {
  #If the exposures file hasn't been created, can get it from the HDP output
  mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
  mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
  key_table=read.table(paste0(HDP_folder,"/key_table.txt"))
  
  exposures_df=create_exposures_df(HDP_multi = mut_example_multi,trinuc_mut_mat = mutations,key_table=key_table,sep="-")
  write.csv(exposures_df,exposures_file,row.names = F)
}

tree=read.tree(get_file_paths_and_project(dataset="EM",Sample_ID = "PX001_2_01",input_data_dir = lesion_seg_input_dir)$tree_file_path)

#Combine the lesion segregation results data frames into a single df and add the exposure information
lesion_seg_res_list=lapply(lesion_seg_results_list,function(list) {
  if(class(list)=="list"){
    if(list$data_set==data_set) {
      df=left_join(list$lesion_segregation,list$lesion_segregation_wald,by="sample_name")%>%left_join(list$lesion_segregation_rl20,by="sample_name")
      df$node=as.numeric(str_split(df$sample_name,pattern="_",simplify=T)[,4])
      df$Sample_ID=apply(str_split(df$sample_name,pattern="_",simplify=T)[,1:3],1,paste0,collapse="_")
      df<-left_join(df,exposures_df,by=c("Sample_ID","node"))
      return(df) 
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
})

#Now review the chemo.sig.1 proportions samples with lesion segregation compared to those without
#Chemo.sig.1 is N2 in this exposures dataframe
p.ls.1<-lesion_seg_res_list$PX001_2_01%>%
  filter(!is.na(rl20))%>%
  mutate(Result=factor(ifelse(rl20>5&fdr.y<0.05,"Lesion\nsegregation","No lesion\nsegregation")))%>%
  ggplot(aes(x=Result,y=N2,fill=Result))+
  geom_violin(col=NA)+
  geom_jitter(col="gray",width=0.2,alpha=0.3)+
  #scale_x_continuous()+
  theme_classic()+
  my_theme+
  labs(x="",y="Proportion of chemo signature")+
  theme(axis.title.x = element_blank(),legend.position = "none")

ggsave(p.ls.1,filename=paste0(plots_dir,"Fig5b.pdf"),width=2,height=2)


## Generate Figure 5c ------
#For lung samples correlation of lesion segregation with APOBEC proportion
data_set="KY"
exposures_file_path=paste0(HDP_dir,data_set,"/exposures.csv")

if(file.exists(exposures_file)) {
  exposures_df<-read_csv(exposures_file)
} else {
  #If the exposures file hasn't been created, can get it from the HDP output
  mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
  mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
  key_table=read.table(paste0(HDP_folder,"/key_table.txt"))
  
  exposures_df=create_exposures_df(HDP_multi = mut_example_multi,trinuc_mut_mat = mutations,key_table=key_table,sep="-")
  write.csv(exposures_df,exposures_file,row.names = F)
}

#APOBEC is N3
KY_lesion_seg_results_list=lapply(lesion_seg_results_list,function(list) {
  if(class(list)=="list"){
    if(list$data_set=="KY"){
      df=left_join(list$lesion_segregation,list$lesion_segregation_wald,by="sample_name")%>%left_join(list$lesion_segregation_rl20,by="sample_name")
      df$node=as.numeric(str_split(df$sample_name,pattern="_",simplify=T)[,2])
      df$Sample_ID=str_split(df$sample_name,pattern="_",simplify=T)[,1]
      df<-left_join(df,exposures_df,by=c("Sample_ID","node"))
      return(df) 
    }else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
})

KY_ls_df=dplyr::bind_rows(KY_lesion_seg_results_list)

p.ls.2<-KY_ls_df%>%
  filter(!is.na(rl20))%>%
  mutate(Result=factor(ifelse(rl20>5&fdr.y<0.05,"Lesion\nsegregation","No lesion\nsegregation")))%>%
  ggplot(aes(x=Result,y=N2,fill=Result))+
  geom_violin(col=NA)+
  geom_jitter(col="gray",width=0.2,alpha=0.25)+
  theme_classic()+
  my_theme+
  labs(x="",y="Proportion of APOBEC muts")+
  theme(axis.title.x = element_blank(),legend.position = "none")

ggsave(p.ls.2,filename=paste0(plots_dir,"Fig5c.pdf"),width=2,height=2)


## Generate Figure 5d ------
### Plot the lesion segregation plots for the top branches ----
for(i in 1:nrow(lesion_seg_df)) {
  print(i)
  if(lesion_seg_df$data_set[i]=="SN") {next}
  if(lesion_seg_df$all_pos[i]) {
    print(lesion_seg_df$Sample_ID[i])
    sample_info=get_file_paths_and_project(dataset=lesion_seg_df$data_set[i],Sample_ID=lesion_seg_df$Sample_ID[i],input_data_dir = lesion_seg_input_dir)
    load(sample_info$filtered_muts_path)
    details<-filtered_muts$COMB_mats.tree.build$mat
    tree=read.tree(sample_info$tree_file_path)
    
    #Find nodes that are positive by any of the 3 metrics
    sample_results<-lesion_seg_results_list[[lesion_seg_df$Sample_ID[i]]]
    
    #Add a node column
    sample_results<-lapply(sample_results, function(df) {
      if(class(df)!="data.frame") {stop(return(df))}
      df$node<-sapply(df$sample_name,function(sample_name) {tail(stringr::str_split(sample_name,pattern="_")[[1]],n=1)})
      return(df)
    })
    #binom test nodes
    binom_pos<-sample_results$lesion_segregation%>%filter(fdr<0.05)%>%pull(sample_name)
    if(length(binom_pos)>0) {binom_pos_nodes<-unlist(sapply(str_split(binom_pos,pattern = "_"),function(x) {as.numeric(tail(x,n=1))}))} else {binom_pos_nodes<-NULL}
    
    #wald test nodes
    wald_pos<-sample_results$lesion_segregation_wald%>%filter(fdr<0.05)%>%pull(sample_name)
    if(length(wald_pos)>0) {wald_pos_nodes<-unlist(sapply(str_split(wald_pos,pattern = "_"),function(x) {x<-x[!is.na(x)];as.numeric(tail(x,n=1))}))} else {wald_pos_nodes<-NULL}
    
    #rl20 pos nodes
    rl20_pos<-sample_results$lesion_segregation_rl20%>%filter(rl20>5)%>%pull(sample_name)
    if(length(rl20_pos)>0) {rl20_pos_nodes<-unlist(sapply(str_split(rl20_pos,pattern = "_"),function(x) {x<-x[!is.na(x)];as.numeric(tail(x,n=1))}))} else {rl20_pos_nodes<-NULL}
    
    #Union of all 3 approaches
    #all_pos=unique(c(binom_pos_nodes,wald_pos_nodes,rl20_pos_nodes))
    all_pos=intersect(unique(binom_pos_nodes),intersect(unique(wald_pos_nodes),unique(rl20_pos_nodes)))
    print(all_pos)
    
    #Plot the lesion segregation plot for any branches that are significant by all measures
    if(!is.null(all_pos)&length(all_pos)>0){
      print("Importing the vcfs")
      system(paste0("mkdir -p ",output_dir,"/temp"))
      VCF_files=paste0(output_dir,"/temp/",paste(lesion_seg_df$Sample_ID[i],all_pos,sep="_"),".vcf")
      if(!all(file.exists(VCF_files))) {
        for(node in all_pos) {
          write.vcf(details,
                    vcf_path=paste0(output_dir,"/temp/",lesion_seg_df$Sample_ID[i],"_",node,".vcf"),
                    select_vector = details$node==node,
                    from_mut_ref=T,
                    vcf_header_path = vcf_header_path)
        }
      }
      lesion_seg_vcfs=read_vcfs_as_granges(vcf_files=VCF_files,sample_names = paste(lesion_seg_df$Sample_ID[i],all_pos,sep="_"),genome=ref_genome,predefined_dbs_mbs=F)
      for(j in 1:length(lesion_seg_vcfs)){
        print(paste("Plotting lesion segregation plot",j),sample_name=paste(lesion_seg_df$Sample_ID[i],all_pos[j],sep=","))
        system(paste0("mkdir -p ",paste0(plots_dir,"lesion_segregation_plots"))) #Create a subdirectory within the plots directory for these
        p<-plot_lesion_segregation_MOD(lesion_seg_vcfs[[j]])
        p<-p+ggtitle(paste0(lesion_seg_df$Sample_ID[i],", branch ",all_pos[j]))+theme(panel.grid.minor=element_blank(),title = element_text(size=8))
        ggsave(paste0(plots_dir,"lesion_segregation_plots/",lesion_seg_df$Sample_ID[i],"_",all_pos[j],".pdf"),p,width=7,height=2.5)
      }
    }
  } else {
    next
  }
}




##NOW FOR LIVER SAMPLES
#For liver samples correlate with new signature proportion ("N")
exposures_df_liver=read.delim("/lustre/scratch119/realdata/mdt1/team154/ms56/lesion_segregation/x.hdp.prop.liver.paper.tsv",header = T)
exposures_df_liver$ID<-gsub("C","",exposures_df_liver$ID)
exposures_df_liver<-exposures_df_liver%>%
  tidyr::separate(col = ID,into=c("Sample_ID","node"),sep="_")%>%
  mutate(node=as.numeric(node))

SN_lesion_seg_results_list=lapply(lesion_seg_results_list,function(list) {
  if(class(list)=="list"){
    if(list$data_set=="SN"){
      df=left_join(list$lesion_segregation,list$lesion_segregation_wald,by="sample_name")%>%left_join(list$lesion_segregation_rl20,by="sample_name")
      df$node=as.numeric(str_split(df$sample_name,pattern="_",simplify=T)[,2])
      df$Sample_ID=str_split(df$sample_name,pattern="_",simplify=T)[,1]
      df<-left_join(df,exposures_df_liver,by=c("Sample_ID","node"))
      return(df) 
    }else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
})
SN_ls_df=dplyr::bind_rows(SN_lesion_seg_results_list)

p.ls.3<-SN_ls_df%>%
  filter(!is.na(rl20))%>%
  mutate(Result=factor(ifelse(rl20>5&fdr.y<0.05,"Lesion\nsegregation","No lesion\nsegregation")))%>%
  ggplot(aes(x=Result,y=N5,fill=Result))+
  geom_violin(col=NA)+
  geom_jitter(col="gray",width=0.2)+
  #scale_x_continuous()+
  theme_classic()+
  my_theme+
  labs(x="Lesion segregation result",y="Proportion of novel signature mutations")

SN_ls_df%>%
  filter(!is.na(rl20))%>%
  mutate(Result=factor(ifelse(rl20>5&fdr.y<0.05,"Lesion\nsegregation","No lesion\nsegregation")))%>%
  dplyr::filter(Result=="Lesion\nsegregation")
