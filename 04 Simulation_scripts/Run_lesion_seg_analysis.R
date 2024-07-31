#CALCULATING LESION SEGREGATION ACROSS BRANCHES
library(GenomicRanges)
library(IRanges)
library("Rsamtools")
library("MASS")
library(stringr)
library(dplyr)
library(tidyr)
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome,character.only=TRUE)
options(stringsAsFactors = FALSE)

my_working_directory<-getwd()

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files[-2],source)


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

write.vcf=function(details,vcf_path,select_vector=NULL,from_mut_ref=F,mut_ref_col="mut_ref",vcf_header_path="~/Documents/vcfHeader.txt") {
  vcf=create_vcf_files(mat=details,select_vector=select_vector,from_mut_ref = from_mut_ref,mut_ref_col=mut_ref_col)
  write.table(vcf,sep = "\t", quote = FALSE,file=paste0(vcf_path,".temp"),row.names = F)
  system(paste0("cat ",vcf_header_path," ",vcf_path,".temp > ",vcf_path))
  system(paste0("rm ",vcf_path,".temp"))
}

write_branch_muts_as_vcfs=function(details,Sample_ID,output_dir,vcf_header_path="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt") {
  out=sapply(unique(details$node),function(node) {
    write.vcf(details,
              vcf_path=paste0(output_dir,"/",Sample_ID,"_",node,".vcf"),
              select_vector = details$node==node,
              from_mut_ref=T,
              vcf_header_path = vcf_header_path)
  })
  return(NULL)
}

# Define autosomal chromosomes
chromosomes <- seqnames(get(ref_genome))[1:22]

lesion_seg_input_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data"
output_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/lesion_segregation_analysis"
system(paste0("mkdir -p ",output_dir))

data_sets=c("NW","MF","EM","KY","PR","MSC_BMT","MSC_fetal")
summary_table_list=lapply(data_sets,function(data_set) {
  print(data_set)
  data_set_samples=readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  lapply(data_set_samples,function(Sample_ID){
    sample_info=get_file_paths_and_project(dataset=data_set,Sample_ID=Sample_ID)
    load(sample_info$filtered_muts_path)
    details<-filtered_muts$COMB_mats.tree.build$mat
    system(paste0("mkdir -p ",output_dir,"/temp"))
    write_branch_muts_as_vcfs(details,Sample_ID,output_dir=paste0(output_dir,"/temp"))
    vcf_files=grep(Sample_ID,list.files(path=paste0(output_dir,"/temp"),pattern = ".vcf",full.names = T),value = T)
    sample_names = gsub(".vcf","",grep(Sample_ID,list.files(path=paste0(output_dir,"/temp"),pattern = ".vcf"),value = T))
    vcfs=read_vcfs_as_granges(vcf_files = vcf_files,sample_names = sample_names,genome=ref_genome)
    
    lesion_segregation <- calculate_lesion_segregation(vcfs, sample_names)
    lesion_segregation_wald <- calculate_lesion_segregation(vcfs, sample_names,test = "wald-wolfowitz")
    lesion_segregation_rl20 <- calculate_lesion_segregation(vcfs,sample_names,test = "rl20",ref_genome = ref_genome,chromosomes = chromosomes)
    
    write.table(lesion_segregation,file = paste0(output_dir,"/",Sample_ID,"_binom",".tsv"),quote = F,row.names = F,sep = "\t")
    write.table(lesion_segregation_wald,file=paste0(output_dir,"/",Sample_ID,"_wald",".tsv"),quote = F,row.names = F,sep = "\t")
    write.table(lesion_segregation_rl20,file = paste0(output_dir,"/",Sample_ID,"_rl20",".tsv"),quote = F,row.names = F,sep = "\t")
    
  })
})
