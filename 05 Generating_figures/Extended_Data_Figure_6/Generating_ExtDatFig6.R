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
# CREATE THE SUMMARY PLOTS ####
#========================================#

all_IDs=unique(mutations%>%filter(data_set!="SN")%>%pull(Sample_ID)) #Can't do this for the liver dataset as these are inferred clones

## Do MAV summaries for the two individuals in the figure -----
# However, can do for any of the individuals

## Generate Extended Data Figures 6a,b ----
for(ID in c("KX003_5_01","KX008_2_01")) {
  
  ## Import the data----
  cat(ID,sep="/n")
  data_set=mutations%>%filter(Sample_ID==ID)%>%pull(data_set)%>%.[1]
  info=get_file_paths_and_project(dataset=data_set,Sample_ID = ID,input_data_dir = lesion_seg_input_dir)
  load(info$filtered_muts_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  matrices<-list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR)
  tree<-read.tree(info$tree_file_path)
  tree_no_ancestral<-drop.tip(tree,tip="Ancestral")

  ## Plot the MAVs ----
  cat("Plotting MAVs",sep="/n")
  
  ## Subset the mutation matrix to get the relevant mutations ----
  individual_MAVs<-mutations%>%
    filter(Class!="FAIL" & Sample_ID==ID & Type=="MAV" & !is.na(lesion_node))%>%
    arrange(lesion_timing)
  
  
  if(nrow(individual_MAVs>0)) {
    #Create the colour heatmap for plotting
    MAV_cols=RColorBrewer::brewer.pal(n=11,"RdBu")[c(2,9)]
    MAV_cols_across_tree<-lapply(1:nrow(individual_MAVs),function(i) {
      cat(i)
      mut1=individual_MAVs$mut_ref1[i];NV1=matrices$NV[mut1,tree_no_ancestral$tip.label]
      mut2=individual_MAVs$mut_ref2[i];NV2=matrices$NV[mut2,tree_no_ancestral$tip.label]
      NR=matrices$NR[mut1,tree_no_ancestral$tip.label]+matrices$NV[mut2,tree_no_ancestral$tip.label]
      this_MAV_cols<-mapply(NV1=NV1,NV2=NV2,NR=NR, function(NV1,NV2,NR) {
        col=NA
        if(NV1>0 & NV2==0) {
          peak_col=MAV_cols[1]
        } else if(NV1==0 & NV2>0) {
          peak_col=MAV_cols[2]
        } else if(NV1>0 & NV2>0) {
          peak_col=colorRampPalette(colors=MAV_cols)(100)[as.integer(round(100*NV2/(NV1+NV2)))]
        } else {
          col="#FFFFFF"
        }
        
        if(is.na(col)) {
          col=colorRampPalette(colors = c("#FFFFFF",peak_col))(100)[as.integer(round(100*calculate_vaf(NV=(NV1+NV2),NR=NR)))]
        }
        return(col)
      })
      
      return(this_MAV_cols)
    })%>%dplyr::bind_rows()%>%as.matrix()
    rownames(MAV_cols_across_tree)<-paste(individual_MAVs$mut_ref1,individual_MAVs$Alt2,sep="/")
    
    if("Ancestral"%in%tree$tip.label) {
      MAV_cols_across_tree<-cbind(MAV_cols_across_tree,matrix(rep("#FFFFFF",nrow(MAV_cols_across_tree)),ncol = 1,dimnames = list(rownames(MAV_cols_across_tree),"Ancestral")))
    }
    
    #Now do the plotting
    MAV_label_cols=colorRampPalette(colors=RColorBrewer::brewer.pal(n=8,name="Dark2"))(nrow(individual_MAVs))
    
    pdf(paste0(plots_dir,"MAV_summary_",ID,".pdf"),width=7,height=7)
    tree=plot_tree(tree = tree,cex.label = 0,lwd=ifelse(ID=="KX004_5_01",0.25,0.5),vspace.reserve = 1,cex.axis=0.5,tck=-0.01)
    add_mut_heatmap(tree=tree,heatmap=MAV_cols_across_tree[,tree$tip.label,drop=F],border="gray",heatmap_bar_height=0.05,cex.label = 0.25,label.cols = MAV_label_cols)
    if(data_set!="SN"){
      temp=sapply(1:nrow(individual_MAVs),function(i) {
        
        #Highlight the lesion node
        info=get_edge_info(tree=tree,details=NULL,node=individual_MAVs$lesion_node[i])
        plotrix::draw.circle(x=info$x,y=info$yb,radius=2,border=NA,col=MAV_label_cols[i]) ##May need to adjust the radius here to make it look sensible
        
        #Plot the lesion path for removed MAVs
        if(individual_MAVs$lesion_node[i]!=individual_MAVs$lesion_repair_node[i]) {
          LN_ancestors<-get_ancestral_nodes(node = individual_MAVs$lesion_node[i],edge=tree$edge)
          LRN_ancestors<-get_ancestral_nodes(node = individual_MAVs$lesion_repair_node[i],edge=tree$edge)
          lesion_path_nodes=LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
          temp=add_annotation(tree=tree,
                              annot_function = highlight_nodes,
                              nodes=lesion_path_nodes,
                              col=MAV_label_cols[i],
                              lwd=ifelse(ID=="KX004_5_01",1,2))
        }
      })
    }
    dev.off()
  }
}

