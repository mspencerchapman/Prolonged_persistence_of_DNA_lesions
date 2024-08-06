#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","ggtree","RColorBrewer","tibble","ape","dichromat","dplyr","tidyr","stringr","readr","phytools","devtools")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("NickWilliamsSanger/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Define custom functions ####
#========================================#

create_exposures_df=function(HDP_multi,trinuc_mut_mat,key_table,minimum_branch_muts=50,sep="_") {
  library(hdp)
  library(tidyr)
  library(dplyr)
  sample_remove=rownames(trinuc_mut_mat)[rowSums(trinuc_mut_mat)<minimum_branch_muts]
  trinuc_mut_mat=trinuc_mut_mat[!rownames(trinuc_mut_mat)%in%sample_remove,]
  key_table=key_table[!key_table$Sample%in%sample_remove,]
  freq=nrow(trinuc_mut_mat)
  
  dp_distn <- comp_dp_distn(HDP_multi)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(trinuc_mut_mat),,drop=FALSE])
  colnames(exposures)=rownames(trinuc_mut_mat)
  rownames(exposures)<-paste0("N",rownames(exposures))
  sigs=rownames(exposures)
  sig_profiles=HDP_multi@comp_categ_distn$mean
  
  exposures_df<-as.data.frame(t(exposures),stringsAsFactors=F)%>%
    tibble::rownames_to_column("branch")%>%
    tidyr::separate(col="branch",into=c("node","Sample_ID"),sep=sep)%>%
    mutate(node=as.numeric(node))
  return(exposures_df)
}

plot_sig_tree=function(tree,sample_ID,exposures,save_path=NULL,...) {
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(ggtree)
  
  #Fortify the
  tree_df=fortify(tree)
  cols=c(brewer.pal(12, "Paired"),"magenta","firebrick")
  
  #Reduce the exposure df to only those for this sample
  exposures<-exposures%>%filter(Sample_ID==sample_ID)
  
  #Define the signature names
  sigs<-exposures%>%dplyr::select(-node,-Sample_ID)%>%colnames()
  
  #Assign colors to signatures
  cols<-cols[1:length(sigs)]
  names(cols)=sigs
  
  #Now generate the plot
  if(!is.null(save_path)) {pdf(save_path,height=20,width=15)}
  
  plot(tree,label.offset=0.01*max(tree_df$x),...)
  for (k in 1:nrow(exposures)){
    n=exposures$node[k]
    x_end=tree_df$x[n]
    x_start=tree_df$x[tree_df$parent[n]]
    x_intv=x_end-x_start
    y=node.height(tree)[n]
    tipnum=sum(tree_df$isTip)
    for (s in sigs){
      x_end=x_start+exposures[[s]][k]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s])
      x_start=x_end
    }
  }
  axisPhylo(side = 1,backward=F)
  legend("topright",title="Signatures", legend=sigs,
         fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  if(!is.null(save_path)) {dev.off()}
}


#========================================#
# Read in the data ####
#========================================#
HDP_folder="."
exposures_file="exposures.csv"
tree_file_path="~/R_work/Prolonged_persistence_of_DNA_lesions/Data/input_data/EM/tree_PX001_2_01_standard_rho01.tree"

tree=read.tree(tree_file_path)

#Create the exposures dataframe from the raw hdp output if not already made
if(file.exists(exposures_file)) {
  exposures_df<-read_csv(exposures_file)
} else {
  #If don't have the exposures df saved, create it from the HDP output as follows
  mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
  mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
  key_table=read.table(paste0(HDP_folder,"/key_table.txt"))
  
  exposures_df=create_exposures_df(HDP_multi = mut_example_multi,trinuc_mut_mat = mutations,key_table=key_table,sep="-")
  write.csv(exposures_df,exposures_file,row.names = F)
}

#========================================#
# Plot the tree ####
#========================================#

plot_sig_tree(tree,sample_ID="KX008_2_01",exposures=exposures_df,show.tip.label=T,cex=0.3)

