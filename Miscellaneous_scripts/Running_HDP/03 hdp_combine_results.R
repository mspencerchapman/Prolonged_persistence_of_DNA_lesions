#!/software/R-3.6.1/bin/Rscript
# HDP Flow III: Combine results
# Tim Coorens, Feb 2020
#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("devtools")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

options(stringsAsFactors = F)

#========================================#
# Set paths for running ####
#========================================#

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
HDP_folder=paste0(root_dir,"/data/HDP")
setwd(HDP_folder)

#========================================#
# Read in the output from the individual chains & combine ####
#========================================#

chlist <- vector("list", 20)
for (i in 1:20){
  if(file.exists(paste0("hdp_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

mut_example_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi) #This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi,"HDP_multi_chain.Rdata")

#========================================#
# Diagnostic plots ####
#========================================#

# Muts attributed plot ----
pdf("muts_attributed.pdf")
plot_comp_size(mut_example_multi, bty="L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

# 96-profile signature plots ----
for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_",i,".pdf"),width=12,height=4)
  
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

# Exposure plot ----
plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"),
                      incl_nonsig=TRUE, cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')

# Signature attribution plot ----
pdf("signature_attribution.pdf",width=10,height=8)
key_table=read.table("key_table.txt")
plot_dp_comp_exposure(mut_example_multi, dpindices=(length(unique(key_table$Patient))+2):length(mut_example_multi@comp_dp_counts), incl_nonsig = T,ylab_exp = 'Signature exposure', leg.title = 'Signature',
                      col=c(RColorBrewer::brewer.pal(12, "Set3"),"magenta","firebrick"))
dev.off()
