library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/lesion_segregation"

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
#sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
resave_plots=F
plots_dir="PVV_MAV_plots/"

mutations_file="mutations.tsv"
mutations_filtered_file="mutations_filtered.tsv"
mutations=read.delim(mutations_filtered_file,stringsAsFactors = F)
all_IDs=unique(mutations%>%filter(data_set!="SN")%>%pull(Sample_ID)) #Can't do this for the liver dataset as these are inferred clones
lesion_seg_input_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/"

#----------------------
###Do overall summaries
#----------------------

for(ID in all_IDs[13:49]) {
  
  cat(ID,sep="/n")
  data_set=mutations%>%filter(Sample_ID==ID)%>%pull(data_set)%>%.[1]
  info=get_file_paths_and_project(dataset=data_set,Sample_ID = ID,input_data_dir = lesion_seg_input_dir)
  load(info$filtered_muts_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  matrices<-list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR)
  tree<-read.tree(info$tree_file_path)
  tree_no_ancestral<-drop.tip(tree,tip="Ancestral")
  #----------------------
  ##Plot the PVVs
  #----------------------
  cat("Plotting PVVs",sep="\n")
  individual_PVVs<-mutations%>%
    filter(Sample_ID==ID & Type=="PVV" & !is.na(lesion_node))%>%
    arrange(lesion_timing)
  
  if(nrow(individual_PVVs)>0) {
    PVV_VAFs_across_tree<-lapply(individual_PVVs$mut_ref1,function(mut) {
      mut_VAF_across_tree<-calculate_vaf(NV=matrices$NV[mut,,drop=F],NR=matrices$NR[mut,,drop=F])[,tree_no_ancestral$tip.label]
    })%>%bind_rows()
    rownames(PVV_VAFs_across_tree)<-individual_PVVs$mut_ref1
    
    base_cols=colorRampPalette(colors=RColorBrewer::brewer.pal(n=8,name = "Dark2"))(nrow(individual_PVVs))
    hm<-matrix(nrow = nrow(PVV_VAFs_across_tree),ncol=ncol(PVV_VAFs_across_tree),dimnames=dimnames(PVV_VAFs_across_tree))
    for(i in 1:nrow(hm)) {
      vaf_vec=as.integer(round(as.matrix(PVV_VAFs_across_tree)[i,]*100))
      col_scale=colorRampPalette(colors=c("white",base_cols[i]))(1+max(vaf_vec))
      hm[i,]<-col_scale[1+vaf_vec]
    }
    if("Ancestral"%in%tree$tip.label) {
      hm<-cbind(hm,matrix(rep("#FFFFFF",nrow(hm)),ncol = 1,dimnames = list(rownames(hm),"Ancestral")))
    }
    
    pdf(paste0(plots_dir,"PVV_summary_",ID,".pdf"),width=7,height=7)
    tree=plot_tree(tree = tree,cex.label = 0,lwd=ifelse(ID=="KX004_5_01",0.25,0.5),vspace.reserve = 3.5,cex.axis = 0.5,tck=-0.01)
    add_mut_heatmap(tree=tree,heatmap=hm[,tree$tip.label,drop=F],border="gray",heatmap_bar_height=0.05,cex.label = 0.25,label.cols = base_cols)
    temp=sapply(1:nrow(PVV_VAFs_across_tree),function(i) {
      LN_ancestors<-get_ancestral_nodes(node = individual_PVVs$lesion_node[i],edge=tree$edge)
      LRN_ancestors<-get_ancestral_nodes(node = individual_PVVs$lesion_repair_node[i],edge=tree$edge)
      lesion_path_nodes=LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
      temp=add_annotation(tree=tree,
                          annot_function = highlight_nodes,
                          nodes=lesion_path_nodes,
                          col=base_cols[i],
                          lwd=ifelse(ID=="KX004_5_01",1,2))
      
    })
    dev.off()
    
  }
  
  #----------------------
  ##Plot the MAVs
  #----------------------
  cat("Plotting MAVs",sep="/n")
  
  #Subset the mutation matrix to get the relevant mutations
  individual_MAVs<-mutations%>%
    filter(Class!="FAIL" & Sample_ID==ID & Type=="MAV" & Mut_type1=="SNV" & Mut_type2=="SNV" & !is.na(lesion_node))%>%
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
        #points(x=info$x,y=info$yb,pch = 20,col=MAV_label_cols[i])
        plotrix::draw.circle(x=info$x,y=info$yb,radius=0.25,border=NA,col=MAV_label_cols[i])
        
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

#----------------------
###Do summaries only of early mutations in early embryogenesis
#----------------------


for(ID in all_IDs[13:49]) {
  
  cat(ID,sep="/n")
  data_set=mutations%>%filter(Sample_ID==ID)%>%pull(data_set)%>%.[1]
  info=get_file_paths_and_project(dataset=data_set,Sample_ID = ID)
  #details_file_path=paste0(data_dir,"input_data/EM/annotated_mut_set_",ID,"_standard_rho01")
  load(info$filtered_muts_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  matrices<-list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR)
  #tree<-read.tree(paste0(data_dir,"input_data/EM/tree_",ID,"_standard_rho01.tree"))
  tree<-read.tree(info$tree_file_path)
  tree_no_ancestral<-drop.tip(tree,tip="Ancestral")
  #----------------------
  ##Plot the PVVs
  #----------------------
  cat("Plotting PVVs",sep="\n")
  individual_PVVs<-mutations%>%
    filter(Sample_ID==ID & Type=="PVV" & !is.na(lesion_node))%>%
    filter(lesion_timing<20)%>%
    arrange(lesion_timing)
  
  if(nrow(individual_PVVs)>0) {
    PVV_VAFs_across_tree<-lapply(individual_PVVs$mut_ref1,function(mut) {
      mut_VAF_across_tree<-calculate_vaf(NV=matrices$NV[mut,,drop=F],NR=matrices$NR[mut,,drop=F])[,tree_no_ancestral$tip.label]
    })%>%bind_rows()
    rownames(PVV_VAFs_across_tree)<-individual_PVVs$mut_ref1
    
    base_cols=colorRampPalette(colors=RColorBrewer::brewer.pal(n=8,name = "Dark2"))(nrow(individual_PVVs))
    hm<-matrix(nrow = nrow(PVV_VAFs_across_tree),ncol=ncol(PVV_VAFs_across_tree),dimnames=dimnames(PVV_VAFs_across_tree))
    for(i in 1:nrow(hm)) {
      vaf_vec=as.integer(round(as.matrix(PVV_VAFs_across_tree)[i,]*100))
      col_scale=colorRampPalette(colors=c("white",base_cols[i]))(1+max(vaf_vec))
      hm[i,]<-col_scale[1+vaf_vec]
    }
    if("Ancestral"%in%tree$tip.label) {
      hm<-cbind(hm,matrix(rep("#FFFFFF",nrow(hm)),ncol = 1,dimnames = list(rownames(hm),"Ancestral")))
    }
    
    pdf(paste0(plots_dir,"Embryonic_PVV_summary_",ID,".pdf"),width=7,height=7)
    tree=plot_tree(tree = squash_tree(tree,cut_off = 20),cex.label = 0,lwd=ifelse(ID=="KX004_5_01",0.25,0.5),vspace.reserve = 3.5,cex.axis = 0.5,tck=-0.01)
    add_mut_heatmap(tree=tree,heatmap=hm[,tree$tip.label,drop=F],border="gray",heatmap_bar_height=0.05,cex.label = 0.25,label.cols = base_cols)
    temp=sapply(1:nrow(PVV_VAFs_across_tree),function(i) {
      LN_ancestors<-get_ancestral_nodes(node = individual_PVVs$lesion_node[i],edge=tree$edge)
      LRN_ancestors<-get_ancestral_nodes(node = individual_PVVs$lesion_repair_node[i],edge=tree$edge)
      lesion_path_nodes=LRN_ancestors[!LRN_ancestors%in%LN_ancestors]
      temp=add_annotation(tree=tree,
                          annot_function = highlight_nodes,
                          nodes=lesion_path_nodes,
                          col=base_cols[i],
                          lwd=ifelse(ID=="KX004_5_01",1,2))
      
    })
    dev.off()
    
  }
  
  #----------------------
  ##Plot the MAVs
  #----------------------
  cat("Plotting MAVs",sep="/n")
  
  #Subset the mutation matrix to get the relevant mutations
  individual_MAVs<-mutations%>%
    filter(Class!="FAIL" & Sample_ID==ID & Type=="MAV" & Mut_type1=="SNV" & Mut_type2=="SNV" & !is.na(lesion_node))%>%
    filter(lesion_timing<20)%>%
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
    
    pdf(paste0(plots_dir,"Embryonic_MAV_summary_",ID,".pdf"),width=7,height=7)
    tree=plot_tree(tree = squash_tree(tree,cut_off = 20),cex.label = 0,lwd=ifelse(ID=="KX004_5_01",0.25,0.5),vspace.reserve = 1,cex.axis=0.5,tck=-0.01)
    add_mut_heatmap(tree=tree,heatmap=MAV_cols_across_tree[,tree$tip.label,drop=F],border="gray",heatmap_bar_height=0.05,cex.label = 0.25,label.cols = MAV_label_cols)
    if(data_set!="SN"){
      temp=sapply(1:nrow(individual_MAVs),function(i) {
        
        #Highlight the lesion node
        info=get_edge_info(tree=tree,details=NULL,node=individual_MAVs$lesion_node[i])
        #points(x=info$x,y=info$yb,pch = 20,col=MAV_label_cols[i])
        plotrix::draw.circle(x=info$x,y=info$yb,radius=0.25,border=NA,col=MAV_label_cols[i])
        
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


#----------------------
#Plot the alternative genotype
#----------------------

#Adjust the add_heatmap function
add_mut_heatmap=function(tree,heatmap,heatvals=NULL,border="white",heatmap_bar_height=0.05,cex.label=2,label.cols=NA){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-heatmap_bar_height*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      if(!is.na(label.cols[1])) {
        text(labels[i],x=0.5,y=0.5*(top+bot),pos = 2,cex = cex.label,col=label.cols[i])
      } else {
        text(labels[i],x=0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
      }
      
    }
    top=bot
  }
  tree
}

source("/lustre/scratch126/casm/team154pc/ms56/my_programs/Prolonged_persistence_functions.R")
mutations<-read.delim("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/mutations_filtered.tsv")
node_confidence_df<-read.delim("/lustre/scratch126/casm/team154pc/ms56/Emily_benchmarking/KX003/vaf_filtered/tree_bootstraps/bootstrap_retained_KX003.tsv")
this_sample_PVVs<-mutations%>%filter(Sample_ID=="KX003_5_01" & Type == "PVV" & Class=="PASS")
this_sample_PVVs%>%
  mutate(LN_confidence=sapply(lesion_node,function(this_node) {node_confidence_df%>%filter(node==this_node)%>%pull(retained)}))%>%
  mutate(LRN_confidence=sapply(lesion_repair_node,function(this_node) {node_confidence_df%>%filter(node==this_node)%>%pull(retained)}))%>%
  dplyr::select(Chrom_pos,lesion_node,lesion_repair_node,LN_confidence,LRN_confidence)

load("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/EM/annotated_mut_set_KX003_5_01_standard_rho01")
tree=ape::read.tree("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/EM/tree_KX003_5_01_standard_rho01.tree")
details<-filtered_muts$COMB_mats.tree.build$mat
LN=472
LRN=473
sub_tree<-extract.clade(phy=tree,node=LN)
PVVs<-this_sample_PVVs%>%filter(lesion_node==LN)%>%pull(mut_ref1)
confirmatory_mutations<-details%>%filter(node==LRN)%>%pull(mut_ref)
combined_muts=c(PVVs,confirmatory_mutations)

create_vaf_mat=function(NV,NR) {NR[NR==0]<-1;return(NV/NR)}
vaf_mat=create_vaf_mat(filtered_muts$COMB_mats.tree.build$NV,filtered_muts$COMB_mats.tree.build$NR)
hm=matrix(NA,nrow=length(combined_muts),ncol=length(sub_tree$tip.label),dimnames = list(combined_muts,sub_tree$tip.label))
PVV_cols=colorRampPalette(colors=c("white","#FF7F00"))(101)
confirmatory_mut_cols=colorRampPalette(colors=c("white","#377EB8"))(101)
for(i in 1:nrow(hm)) {
  if(i %in% 1:length(PVVs)) {palette<-PVV_cols} else {palette <-confirmatory_mut_cols}
  for(j in 1:ncol(hm)) {
    vaf<-vaf_mat[combined_muts[i],sub_tree$tip.label[j]]
    hm[i,j]<-palette[1+round(100*vaf)]
  }
}

pdf("plots/KX003_alternative_tree_comparison.pdf",width=5,height=3)
sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve=0.5)
add_mut_heatmap(tree=sub_tree,heatmap=hm,heatmap_bar_height=0.05,cex.label=0.3)
dev.off()

filtered_muts$COMB_mats.tree.build$NV[combined_muts,sub_tree$tip.label]

