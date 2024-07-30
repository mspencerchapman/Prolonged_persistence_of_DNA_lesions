mutations_filtered_file="~/R_work/Prolonged_persistence_of_DNA_lesions/Data/mutations_filtered.tsv"
write.table(dplyr::bind_rows(MAV_mutations,PVV_mutations),file=mutations_filtered_file,quote=F,sep="\t",row.names = F)

##Use KX008 as example for plotting the PVVs and MAVs
ID="KX008_2_01"

for(ID in c(#"KX003_5_01",
            #"KX004_5_01",
            "KX007_2_01",
            #"KX008_2_01",
            "KX009_1_01",
            "KX010_1_01"
            )) {
  
  cat(ID,sep="/n")
  
  details_file_path=paste0(data_dir,"input_data/EM/annotated_mut_set_",ID,"_standard_rho01")
  load(details_file_path)
  details<-filtered_muts$COMB_mats.tree.build$mat
  matrices<-list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR)
  tree<-read.tree(paste0(data_dir,"input_data/EM/tree_",ID,"_standard_rho01.tree"))
  
  #----------------------
  ##Plot the PVVs
  #----------------------
  cat("Plotting PVVs",sep="\n")
  individual_PVVs<-PVV_mutations%>%filter(Sample_ID==ID & Type=="PVV" & !is.na(lesion_node))
  
  PVV_VAFs_across_tree<-lapply(individual_PVVs$mut_ref1,function(mut) {
    mut_VAF_across_tree<-calculate_vaf(NV=matrices$NV[mut,,drop=F],NR=matrices$NR[mut,,drop=F])[,tree$tip.label]
  })%>%bind_rows()
  rownames(PVV_VAFs_across_tree)<-individual_PVVs$mut_ref1
  
  base_cols=colorRampPalette(colors=RColorBrewer::brewer.pal(n=8,name = "Dark2"))(nrow(individual_PVVs))
  hm<-matrix(nrow = nrow(PVV_VAFs_across_tree),ncol=ncol(PVV_VAFs_across_tree),dimnames=dimnames(PVV_VAFs_across_tree))
  for(i in 1:nrow(hm)) {
    vaf_vec=as.integer(round(as.matrix(PVV_VAFs_across_tree)[i,]*100))
    col_scale=colorRampPalette(colors=c("white",base_cols[i]))(1+max(vaf_vec))
    hm[i,]<-col_scale[1+vaf_vec]
  }
  
  pdf(paste0(plots_dir,"PVV_summary_",ID,".pdf"),width=7,height=7)
  tree=plot_tree(tree = tree,cex.label = 0,lwd=ifelse(ID=="KX004_5_01",0.25,0.5),vspace.reserve = 3.5,cex.axis = 0.5,tck=-0.01)
  add_mut_heatmap(tree=tree,heatmap=hm[,tree$tip.label],border="gray",heatmap_bar_height=0.05,cex.label = 0.25,label.cols = base_cols)
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
  
  #----------------------
  ##Plot the MAVs
  #----------------------
  cat("Plotting MAVs",sep="/n")
  
  #Subset the mutation matrix to get the relevant mutations
  individual_MAVs<-MAV_mutations%>%filter(Class!="FAIL"&Sample_ID==ID & Type=="MAV" & !is.na(lesion_node))
  #individual_MAVs<-MAV_mutations%>%filter(Class=="FAIL"&Sample_ID==ID & Type=="MAV")
  
  if(nrow(individual_MAVs>0)) {
    #Create the colour heatmap for plotting
    MAV_cols=RColorBrewer::brewer.pal(n=11,"RdBu")[c(2,9)]
    MAV_cols_across_tree<-lapply(1:nrow(individual_MAVs),function(i) {
      cat(i)
      mut1=individual_MAVs$mut_ref1[i];NV1=matrices$NV[mut1,tree$tip.label]
      mut2=individual_MAVs$mut_ref2[i];NV2=matrices$NV[mut2,tree$tip.label]
      NR=matrices$NR[mut1,tree$tip.label]+matrices$NV[mut2,tree$tip.label]
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
    
    
    #Now do the plotting
    MAV_label_cols=colorRampPalette(colors=RColorBrewer::brewer.pal(n=8,name="Dark2"))(nrow(individual_MAVs))
    
    pdf(paste0(plots_dir,"MAV_summary_",ID,".pdf"),width=7,height=7)
    tree=plot_tree(tree = tree,cex.label = 0,lwd=ifelse(ID=="KX004_5_01",0.25,0.5),vspace.reserve = 1,cex.axis=0.5,tck=-0.01)
    add_mut_heatmap(tree=tree,heatmap=MAV_cols_across_tree[,tree$tip.label,drop=F],border="gray",heatmap_bar_height=0.05,cex.label = 0.25,label.cols = MAV_label_cols)
    temp=sapply(1:nrow(individual_MAVs),function(i) {
      
      #Highlight the lesion node
      info=get_edge_info(tree=tree,details=NULL,node=individual_MAVs$lesion_node[i])
      #points(x=info$x,y=info$yb,pch = 20,col=MAV_label_cols[i])
      plotrix::draw.circle(x=info$x,y=info$yb,radius=2,border=NA,col=MAV_label_cols[i])
      
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

###
individual_PVVs<-mutations%>%filter(Sample_ID==ID & Type=="PVV" & !is.na(lesion_node))
j=14
pdf(paste0(plots_dir,"PVV_no_possible_phylo_examples.pdf"),width=7,height=4)
par(mfrow=c(1,4))
for(j in 12:15) {
  example_PVV<-individual_PVVs$mut_ref1[j]
  
  LN<-individual_PVVs%>%filter(mut_ref1==example_PVV)%>%pull(lesion_node)
  LRN<-individual_PVVs%>%filter(mut_ref1==example_PVV)%>%pull(lesion_repair_node)
  all_LRN_muts<-details%>%filter(node==LRN)%>%pull(mut_ref)
  num=min(15,length(all_LRN_muts))
  example_mut_refs<-details%>%filter(node==LRN)%>%pull(mut_ref)%>%.[1:num]
  
  sub_tree=extract.clade(phy=tree,node=LN)
  sub_tree=plot_tree(sub_tree,cex.label=0,vspace.reserve = 1,plot_axis = F)
  
  hm<-matrix(nrow = 1+num,ncol=length(sub_tree$tip.label),dimnames=list(c(example_PVV,example_mut_refs),sub_tree$tip.label))
  
  base_cols=RColorBrewer::brewer.pal(n=3,name="Dark2")[2:3]
  vaf_vec<-calculate_vaf(NV=matrices$NV[example_PVV,],NR=matrices$NR[example_PVV,])[,sub_tree$tip.label]
  vaf_col_vec=as.integer(round(vaf_vec*100))
  col_scale=colorRampPalette(colors=c("white",base_cols[1]))(1+max(vaf_col_vec))
  hm[1,]<-col_scale[1+vaf_col_vec]
  
  for(i in 1:length(example_mut_refs)) {
    this_mut<-example_mut_refs[i]
    vaf_vec=calculate_vaf(NV=matrices$NV[this_mut,],NR=matrices$NR[this_mut,])[,sub_tree$tip.label]
    vaf_col_vec=as.integer(round(vaf_vec*100))
    col_scale=colorRampPalette(colors=c("white",base_cols[2]))(1+max(vaf_col_vec))
    hm[(i+1),]<-col_scale[1+vaf_col_vec]
  }
  
  add_mut_heatmap(tree=sub_tree,heatmap=hm[,sub_tree$tip.label],border="gray",heatmap_bar_height=0.05,cex.label = 0.5)
}
dev.off()
