#!/software/R-3.6.1/bin/Rscript
library(GenomicRanges)
library(IRanges)
library("Rsamtools")
library("MASS")
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ape)
options(stringsAsFactors = FALSE)

my_working_directory<-getwd()
nsim=50000

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

Phasing_MAV_file_path="output2/Phasing_results_MAVs_all"
Phasing_PVV_file_path="output2/Phasing_results_PVVs_all"

#Set data file paths
data_sets=c("MSC_chemo","NW","MF","EM","KY","PR","MSC_BMT","MSC_fetal","SN")
out_list=lapply(data_sets,function(data_set) {
  print(data_set)
  files=list.files(paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/output2/",data_set,"/"),pattern = "_mut_table.tsv",full.names = T)
  mut_tables_list=lapply(files,read.delim)
  mut_tables_list=lapply(mut_tables_list,function(df) {
    df$colonies_per_negative_subclade<-as.character(df$colonies_per_negative_subclade)
    df$depth_per_negative_subclade<-as.character(df$depth_per_negative_subclade)
    df$Ref<-as.character(df$Ref);df$Alt1<-as.character(df$Alt1);df$Alt2<-as.character(df$Alt2)
    return(df)})
  mut_tables_df=dplyr::bind_rows(mut_tables_list)
  mut_tables_df$data_set=data_set
  return(mut_tables_df)
})
mutations=dplyr::bind_rows(out_list)
mutations$Ref[mutations$Ref=="TRUE"]<-"T"
mutations$Alt1[mutations$Alt1=="TRUE"]<-"T"
mutations$Alt2[mutations$Alt2=="TRUE"]<-"T"

genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

#Define palettes and themes for ggplot2 figures
palette2=c("#E30613","#1D71B8")
palette6=c("#72e5ef", "#074d65", "#3d99ce", "#c257d3", "#d1add5", "#61356e")
#palette8=c(MutationalPatterns:::COLORS7,"blue")
my_theme=theme_classic(base_family="Helvetica",base_size = 7)+theme(text=element_text(size=7,family="Helvetica"),
                                                                    axis.text=element_text(size=7,family="Helvetica"),
                                                                    axis.text.x=element_text(angle = 90),
                                                                    strip.text = element_text(size=7,family="Helvetica"),
                                                                    legend.text = element_text(size=7,family="Helvetica"))


#Import the MAV phasing info, extract a basic summary
load(Phasing_MAV_file_path)
MAV_phasing_summary=sapply(Phasing_results_MAVs,extract_MAV_phasing_summary)
names(MAV_phasing_summary)<-str_split(names(MAV_phasing_summary),pattern="\\.",simplify = T)[,1]
mutations<-left_join(mutations,data.frame(Chrom_pos=names(MAV_phasing_summary),phasing_summary=MAV_phasing_summary),by="Chrom_pos")
mutations$phasing_summary<-sapply(mutations$phasing_summary,function(x) ifelse(x=="No SNPs appear to be heterozygous","Unable to confirm phasing",x))

Sample_IDs=unique(mutations$Sample_ID)
all_MAV_results_list=lapply(data_sets,function(dataset) {
  data_set_MAVs<-mutations%>%filter(Type=="MAV" & data_set==dataset)
  Sample_IDs=unique(data_set_MAVs$Sample_ID)
  data_set_out=lapply(Sample_IDs,function(sample) {
    print(paste("Starting analysis for sample",sample))
    
    #Now do the simulations
    sample_info=get_file_paths_and_project(dataset,Sample_ID=sample)
    tree=read.tree(sample_info$tree_file_path)
    
    #Get the data into the right format
    MAV_data=data_set_MAVs%>%dplyr::filter(Sample_ID==sample)
    MAV_data$branches_prob=mapply(FUN=function(node1,node2){
      prob1=tree$edge.length[tree$edge[,2]==node1]/sum(tree$edge.length)
      prob2=tree$edge.length[tree$edge[,2]==node2]/sum(tree$edge.length)
      score=-log10(prob1*prob2)
    },node1=MAV_data$Node1, node2=MAV_data$Node2)
    MAV_data$data_type="data"
    
    if(dataset=="SN"){
      MAV_output_list=lapply(1:nsim, function(i) {
        nodes=base::sample(x = tree$edge[,2],size = 2,replace = F,prob = tree$edge.length)
        ancestor_1=get_ancestor_node(nodes[1],tree=tree)
        ancestor_2=get_ancestor_node(nodes[2],tree=tree)
        if(ancestor_1%in%c(ancestor_2,get_all_node_children(ancestor_2,tree))|ancestor_2%in%c(ancestor_1,get_all_node_children(ancestor_1,tree))){
          ROOT=tree$edge[1,1]
          if(ancestor_1==ROOT&!ancestor_2%in%c(nodes[1],get_all_node_children(nodes[1],tree))|ancestor_2==ROOT&!ancestor_1%in%c(nodes[2],get_all_node_children(nodes[2],tree))){
            Class<-"FAIL"
          } else {
            Class<-"PASS"
          }
        } else {
          Class<-"FAIL"
        }
        return(data.frame(Class=Class,
                          Node1=nodes[1],
                          Node2=nodes[2],
                          no_of_cell_divisions=NA,
                          lesion_node=NA,
                          lesion_timing=NA,
                          lesion_repair_node=NA,
                          lesion_repair_timing=NA,
                          lesion_duration=NA))
      })
    }else {
      MAV_output_list=lapply(1:nsim, function(i) {
        nodes=base::sample(x = tree$edge[,2],size = 2,replace = F,prob = tree$edge.length)
        ancestor_1=get_ancestor_node(nodes[1],tree=tree)
        ancestor_2=get_ancestor_node(nodes[2],tree=tree)
        if((nodes[1]!=nodes[2])&(ancestor_1==ancestor_2)){
          lesion_duration=0
          lesion_node=lesion_repair_node=ancestor_1
          minimum_cell_divisions=1
          lesion_node_res=list(initial_lesion_node=lesion_node,Filter="PASS",Class="simple")
        }else if(ancestor_1%in%get_all_node_children(ancestor_2,tree)|ancestor_2%in%get_all_node_children(ancestor_1,tree)) {
          lesion_node_res=find_MAV_lesion_node(node1=nodes[1],node2=nodes[2],tree)
          lesion_repair_node=ifelse(nodeheight(tree,ancestor_1)>nodeheight(tree,ancestor_2),ancestor_1,ancestor_2)
          ancestor_node_heights=sapply(c(ancestor_1,ancestor_2),function(node) nodeheight(tree=tree,node=node))
          lesion_duration=max(ancestor_node_heights)-nodeheight(tree=tree,node=lesion_node_res$initial_lesion_node)
          if(lesion_node_res$initial_lesion_node==lesion_repair_node){
            lesion_node_res$Class<-"simple"
            minimum_cell_divisions=1
            lesion_duration=0
          } else {
            cell_divs=1
            repeat{an_node=get_ancestor_node(lesion_repair_node,tree,cell_divs);if(an_node==lesion_node_res$initial_lesion_node){break}else{cell_divs=cell_divs+1}}
            minimum_cell_divisions=cell_divs+1
          }
        } else {
          lesion_node=lesion_repair_node=NA
          lesion_duration=NA
          minimum_cell_divisions=NA
          lesion_node_res=list(initial_lesion_node=NA,Filter="FAIL",Class="FAIL")
        }
        return(data.frame(Class=lesion_node_res$Class,
                          Node1=nodes[1],
                          Node2=nodes[2],
                          no_of_cell_divisions=minimum_cell_divisions,
                          lesion_node=lesion_node_res$initial_lesion_node,
                          lesion_timing=ifelse(is.na(lesion_node_res$initial_lesion_node),NA,nodeheight(tree,lesion_node_res$initial_lesion_node)),
                          lesion_repair_node=lesion_repair_node,
                          lesion_repair_timing=ifelse(is.na(lesion_repair_node),NA,nodeheight(tree=tree,node=lesion_repair_node)),
                          lesion_duration=lesion_duration))
      })
    }
    
    MAV_out=bind_rows(MAV_output_list)
    
    #Calculate the "branches_prob": -log10 that a mutation occurring twice by chance would be on these two branches
    MAV_out$branches_prob=mapply(FUN=function(node1,node2){
      prob1=tree$edge.length[tree$edge[,2]==node1]/sum(tree$edge.length)
      prob2=tree$edge.length[tree$edge[,2]==node2]/sum(tree$edge.length)
      score=-log10(prob1*prob2)
    },node1=MAV_out$Node1, node2=MAV_out$Node2)
    MAV_out$data_type="simulation"
    MAV_out$Sample_ID=sample
    MAV_out$data_set=dataset

    #Now generate the plots
    MAV_data_removed=MAV_data%>%filter(Class=="removed")%>%mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Unable to confirm phasing","Non-matching phasing confirmed")))
    
    myColors <- c("#377EB8","grey","#E41A1C")
    names(myColors) <- levels(MAV_data_removed$phasing_summary)
    # 
    # p1<-MAV_out%>%
    #   mutate(data_type="simulation")%>%
    #   #bind_rows(MAV_data)%>%
    #   filter(Class=="removed")%>%
    #   ggplot(aes(y=branches_prob,x=lesion_duration)) +
    #   stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white",size=0.15) +
    #   facet_grid(rows=vars(Class))+
    #   my_theme+
    #   #theme(legend.position = "none")+
    #   geom_point(data=MAV_data_removed,aes(col=phasing_summary),size=3,alpha=0.5)+
    #   scale_color_manual(name = "phasing_summary",values = myColors)+
    #   scale_fill_continuous(guide='none')+
    #   labs(title=paste0(sample,": Data vs simulation\n Distribution of removed MAVS"),
    #        x="Lesion duration",
    #        y="-log10 probability of branches involved by chance")
    # 
    # p2<-MAV_out%>%
    #   mutate(data_type="simulation",phasing_summary="Unable to confirm phasing")%>%
    #   bind_rows(MAV_data)%>%
    #   mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Unable to confirm phasing","Non-matching phasing confirmed")))%>%
    #   filter(Class=="simple")%>%
    #   ggplot(aes(x=data_type,y=branches_prob,fill=data_type))+
    #   #geom_boxplot(col="black")+
    #   geom_violin()+
    #   geom_jitter(aes(col=phasing_summary),width=0.1,size=2,alpha=0.75)+
    #   scale_color_manual(name = "phasing_summary",values = myColors)+
    #   scale_fill_discrete(guide='none')+
    #   labs(title=paste0(sample,": Data vs simulation\nProbability of simple MAV branches by chance"),
    #        col="Data type",
    #        fill="Data type",
    #        x="",
    #        y="-log10 probability of branches involved by chance")+
    #   my_theme
    # 
    # MAV_data_props<-MAV_data%>%mutate(Class=factor(Class,levels=c("simple","removed","FAIL")))%>%count(Class)%>%complete(Class,fill=list(n=0))%>%mutate(prop=n/sum(n),data_type="data")
    # MAV_out_props<-MAV_out%>%mutate(Class=factor(Class,levels=c("simple","removed","FAIL")))%>%count(Class)%>%complete(Class,fill=list(n=0))%>%mutate(prop=n/sum(n),data_type="simulation")
    # 
    # p3<-bind_rows(MAV_data_props,MAV_out_props)%>%
    #   ggplot(aes(x=Class,y=prop,fill=data_type))+
    #   geom_bar(stat="identity",position="dodge",col="black",size=0.15)+
    #   my_theme+
    #   labs(fill="Data type",
    #        x="MAV class",
    #        y="Proportion of total MAVs",
    #        title=paste0(sample,": Data vs simulation\n MAV class proportions"))
    # 
    # p4<-bind_rows(MAV_data_props,MAV_out_props)%>%
    #   ggplot(aes(x=data_type,y=prop,fill=Class))+
    #   geom_bar(stat="identity",position="stack",col="black",size=0.15)+
    #   my_theme+
    #   scale_fill_brewer(palette = "Set2")+
    #   labs(fill="MAV Class",
    #        x="Data type",
    #        y="Proportion of total MAVs",
    #        title=paste0(sample,": Data vs simulation\n MAV class proportions"))
    # 
    # 
    # plot_comb<-arrangeGrob(p4,p1,p2,ncol=3,widths = c(1,3,2))
    # 
    MAV_comb=dplyr::bind_rows(MAV_data,MAV_out)
    return(MAV_comb)
    #ggsave(filename=paste0("output2/MAV_vs_sim_plots/",sample,"_MAV_vs_sim_plots.pdf"),plot=plot_comb,device="pdf",width=12,height=3)
  })
  data_set_df=dplyr::bind_rows(data_set_out)
  return(data_set_df)
})

all_MAV_results_df=dplyr::bind_rows(all_MAV_results_list)
write.table(all_MAV_results_df,file="output2/MAV_sim_res/MAV_sim_results.tsv",quote=F,sep="\t",row.names=F)
