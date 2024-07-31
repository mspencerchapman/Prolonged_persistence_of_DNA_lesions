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
# START THE INDEPENDENT MAV SIMULATIONS ####
#========================================#
#Either import the existing simulation output, or can re-run if desired (but takes a while)

Sample_IDs=unique(mutations$Sample_ID)
MAV_simulation_results_file=paste0(root_dir,"/Data/simulation_results/MAV_sim_results.tsv")

if(file.exists(MAV_simulation_results_file)) {
  MAV_sim_dat<-data.table::fread(MAV_simulation_results_file)%>%filter(Sample_ID%in%ref_table$Sample_ID)
} else {
  #To rerun the simulations
  nsim=50000
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
      
      MAV_comb=dplyr::bind_rows(MAV_data,MAV_out)
      return(MAV_comb)
    })
    data_set_df=dplyr::bind_rows(data_set_out)
    return(data_set_df)
  })
  
  all_MAV_results_df=dplyr::bind_rows(all_MAV_results_list)
  write.table(all_MAV_results_df,file=MAV_simulation_results_file,quote=F,sep="\t",row.names=F)
}

#========================================#
# Now do the independent MAV simulation plotting ####
#========================================#

## Generate Extended Data Fig. 3a ------
# plot PASS MAV for the SN data_set

rename_LIVER_MAV_vec=c("Related","Unrelated")
names(rename_LIVER_MAV_vec)=c("PASS","FAIL")

MAVsims.plot.SN<-MAV_sim_dat%>%
  filter(data_set=="SN")%>%
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,PASS=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+PASS),PASS_p=PASS/(FAIL+PASS))%>%
  dplyr::select(-PASS,-FAIL)%>%
  gather(-Sample_ID,-data_type,key="Class",value="Proportion")%>%
  mutate(Class=gsub("_p","",Class))%>%
  mutate(Class=factor(rename_LIVER_MAV_vec[Class],levels=rename_LIVER_MAV_vec))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(data_type=rename_data_type_vec[data_type],cat="Liver")%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values = c("#66C2A5","lightgray"))+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,0.5)) +
  my_theme+
  facet_grid(cat~data_type,scales="free",space = "free")+
  coord_flip()+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5))+
  labs(x="",y="Proportion",fill="")+
  theme(legend.position = "bottom",legend.margin=margin(0,0,0,0),legend.box.margin = margin(0,0,0,0,unit="mm"))

ggsave(MAVsims.plot.SN,filename = paste0(plots_dir,"ExtDatFig3a.pdf"),width=3.5,height=4.5)