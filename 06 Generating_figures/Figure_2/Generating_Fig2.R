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

#Do the first plot of class proportions across all the samples
MAV_props_df<-MAV_sim_dat%>%
  filter(data_set!="PR"&data_set!="SN"&!data_set=="MSC_chemo")%>%
  mutate(Class=ifelse(Class=="removed"&no_of_cell_divisions>2&data_type=="simulation","FAIL",Class))%>%
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,removed=0,simple=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+removed+simple),removed_p=removed/(FAIL+removed+simple),simple_p=simple/(FAIL+removed+simple))%>%
  dplyr::select(-removed,-FAIL,-simple)%>%
  gather(-Sample_ID,-data_type,key="Class",value="Proportion")%>%
  mutate(Class=gsub("_p","",Class))

rename_MAV_vec=c("Simple","Separated","Unrelated")
names(rename_MAV_vec)=c("simple","removed","FAIL")
rename_data_type_vec=c("Observed data","Simulations")
names(rename_data_type_vec)=c("data","simulation")
rename_cat_vec=c("Adult HSPC","C","F","Bronchial")
names(rename_cat_vec)=c("Adult_HSPC","Chemo_HSPC","Foetal_HSPC","Bronchial")

relationship_col_vec=c("#66C2A5","#FC8D62","lightgray")
names(relationship_col_vec)<-rename_MAV_vec

## Generate Figure 2a ------
MAVsims.p4<-MAV_props_df%>%
  left_join(summary_table_df%>%dplyr::select(Sample_ID,cat))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Class=factor(rename_MAV_vec[Class],levels=rev(rename_MAV_vec)))%>%
  mutate(data_type=rename_data_type_vec[data_type],cat=factor(rename_cat_vec[cat],levels=rename_cat_vec))%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(cat~data_type,scales="free",space = "free")+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 10),values = relationship_col_vec)+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,0.5)) +
  my_theme+
  theme(axis.text.x = element_text(size=6,angle=90,vjust = +0.5))+
  labs(x="")+
  coord_flip()+
  theme(legend.position = "bottom",legend.margin=margin(0,0,0,0),legend.box.margin = margin(0,0,0,0,unit="mm"))

ggsave(MAVsims.p4,filename = paste0(plots_dir,"Fig2a.pdf"),width=3.5,height=5.3)


#Now plot PASS MAV for the SN data_set

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


ggsave(MAVsims.plot.SN,filename = paste0(plots_dir,"ExtDatFig3a.pdf"),width=3.5,height=5.3)

#========================================#
#Calculate a WEIGHTED combined MAV proportions across samples ####
#========================================#
#Plot overall proportions of classes across samples, data vs simulation
#To be an accurate representation of the data, simulation proportions should be weighted by the total number of MAVs recorded in the data in that sample

## 1. calculate the raw overall proportions in each sample - simulations vs data ----
MAV_prop_summary<-MAV_sim_dat%>%
  filter(data_set!="PR"&data_set!="SN"&!data_set=="MSC_chemo")%>%
  mutate(n_neg=ifelse(Class=="removed"&data_type=="data",sapply(colonies_per_negative_subclade,function(x) length(unlist(strsplit(x[!is.na(x)],",")))),NA))%>%
  mutate(Class=ifelse(Class=="removed"&n_neg>1&data_type=="data","FAIL",Class))%>% #Fail the MAVs in the data with more than 1 x negative subclade
  mutate(Class=ifelse(Class=="removed"&no_of_cell_divisions>2&data_type=="simulation","FAIL",Class))%>% #Fail the removed MAVs in the simulations that would have more than 1 x negative subclade
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,removed=0,simple=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+removed+simple),removed_p=removed/(FAIL+removed+simple),simple_p=simple/(FAIL+removed+simple))%>%
  dplyr::select(-removed,-FAIL,-simple)

## 2. Calculate weights, based on the absolute number of MAVs found in each sample (regardless of class) - otherwise simulations from all samples will be weighted equally even though contribute very few MAVs to the actual data set ----
MAV_weights=MAV_sim_dat%>%filter(data_set!="PR"&data_set!="SN"&!data_set=="MSC_chemo")%>%filter(data_type=="data")%>%group_by(Sample_ID)%>%summarise(n=n())

## 3. Calculate the proportions of different classes in the simulations - including the weighting ----
MAV_sim_props<-MAV_prop_summary%>%
  filter(data_type=="simulation")%>%
  group_by(.)%>%
  summarise(data_type="simulation",FAIL=weighted.mean(FAIL_p,w=MAV_weights$n),removed=weighted.mean(removed_p,w=MAV_weights$n),simple=weighted.mean(simple_p,w=MAV_weights$n))%>%
  gather(-data_type,key="Class",value="Proportion")

## 4. Calculate the proportions of different classes in the data - including the weighting ----
MAV_dat_props<-MAV_sim_dat%>%
  filter(data_set!="PR"&data_set!="SN"&!data_set=="MSC_chemo")%>%
  mutate(n_neg=ifelse(Class=="removed"&data_type=="data",sapply(colonies_per_negative_subclade,function(x) length(unlist(strsplit(x[!is.na(x)],",")))),NA))%>%
  mutate(Class=ifelse(Class=="removed"&n_neg>1&data_type=="data","FAIL",Class))%>%
  mutate(Class=ifelse(Class=="removed"&no_of_cell_divisions>2&data_type=="simulation","FAIL",Class))%>%
  group_by(Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,removed=0,simple=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+removed+simple),removed_p=removed/(FAIL+removed+simple),simple_p=simple/(FAIL+removed+simple))%>%
  dplyr::select(-removed,-FAIL,-simple)%>%
  gather(-data_type,key="Class",value="Proportion")%>%
  mutate(Class=gsub("_p","",Class))%>%
  filter(data_type=="data")

## Generate Figure 2b ------
MAVsims.p5<-rbind(MAV_dat_props,MAV_sim_props)%>%
  mutate(Class=factor(rename_MAV_vec[Class],levels=rename_MAV_vec))%>%
  mutate(data_type=str_wrap(rename_data_type_vec[data_type],width=10))%>%
  ggplot(aes(x=data_type,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 10),values = relationship_col_vec)+
  my_theme+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.text.x = element_text(angle=90),vjust = +0.5)+
  labs(x="",y="Proportions",fill="")
ggsave(MAVsims.p5,filename = paste0(plots_dir,"Fig2b.pdf"),width=1.7,height=2)

#========================================#
# Review the phasing results of the MAVs ####
#========================================#
MAV_mutations<-mutations%>%
  filter(Type=="MAV")%>%
  mutate(phasing_summary=factor(phasing_summary,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

myPhasingColors <- c("#478EB0","#A35563","grey")
names(myPhasingColors) <- levels(MAV_mutations$phasing_summary)

## Generate Figure 2c ------
rename_categories=c("Simple","Separated","Unrelated\n(controls)")
names(rename_categories)=c("simple","removed","FAIL")

rename_phasing=c("Matching phasing","Opposite phasing","Unable to confirm")
names(rename_phasing)=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")

myPhasingColors <- c("#478EB0","#A35563","grey")
names(myPhasingColors) <- rename_phasing

p.MAV.3.2<-MAV_mutations%>%
  filter(data_set!="SN")%>%
  mutate(Original_Class=ifelse(n_neg>1,"FAIL",Original_Class))%>%
  #mutate(phasing_summary=factor(phasing_summary))%>%
  group_by(Original_Class,phasing_summary)%>%
  dplyr::count(Original_Class,phasing_summary)%>%
  ungroup()%>%
  complete(Original_Class,phasing_summary,fill=list(n=0))%>%
  mutate(Original_Class=factor(rename_categories[Original_Class],levels=rename_categories))%>%
  mutate(phasing_summary=factor(rename_phasing[phasing_summary],levels=rename_phasing))%>%
  filter(phasing_summary%in%rename_phasing[1:2])%>%
  ggplot(aes(x=Original_Class,y=n,fill=phasing_summary))+
  geom_bar(stat="identity",position="dodge",col="black",linewidth=0.2)+
  scale_fill_manual(name = "MAV phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  guides(fill=guide_legend(title="MAV phasing",position="top",ncol=1))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=45,vjust=0.9,hjust=0.9),axis.title.x = element_blank(),legend.title = element_text(face="bold"))+
  labs(x="",
       y="Number of mutations")

ggsave(p.MAV.3.2,filename=paste0(plots_dir,"Fig2c.pdf"),width=2.5,height=2)

#========================================#
# Review the phasing results of the PVVs ####
#========================================#
PVV_mutations<-mutations%>%
  filter(Type=="PVV")%>%
  mutate(PVV_pos_clade_phasing=factor(PVV_pos_clade_phasing,levels=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")))

myPhasingColors <- c("#478EB0","#A35563","grey")
names(myPhasingColors) <- levels(PVV_mutations$PVV_pos_clade_phasing)

## Generate Figure 2c ------
rename_categories=c("Simple","Separated","Unrelated\n(controls)")
names(rename_categories)=c("simple","removed","FAIL")

rename_phasing=c("Matching phasing","Opposite phasing","Unable to confirm")
names(rename_phasing)=c("Same phasing confirmed","Non-matching phasing confirmed","Unable to confirm phasing")

myPhasingColors <- c("#478EB0","#A35563","grey")
names(myPhasingColors) <- rename_phasing

p.PVV<-PVV_mutations%>%
  filter(data_set!="SN")%>%
  mutate(Original_Class=ifelse(is.na(lesion_node)|ASCAT_result=="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected","FAIL","PASS"))%>%
  mutate(Class=ifelse(max_neg_clade_depth<13|is.na(max_neg_clade_depth)|is.na(lesion_node),"FAIL",Original_Class))%>%
  mutate(Class=ifelse(ASCAT_result!="LOH"& basic_result!="LOH",Class,"FAIL"))%>%
  group_by(Class,PVV_pos_clade_phasing)%>%
  dplyr::count(cat,Class,PVV_pos_clade_phasing)%>%
  ungroup()%>%
  complete(cat,Class,PVV_pos_clade_phasing,fill=list(n=0))%>%
  mutate(PVV_pos_clade_phasing=factor(rename_phasing[PVV_pos_clade_phasing],levels=rename_phasing))%>%
  filter(PVV_pos_clade_phasing%in%rename_phasing[1:2])%>%
  ggplot(aes(x=cat,y=n,fill=PVV_pos_clade_phasing))+
  geom_bar(stat="identity",position="dodge",col="black",linewidth=0.2)+
  scale_fill_manual(name = "PVV phasing",values = myPhasingColors,labels = function(x) str_wrap(x, width = 18))+
  guides(fill=guide_legend(title="PVV phasing",position="top",ncol=1))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=45,vjust=0.9,hjust=0.9),axis.title.x = element_blank(),legend.title = element_text(face="bold"))+
  labs(x="",
       y="Number of mutations")

ggsave(p.PVV,filename=paste0(plots_dir,"Fig2d.pdf"),width=2.5,height=3)

#========================================#
# Review the ASCAT LOH results of the PVVs ####
#========================================#

## Generate Figure 2e ------
myLOHColors<-c(rev(RColorBrewer::brewer.pal(3,"Accent")),"grey")
names(myLOHColors)<-c("LOH","Unable to confirm","No LOH","Homozygous in all clades")

p.PVV.LOH<-PVV_mutations%>%
  mutate(Original_Class=ifelse(is.na(lesion_node)|ASCAT_result=="Need at least 2 pure positive AND 1 pure negative clade to be within the lesion node: one or both not detected","FAIL","PASS"))%>%
  mutate(Class=ifelse(max_neg_clade_depth<13|is.na(max_neg_clade_depth)|is.na(lesion_node),"FAIL",Original_Class))%>%
  mutate(ASCAT_result=factor(ASCAT_result))%>%filter(Class!="FAIL")%>%
  mutate(cat=stringr::str_replace(cat,pattern = "_",replacement = "\n"))%>%
  group_by(ASCAT_result)%>%
  dplyr::count(ASCAT_result,cat)%>%
  ggplot(aes(x=cat,y=n,fill=ASCAT_result))+
  geom_bar(stat="identity",position="stack",col=NA,size=0.2)+
  scale_fill_manual(name = "ASCAT LOH result",values=myLOHColors,labels = function(x) str_wrap(x, width = 9))+
  guides(fill=guide_legend(title = "Copy number at PVVs",position = "inside"))+
  theme_classic()+
  my_theme+
  theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=45,hjust=0.8,vjust=0.8),legend.position.inside=c(0.7,0.8),legend.title = element_text(face="bold"))+
  labs(y="Number of mutations")

ggsave(p.PVV.LOH,filename=paste0(plots_dir,"Fig2e.pdf"),width=2.5,height=3)

#========================================#
# Now do the independent PVV simulation plotting ####
#========================================#

#Raw PVV simulation results - see the script raw_analysis_scripts/PVV_simulation_2.R for how these files are generated.
PVV_res=read.delim(paste0(root_dir,"/Data/simulation_results/PVV_sim_results.tsv"),stringsAsFactors=F)
Sample_IDs=unique(PVV_res$Sample_ID)
PVV_outcomes=read.delim(paste0(root_dir,"/Data/simulation_results/PVV_sim_outcomes.tsv"),stringsAsFactors=F)
PVV_outcomes$Sample_ID=rep(Sample_IDs,each=4)

PVV_outcomes<-PVV_outcomes%>%filter(Sample_ID%in%ref_table$Sample_ID)
PVV_outcomes$cat=sapply(PVV_outcomes$Sample_ID, function(sample) {ref_table$Category[ref_table$Sample_ID==sample]})
PVV_res<-PVV_res%>%filter(Sample_ID%in%ref_table$Sample_ID)

PVV_class_rename_vec=c("Assigned terminal branch","Not called PVV","Failed PVV assessment","Called as PVV")
names(PVV_class_rename_vec)=c("Assigned to terminal branch","Failed PVV assessment","FAIL","PASS")

## Generate Figure 2f ------
p.PVV.12.1<-PVV_outcomes%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Class=factor(PVV_class_rename_vec[Class],levels=PVV_class_rename_vec))%>%
  mutate(cat=factor(rename_cat_vec[cat],levels=rename_cat_vec))%>%
  mutate(data_type="Simulations")%>%
  ggplot(aes(x=Sample_ID,y=prop,fill=Class))+
  geom_bar(stat="identity")+
  facet_grid(cat~data_type,scales="free",space = "free")+
  guides(fill=guide_legend(nrow=1,title = ""))+
  coord_flip()+
  scale_fill_manual(values=c("#A2BBD4","#478EB0","#5AADC5","#D71920"))+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,0.5)) +
  labs(x="",y="Proportion")+
  theme_classic()+
  my_theme+
  theme(strip.text.y = element_text(size=8),strip.background = element_rect(fill=NA),legend.position = "bottom")
ggsave(p.PVV.12.1,filename = paste0(plots_dir,"Fig2f.pdf"),width=2,height=5.3)

## Generate Figure 2g ------
p.PVV.12.2<-PVV_res%>%
  filter(!is.na(Class))%>%
  mutate(Class=ifelse(Class=="PASS"&no_of_cell_divisions>6,"FAIL",Class))%>%
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,PASS=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+PASS),PASS_p=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL,-PASS)%>%
  gather(-Sample_ID,-data_type,key="Class",value="Proportion")%>%
  mutate(Class=gsub("_p","",Class))%>%
  mutate(Class=factor(PVV_class_rename_vec[Class],levels=PVV_class_rename_vec[3:4]))%>%
  left_join(summary_table_df%>%dplyr::select(Sample_ID,cat))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(data_type=rename_data_type_vec[data_type],cat=factor(rename_cat_vec[cat],levels=rename_cat_vec))%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,0.5)) +
  facet_grid(cat~data_type,scales="free",space = "free")+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 15),values=c("#DEEBf7","#C41E3A"))+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5))+labs(x="")+
  coord_flip()+
  theme(legend.position = "bottom",legend.margin=margin(0,0,0,0),legend.box.margin = margin(0,0,0,0,unit="mm"))


ggsave(p.PVV.12.2,filename = paste0(plots_dir,"Fig2g.pdf"),width=3.5,height=5.3)

#========================================#
#Calculate a WEIGHTED combined PVV proportions across samples ####
#========================================#
## 1. calculate the raw overall proportions in each sample - simulations vs data ----
PVV_prop_summary<-PVV_res%>%
  filter(data_set!="PR" & Sample_ID%in%ref_table$Sample_ID)%>%
  mutate(Class=ifelse(Class=="PASS"&no_of_cell_divisions>6,"FAIL",Class))%>%
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,PASS=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+PASS),PASS_p=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL,-PASS)

## 2. Calculate weights, based on the absolute number of MAVs found in each sample (regardless of class) - otherwise simulations from all samples will be weighted equally even though contribute very few MAVs to the actual data set ----
weights=PVV_res%>%filter(data_set!="PR"& Sample_ID%in%ref_table$Sample_ID)%>%filter(data_type=="data")%>%group_by(Sample_ID)%>%summarise(n=n())

## 3. Calculate the proportions of different classes in the simulations - including the weighting ----
PVV_sim_props<-PVV_prop_summary%>%
  filter(data_type=="simulation")%>%
  group_by(.)%>%
  summarise(data_type="simulation",FAIL=weighted.mean(FAIL_p,w=weights$n),PASS=weighted.mean(PASS_p,w=weights$n))%>%
  gather(-data_type,key="Class",value="Proportion")

## 4. Calculate the proportions of different classes in the data - including the weighting ----
PVV_dat_props<-PVV_res%>%
  filter(data_set!="PR" & Sample_ID%in%ref_table$Sample_ID)%>%
  group_by(Class,data_type)%>%
  mutate(Class=ifelse(Class=="PASS"&no_of_cell_divisions>6,"FAIL",Class))%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,PASS=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+PASS),PASS_p=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL,-PASS)%>%
  gather(-data_type,key="Class",value="Proportion")%>%
  mutate(Class=gsub("_p","",Class))%>%
  mutate(Class=factor(Class,levels=c("PASS","FAIL")))%>%
  filter(data_type=="data")

## Generate Figure 2h ------
rename_PVV_vec=c("Failed PVV assessment","Called as PVV")
names(rename_PVV_vec)=c("FAIL","PASS")
p.PVV.12.3<-rbind(PVV_dat_props,PVV_sim_props)%>%
  mutate(Class=factor(rename_PVV_vec[as.character(Class)],levels=rename_PVV_vec))%>%
  mutate(data_type=str_wrap(rename_data_type_vec[data_type],width=10))%>%
  ggplot(aes(x=data_type,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 15),values=c("#DEEBf7","#C41E3A"))+
  my_theme+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.text.x = element_text(angle=45,vjust = +0.3,hjust=0.5))+
  labs(x="",y="Proportions",fill="")
ggsave(p.PVV.12.3,filename = paste0(plots_dir,"Fig2h.pdf"),width=2.3,height=2.5)


#========================================#
# Now do the SOMATIC REVERSION SIMULATION plotting ####
#========================================#

sr_files=list.files(paste0(root_dir,"Data/simulation_results"),pattern="reversion",full.names=T)
list=lapply(sr_files,read.delim,stringsAsFactors=F)
dat=dplyr::bind_rows(list)
dat$mut_ref1=as.character(dat$mut_ref1)

nsim=2000
p.SR.1<-dat%>%
  group_by(Sample_ID,result)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "result",values_from="n")%>%
  mutate(FAIL=replace_na(FAIL,replace=0),PASS=replace_na(PASS,replace=0))%>%
  mutate(FAIL_p=FAIL/nsim,PASS_p=PASS/nsim)%>%
  mutate(Not_detected_p=1-FAIL_p-PASS_p)%>%
  dplyr::select(-FAIL,-PASS)%>%
  gather(key="Class",value="Proportion",-Sample_ID)%>%
  mutate(Class=gsub("_p","",Class))%>%
  mutate(Class=factor(Class,levels=c("Not_detected","FAIL","PASS")))%>%
  left_join(summary_table_df%>%dplyr::select(Sample_ID,cat))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1],cat=factor(rename_cat_vec[cat],levels=rename_cat_vec),data_type="Simulations")%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,0.5))+
  scale_fill_manual(values = c("#478EB0","#5AADC5","#D71920"))+
  facet_grid(cat~data_type,scales="free",space = "free")+
  my_theme+
  coord_flip()+
  theme(axis.title.y=element_blank(),legend.position="bottom")

ggsave(p.SR.1,filename = paste0(plots_dir,"Somatic_reversion_props.pdf"),width=3,height=2.5)
