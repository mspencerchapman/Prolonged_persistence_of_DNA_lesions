#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("ggplot2","dplyr","tidyr","gridExtra","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr")
bioconductor_packages=c("GenomicRanges","IRanges","MutationalPatterns","MASS","Rsamtools")

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
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

if(!require("phytools", character.only=T,quietly = T, warn.conflicts = F)){
  remotes::install_github("liamrevell/phytools")
  library(phytools)
}


#----------------------------------
# Set the ggplot2 theme and palettes for plotting
#----------------------------------

palette2=c("#E30613","#1D71B8")
palette6=c("#72e5ef", "#074d65", "#3d99ce", "#c257d3", "#d1add5", "#61356e")
palette8=c("#2EBAED" ,"#000000" ,"#DE1C14", "#E98C7B", "#D4D2D2" ,"#ADCC54" ,"#F0D0CE","blue")
my_theme=theme_classic(base_family="Helvetica")+theme(text=element_text(size=7,family="Helvetica"),
                                                      axis.text=element_text(size=5,family="Helvetica"),
                                                      strip.text = element_text(size=6,family="Helvetica"),
                                                      legend.key.height = unit(0.25, 'cm'),
                                                      legend.key.width = unit(0.25, 'cm'),
                                                      legend.title = element_text(size=6),
                                                      legend.text = element_text(size=5,family="Helvetica"))

#----------------------------------
# DO THE VARIOUS PLOTS RELATING TO MAV SIMULATIONS
#----------------------------------

functions_script_path="/lustre/scratch126/casm/team154pc/ms56/my_programs/Prolonged_persistence_functions.R"
root_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation"

source(functions_script_path)
setwd(root_dir)
ref_table=read.csv("Individual_ref.csv")
MAV_sim_dat<-data.table::fread("output2/MAV_sim_res/MAV_sim_results.tsv")%>%filter(Sample_ID%in%ref_table$Sample_ID)

#Do the first plot of class proportions across all the samples
MAV_props_df<-MAV_sim_dat%>%
  filter(data_set!="PR"&data_set!="SN"&!data_set=="MSC_chemo")%>%
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,removed=0,simple=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+removed+simple),removed_p=removed/(FAIL+removed+simple),simple_p=simple/(FAIL+removed+simple))%>%
  dplyr::select(-removed,-FAIL,-simple)%>%
  gather(-Sample_ID,-data_type,key="Class",value="Proportion")%>%
  mutate(Class=gsub("_p","",Class))

rename_MAV_vec=c("Simple","Unrelated","Separated")
names(rename_MAV_vec)=c("simple","FAIL","removed")
rename_data_type_vec=c("Observed data","Simulations")
names(rename_data_type_vec)=c("data","simulation")
MAVsims.p1<-MAV_props_df%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Class=factor(rename_MAV_vec[Class],levels=rename_MAV_vec))%>%
  mutate(data_type=rename_data_type_vec[data_type])%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  facet_wrap(~data_type,nrow=1)+
  scale_fill_manual(values = c("#66C2A5","lightgray","#FC8D62"))+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5)) 

ggsave(MAVsims.p1,filename = "plots/MAV_class_proportions_simulation.pdf",width=10,height=2)

#Now plot PASS MAV for the SN data_set
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
  mutate(Class=factor(Class,levels=c("PASS","FAIL")))%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  filter(Class=="PASS")%>%
  mutate(data_type=rename_data_type_vec[data_type])%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=data_type))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_brewer(palette = "Set2")+
  my_theme+
  #scale_y_log10()+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5))+
  labs(x="Sample ID",y="PASS MAV proportion",fill="")

ggsave(MAVsims.plot.SN,filename = "plots/MAV_SN_PASS_proportions_simulation.pdf",width=7,height=2)

#Now only include removed MAVS in the data that would pass the filtering i.e. would have no more than 1 negative clade.
#This would be only those with 2 cell divisions (≥3 would be filtered)
MAV_sim_dat%>%
  filter(data_set!="PR"&Class=="removed"&data_type=="simulation")%>%
  group_by(no_of_cell_divisions)%>%
  summarise(n=n())

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

MAVsims.p4<-MAV_props_df%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Class=factor(rename_MAV_vec[Class],levels=rename_MAV_vec))%>%
  mutate(data_type=rename_data_type_vec[data_type])%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  facet_wrap(~data_type,nrow=1)+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 10),values = c("#66C2A5","lightgray","#FC8D62"))+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5))+
  labs(x="")

ggsave(MAVsims.p4,filename = "plots/MAV_removed_nonfiltered_proportions_simulation.pdf",width=8,height=2)

#Do WEIGHTED combined MAV proportions across samples
#Plot overall proportions of classes across samples, data vs simulation
#To be an accurate representation of the data, simulation proportions should be weighted by the total number of MAVs recorded in the data in that sample

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

#Calculate weights, based on the absolute number of MAVs found in each sample (regardless of class) - otherwise simulations from all samples will be weighted equally even though contribute
#very few MAVs to the actual data set
MAV_weights=MAV_sim_dat%>%filter(data_set!="PR"&data_set!="SN"&!data_set=="MSC_chemo")%>%filter(data_type=="data")%>%group_by(Sample_ID)%>%summarise(n=n())

MAV_sim_props<-MAV_prop_summary%>%
  filter(data_type=="simulation")%>%
  group_by(.)%>%
  summarise(data_type="simulation",FAIL=weighted.mean(FAIL_p,w=MAV_weights$n),removed=weighted.mean(removed_p,w=MAV_weights$n),simple=weighted.mean(simple_p,w=MAV_weights$n))%>%
  gather(-data_type,key="Class",value="Proportion")

MAV_dat_props<-MAV_sim_dat%>%
  filter(data_set!="PR" & data_set!="SN")%>%
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

MAVsims.p5<-rbind(MAV_dat_props,MAV_sim_props)%>%
  mutate(Class=factor(rename_MAV_vec[Class],levels=rename_MAV_vec))%>%
  mutate(data_type=str_wrap(rename_data_type_vec[data_type],width=10))%>%
  ggplot(aes(x=data_type,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 10),values = c("#66C2A5","lightgray","#FC8D62"))+
  my_theme+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.text.x = element_text(angle=90),vjust = +0.5)+
  labs(x="",y="Proportions",fill="")
ggsave(MAVsims.p5,filename = "plots/MAV_combined_props_weighted.pdf",width=1.7,height=2)

#----------------------------------
# DO THE VARIOUS PLOTS RELATING TO PVV SIMULATIONS
#----------------------------------

res=read.delim("output2/PVV_sim_res/PVV_sim_results.tsv",stringsAsFactors=F)
Sample_IDs=unique(res$Sample_ID)
outcomes=read.delim("output2/PVV_sim_res/PVV_sim_outcomes.tsv",stringsAsFactors=F)
outcomes$Sample_ID=rep(Sample_IDs,each=4)

outcomes<-outcomes%>%filter(Sample_ID%in%ref_table$Sample_ID)
outcomes$cat=sapply(outcomes$Sample_ID, function(sample) {ref_table$Category[ref_table$Sample_ID==sample]})
res<-res%>%filter(Sample_ID%in%ref_table$Sample_ID)

PVV_class_rename_vec=c("Assigned terminal branch","Not called PVV","Failed PVV assessment","Called as PVV")
names(PVV_class_rename_vec)=c("Assigned to terminal branch","Failed PVV assessment","FAIL","PASS")

p.PVV.12.1<-outcomes%>%
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(Class=factor(PVV_class_rename_vec[Class],levels=PVV_class_rename_vec))%>%
  ggplot(aes(x=Sample_ID,y=prop,fill=Class))+
  geom_bar(stat="identity")+
  facet_grid(rows=vars(cat),scales="free",space = "free")+
  coord_flip()+
  scale_fill_manual(values=c("#A7C7E7","#2171B5","#87CEEB","#C41E3A"))+
  labs(x="",y="Proportion")+
  my_theme+
  theme(strip.text.x = element_text(size=5),legend.position = "top")
ggsave(p.PVV.12.1,filename = "plots/PVV_sim_outcomes.pdf",width=3,height=4)

p.PVV.12.2<-res%>%
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
  mutate(Sample_ID=stringr::str_split(Sample_ID,pattern="_",simplify=T)[,1])%>%
  mutate(data_type=rename_data_type_vec[data_type])%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  facet_wrap(~data_type,nrow=1)+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 15),values=c("#DEEBf7","#C41E3A"))+
  my_theme+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5))+labs(x="")

ggsave(p.PVV.12.2,filename = "plots/PVV_class_proportions_simulation.pdf",width=8,height=2)

#Calculate weights, based on the absolute number of PVVs found in each sample - otherwise simulations from all samples will be weighted equally even though contribute
#very few PVVs to the actual data set
PVV_prop_summary<-res%>%
  filter(data_set!="PR" & Sample_ID%in%ref_table$Sample_ID)%>%
  mutate(Class=ifelse(Class=="PASS"&no_of_cell_divisions>6,"FAIL",Class))%>%
  group_by(Sample_ID,Class,data_type)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "Class",values_from="n")%>%
  replace_na(list(FAIL=0,PASS=0))%>%
  mutate(FAIL_p=FAIL/(FAIL+PASS),PASS_p=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL,-PASS)

weights=res%>%filter(data_set!="PR"& Sample_ID%in%ref_table$Sample_ID)%>%filter(data_type=="data")%>%group_by(Sample_ID)%>%summarise(n=n())

sim_props<-PVV_prop_summary%>%
  filter(data_type=="simulation")%>%
  group_by(.)%>%
  summarise(data_type="simulation",FAIL=weighted.mean(FAIL_p,w=weights$n),PASS=weighted.mean(PASS_p,w=weights$n))%>%
  gather(-data_type,key="Class",value="Proportion")

dat_props<-res%>%
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

p.PVV.12.3<-rbind(dat_props,sim_props)%>%
  mutate(Class=factor(rename_vec[as.character(Class)],levels=rename_vec))%>%
  mutate(data_type=str_wrap(rename_data_type_vec[data_type],width=10))%>%
  ggplot(aes(x=data_type,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(labels=function(x) str_wrap(x, width = 15),values=c("#DEEBf7","#C41E3A"))+
  my_theme+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.text.x = element_text(angle=90,vjust = +0.5))+
  labs(x="",y="Proportions",fill="")
ggsave(p.PVV.12.3,filename = "plots/PVV_combined_props_weighted.pdf",width=1.7,height=2)


#----------------------------------
# DO THE VARIOUS PLOTS RELATING TO SOMATIC REVERSION SIMULATIONS
#----------------------------------

sr_files=list.files("simulation_results",pattern="reversion",full.names=T)
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
  select(-FAIL,-PASS)%>%
  gather(key="Class",value="Proportion",-Sample_ID)%>%
  mutate(Class=gsub("_p","",Class))%>%
  mutate(Class=factor(Class,levels=c("Not_detected","FAIL","PASS")))%>%
  ggplot(aes(x=Sample_ID,y=Proportion,fill=Class))+
  geom_bar(stat="identity",position="stack")+
  my_theme+
  coord_flip()

ggsave(p.SR.1,filename = "SR_test.pdf",width=3,height=2.5)

#Summary of median, mean and range of the PASS proportions
dat%>%
  group_by(Sample_ID,result)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "result",values_from="n")%>%
  mutate(FAIL=replace_na(FAIL,replace=0),PASS=replace_na(PASS,replace=0))%>%
  mutate(FAIL_p=FAIL/nsim,PASS_p=PASS/nsim)%>%
  mutate(Not_detected_p=1-FAIL_p-PASS_p)%>%
  select(-FAIL,-PASS)%>%
  gather(key="Class",value="Proportion",-Sample_ID)%>%
  mutate(Class=gsub("_p","",Class))%>%
  mutate(Class=factor(Class,levels=c("Not_detected","FAIL","PASS")))%>%
  filter(Class=="PASS")%>%
  summarise(mean=mean(Proportion),median=median(Proportion),min=min(Proportion),max=max(Proportion))

#
p.SR.2<-dat%>%
  ggplot(aes(x=no_of_cell_divisions))+
  geom_histogram()+
  facet_grid(rows=vars(Sample_ID))+
  my_theme


p.SR.2<-bind_rows(dat,PVV_mutations)%>%
  mutate(data_type=replace_na(data_type,replace="data"))%>%
  mutate(no_of_cell_divisions=factor(no_of_cell_divisions,levels=(2:12)))%>%
  group_by(Sample_ID,data_type,no_of_cell_divisions)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = "no_of_cell_divisions",values_from="n")%>%
  replace_na(.,replace=0)
ggplot(aes(x=no_of_cell_divisions))+
  geom_histogram()+
  facet_grid(Sample_ID~data_type)+
  my_theme
ggsave(p.SR.2,filename = "SR_test2.pdf",width=5,height=10)

