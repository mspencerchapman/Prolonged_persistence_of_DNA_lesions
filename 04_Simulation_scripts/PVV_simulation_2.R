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

#Total number of times to simulate a mutation occurring twice independently for each phylogeny
nsim=10000
re_run=F

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

Phasing_MAV_file_path="output2/Phasing_results_MAVs_all"
Phasing_PVV_file_path="output2/Phasing_results_PVVs_all"

#Set data file paths
data_sets=c("MSC_chemo","NW","MF","EM","KY","MSC_BMT","MSC_fetal")
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

#Add the phasing information
load(Phasing_PVV_file_path)
PVV_phasing_summary=sapply(Phasing_results_PVVs,extract_PVV_pos_clade_phasing_summary)
names(PVV_phasing_summary)<-str_split(names(PVV_phasing_summary),pattern="\\.",simplify=T)[,1]
table(PVV_phasing_summary)

mutations<-left_join(mutations,data.frame(Chrom_pos=names(PVV_phasing_summary),PVV_pos_clade_phasing=PVV_phasing_summary),by="Chrom_pos")

#Review the negative subclades - is there evidence of loss of the mutant allele? Can we prove it either way?
PVV_both_alleles_summary=sapply(Phasing_results_PVVs,extract_PVV_neg_clade_phasing_summary)
Neg_clades_phasing_df=data.frame(Chrom_pos=names(Phasing_results_PVVs),result=sapply(PVV_both_alleles_summary,function(x) return(x[1])))
#Simplify the output to only "No loss of mutant allele", "Unable to confirm", "Likely LOH" or "Positive clade have non-matching phasing"
Neg_clades_phasing_df$basic_result=sapply(Neg_clades_phasing_df$result,function(result) {
  if(result%in%c("Alt allele reads present","Both alleles confirmed with reference allele","Both alleles confirmed with reference allele in at least one subclade")){
    return("No loss of mutant allele")
  } else if(grepl("with maximum",result)){
    n_reads<-as.numeric(sapply(str_extract_all(result,"\\d"),paste,collapse=""))
    if(n_reads<=4) {return("Unable to confirm")} else if(n_reads>4) {return("Likely LOH")}
  } else if(result=="Positive clades have non-matching phasing"){return(NA)} else {return("Unable to confirm")}
})
#Join this to the mutations dataframe
mutations=left_join(mutations,Neg_clades_phasing_df[,c("Chrom_pos","basic_result")],by="Chrom_pos")
mutations$PVV_pos_clade_phasing<-sapply(mutations$PVV_pos_clade_phasing,function(x) {
  if(grepl("Non-matching phasing confirmed",x)) {
    res<-"Non-matching phasing confirmed"
  } else if(grepl("Same phasing confirmed",x)){
      res<-"Same phasing confirmed"
  } else if(is.na(x)) {
    res<-NA
  } else {
    res<-"Unable to confirm"
  }
  return(res)
})
mutations$PVV_pos_clade_phasing<-factor(mutations$PVV_pos_clade_phasing,levels=c("Same phasing confirmed","Unable to confirm","Non-matching phasing confirmed"))
myColors <- c("#33A02C","grey","#FB9A99")
names(myColors) <- levels(mutations$PVV_pos_clade_phasing)

##NOW START THE ANALYSIS SAMPLE BY SAMPLE
simulation_output_dir="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/simulation_results"
Sample_IDs=unique(mutations$Sample_ID)
all_PVV_res_list=lapply(data_sets,function(dataset) {
  data_set_PVVs<-mutations%>%filter(Type=="PVV" & data_set==dataset)
  Sample_IDs=unique(data_set_PVVs$Sample_ID)
  data_set_out=lapply(Sample_IDs,function(sample) {
    print(paste("Starting analysis for sample",sample))
    
    #Get the necessary files for the analysis
    sample_info=get_file_paths_and_project(dataset,Sample_ID=sample)
    tree=read.tree(sample_info$tree_file_path)
    load(sample_info$filtered_muts_path)
    
    #Get the data into the right format
    PVV_data=data_set_PVVs%>%dplyr::filter(Sample_ID==sample)
    PVV_data$data_type="data"
    PVV_data$Class<-ifelse(is.na(PVV_data$lesion_node),"FAIL","PASS")
    
    #Assuming each positive clade as independent mutation events
    PVV_data$by_chance_prob=sapply(1:nrow(PVV_data),function(i) {
      print(i)
      lesion_node=PVV_data$lesion_node[i]
      if(is.na(lesion_node)|!lesion_node%in%tree$edge[,1]){stop(return(NA))}
      mut=PVV_data$mut_ref1[i]
      pure_subclades=get_pure_subclades(mut1 = mut,lesion_node = lesion_node,matrices = list(NV=filtered_muts$COMB_mats.tree.build$NV,NR=filtered_muts$COMB_mats.tree.build$NR),tree = tree)
      if(any(pure_subclades=="Not PVV")){stop(return(NA))}
      pure_positive=pure_subclades[which(names(pure_subclades)=="pure_positive")]
      by_chance=prod(sapply(pure_positive,function(node) tree$edge.length[tree$edge[,2]==node]/sum(tree$edge.length)))
      score=-log10(by_chance)
      return(score)
    })
    PVV_data$by_chance_prob=as.numeric(PVV_data$by_chance_prob)
    
    #NOW START THE SIMULATIONS
    sim_res_file_path=paste0(simulation_output_dir,"/",sample,"_poor_fit_sim_res.RDS")
    
    if(file.exists(sim_res_file_path)&!re_run) {
      res<-readRDS(sim_res_file_path)
    } else {
      #Set up the simulation read count matrices
      sim=list()
      sim$NV=sim$NR=matrix(0,nrow=nsim,ncol=length(tree$tip.label))
      colnames(sim$NV)=colnames(sim$NR)=tree$tip.label
      sim$details=data.frame(mut_ref=1:nsim,node=NA,pval=NA)
      sim$sim_nodes<-data.frame(node1=numeric(nsim),node2=numeric(nsim))
      mean_depth=mean(rowMeans(filtered_muts$COMB_mats.tree.build$NR))
      
      #Simulate read counts
      for(i in 1:nsim) {
        sim$NR[i,]=rpois(length(tree$tip.label),lambda = mean_depth)
        nodes=base::sample(x = tree$edge[,2],size = 2,replace = F,prob = tree$edge.length)
        pos_samples=unlist(sapply(nodes,function(node) getTips(tree,node)))
        neg_samples=tree$tip.label[!tree$tip.label%in%pos_samples]
        sim$NV[i,pos_samples]=rbinom(n=length(pos_samples),size=sim$NR[i,pos_samples],prob = 0.5)
        sim$NV[i,neg_samples]=rbinom(n=length(neg_samples),size=sim$NR[i,neg_samples],prob = 1e-6)
        sim$sim_nodes[i,]<-nodes
      }
      
      #Run treemut - the majority are assigned to a private branch and will get binned in the next section (>80% depending on tree structure)
      include_ancestral_tip=ifelse(dataset=="MF",T,F)
      df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
      if(!include_ancestral_tip) {
        p.error = rep(0.01, ncol(sim$NV))
      } else if (include_ancestral_tip) { #If the tree includes an ancestral tip, add in a dummy ancestral sample to the NV, NR and p.error objects
        sim$NV<-cbind(sim$NV,"Ancestral"=rep(0,nrow(sim$NV)))
        sim$NR<-cbind(sim$NR,"Ancestral"=rep(10,nrow(sim$NR)))
        p.error = c(rep(0.01, ncol(sim$NV)),1e-6)
      }
      res = assign_to_tree(sim$NV[,df$samples],sim$NR[,df$samples], df, error_rate = p.error) #Get res (results!) object
      sim$details$node <- tree$edge[res$summary$edge_ml,2]
      sim$details$pval<-res$summary$pval
      
      #Now call PVVs as per the data
      mutations_to_test=which(!sim$details$node%in%1:length(tree$tip.label))
      length(mutations_to_test)
      
      print("Assessing mutations for whether they are truly phylogeny breaking")
      MC_CORES=1
      remove_duplicates=ifelse(dataset=="MSC_BMT",T,F)
      
      if(remove_duplicates) {
        pseudo_terminal_nodes=sapply(tree$edge[,2][!tree$edge[,2]%in%1:length(tree$tip.label)],function(node) {
          node_height=nodeheight(tree = tree,node=node)
          samples=get_edge_info(tree,details=NULL,node)$samples
          sample_heights=nodeHeights(tree)[tree$edge[,2]%in%which(tree$tip.label%in%samples),2]
          
          if(all((sample_heights-node_height)<30)){ #This is the method of determining the duplicates - may need to alter the "30" for low mutation burden samples
            return(node)
          }else{
            return(NA)
          }
        })
        pseudo_terminal_nodes=pseudo_terminal_nodes[!is.na(pseudo_terminal_nodes)]
        duplicate_samples=unlist(sapply(pseudo_terminal_nodes,get_all_node_children,tree))
      } else {
        duplicate_samples<-NULL
      }
      
      if(remove_duplicates) {
        mutations_to_test=sim$details$mut_ref[!sim$details$node%in%c(1:length(tree$tip.label),pseudo_terminal_nodes)]
      } else {
        mutations_to_test=sim$details$mut_ref[!sim$details$node%in%1:length(tree$tip.label)]
      }
      
      filter_output_df=create_PVV_filter_table(mutations_to_test,details=sim$details,tree,matrices=list(NV=sim$NV,NR=sim$NR),look_back="all",remove_duplicates = remove_duplicates,duplicate_samples=duplicate_samples,MC_CORES=MC_CORES)
      
      #Then apply cutoffs to these parameters to decide which of the poor fit mutations are "TRUE" phylogeny breakers (as per data analysis [take 2])
      filter_res=sapply(1:nrow(filter_output_df),
                        function(i) (filter_output_df[i,"pos_rho"]>0.1 & filter_output_df[i,"pos_test"])| #original cut-off is 0.15 for the pos_rho
                          (filter_output_df[i,"neg_rho"]>0.1 & filter_output_df[i,"neg_test"])) #had lower cut off of 0.1 if had strong pval result also
      
      print(paste0(sum(filter_res)," mutations pass the filtering to be phylogeny breaking mutations"))
      
      if(sum(filter_res)>0) {
        #Create a summary table of these mutations
        poor_fit_list=lapply(mutations_to_test[filter_res],function(mut) {
          df=data.frame(mut_ref=mut,
                        node=sim$details$node[mut],
                        Sample_ID=sample,
                        lesion_node=NA,
                        lesion_timing=NA,
                        lesion_repair_node=NA,
                        lesion_repair_timing=NA,
                        lesion_duration=NA,
                        no_of_cell_divisions=NA
          )
          return(df)
        })
        poor_fit_df=dplyr::bind_rows(poor_fit_list)
        
        #Work out the"initial lesion" node, and the "latest normal clade" nodes
        print(paste0("Starting assessment of the 'lesion node' and 'lesion repair node' for the PVVs"))
        for(i in 1:nrow(poor_fit_df)) {
          print(i)
          allocated_node=poor_fit_df$node[i]
          mut=poor_fit_df$mut_ref[i]
          print(mut)
          
          initial_lesion_node=find_PVV_lesion_node(mut=mut,allocated_node=poor_fit_df$node[i],pos_test = filter_output_df$pos_test[filter_output_df$mut==mut],neg_test=filter_output_df$neg_test[filter_output_df$mut==mut],tree,matrices=list(NV=sim$NV,NR=sim$NR))
          initial_lesion_timing=nodeheight(tree,initial_lesion_node)
          pure_subclades=get_pure_subclades(mut1 = mut, lesion_node = initial_lesion_node,tree=tree,matrices=list(NV=sim$NV,NR=sim$NR))
          if(any(pure_subclades=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")|
             any(pure_subclades=="Not PVV")){next}
          #Now get the earliest time that the lesion may have been repaired
          #1. get daughter nodes of lesion node
          
          #Update the df
          poor_fit_df$lesion_node[i]<-initial_lesion_node
          poor_fit_df$lesion_timing[i]<-nodeheight(tree,initial_lesion_node)
          poor_fit_df$lesion_repair_node[i]<-get_ancestor_node(node = tail(pure_subclades,n=1),tree)
          poor_fit_df$lesion_repair_timing[i]<-nodeheight(tree,poor_fit_df$lesion_repair_node[i])
          poor_fit_df$lesion_duration[i]<-poor_fit_df$lesion_repair_timing[i]-poor_fit_df$lesion_timing[i]
          poor_fit_df$no_of_cell_divisions[i]<-length(pure_subclades)
        }
      }
      
      poor_fit_df$Class=ifelse(is.na(poor_fit_df$lesion_node),"FAIL","PASS")
      table(poor_fit_df$Class)
      
      #Assuming each positive clade as independent mutation events
      poor_fit_df$data_type="simulation"
      poor_fit_df$data_set=dataset
      
      poor_fit_df$by_chance_prob=sapply(1:nrow(poor_fit_df),function(i) {
        lesion_node=poor_fit_df$lesion_node[i]
        if(is.na(lesion_node)){stop(return(NA))}
        mut=poor_fit_df$mut_ref[i]
        pure_subclades=get_pure_subclades(mut1 = mut,lesion_node = lesion_node,matrices = list(NV=sim$NV,NR=sim$NR),tree = tree)
        pure_positive=pure_subclades[which(names(pure_subclades)=="pure_positive")]
        by_chance=prod(sapply(pure_positive,function(node) tree$edge.length[tree$edge[,2]==node]/sum(tree$edge.length)))
        score=-log10(by_chance)
        return(score)
      })
      
      n_failed_PVV_assessment=length(mutations_to_test)-sum(filter_res)
      
      poor_fit_df$PVV_pos_clade_phasing="Unable to confirm"
      sim_outcome_df<-poor_fit_df%>%
        group_by(Class)%>%
        summarise(n=n())%>%
        bind_rows(data.frame(Class=c("Assigned to terminal branch","Failed PVV assessment"),n=c(nsim-length(mutations_to_test),n_failed_PVV_assessment)))%>%
        mutate(prop=n/sum(n),Class=factor(Class,levels=c("Assigned to terminal branch","Failed PVV assessment","FAIL","PASS")))
      
      res=list(sim_outcomes=sim_outcome_df,PVV_data_and_sim=dplyr::bind_rows(poor_fit_df,PVV_data))
      saveRDS(object=res,file = sim_res_file_path)
    }

    

    # 
    # p1<-poor_fit_df%>%
    #   group_by(Class)%>%
    #   summarise(n=n())%>%
    #   bind_rows(data.frame(Class=c("Assigned to terminal branch","Failed PVV assessment"),n=c(nsim-length(mutations_to_test),n_failed_PVV_assessment)))%>%
    #   mutate(prop=n/sum(n),Class=factor(Class,levels=c("Assigned to terminal branch","Failed PVV assessment","FAIL","PASS")))%>%
    #   ggplot(aes(x=Class,y=n,fill=Class))+
    #   geom_bar(stat="identity")+
    #   scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 12, simplify = FALSE), paste, collapse="\n"))+
    #   geom_text(aes(label=n,y=n+100),size=3)+
    #   my_theme+
    #   labs(title="Outcome of simulated independent\n mutations at same site")
    # 
    # p2<-bind_rows(poor_fit_df,PVV_data)%>%
    #   filter(Class=="PASS")%>%
    #   ggplot(aes(x=lesion_duration,y=by_chance_prob)) +
    #   geom_point()+
    #   facet_grid(cols=vars(data_type))+
    #   my_theme
    # 
    # p3<-poor_fit_df%>%
    #   filter(Class=="PASS")%>%
    #   ggplot(aes(x=lesion_duration,y=by_chance_prob)) +
    #   stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white",size=0.1)+
    #   geom_point(data=PVV_data,aes(col=PVV_pos_clade_phasing),size=3,alpha=0.75)+
    #   scale_color_manual(name = "PVV_pos_clade_phasing",values = myColors)+
    #   scale_fill_continuous(guide='none')+
    #   my_theme
    # 
    # p4<-poor_fit_df%>%
    #   filter(Class=="PASS")%>%#
    #   ggplot(aes(x=lesion_duration,y=no_of_cell_divisions)) +
    #   stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white",size=0.1)+
    #   my_theme+
    #   geom_jitter(data=PVV_data,aes(col=PVV_pos_clade_phasing),size=3)+
    #   scale_color_manual(name = "PVV_pos_clade_phasing",values = myColors)+
    #   scale_fill_continuous(guide='none')
    # 
    # PVV_data_probs<-PVV_data%>%mutate(Class=factor(Class,levels=c("FAIL","PASS")))%>%count(Class)%>%complete(Class,fill=list(n=0))%>%mutate(prop=n/sum(n),data_type="data")
    # PVV_sim_probs<-poor_fit_df%>%filter(Class!="Not PVV")%>%mutate(Class=factor(Class,levels=c("FAIL","PASS")))%>%count(Class)%>%complete(Class,fill=list(n=0))%>%mutate(prop=n/sum(n),data_type="simulation")
    # 
    # p5<-bind_rows(PVV_data_probs,PVV_sim_probs)%>%
    #   ggplot(aes(x=data_type,y=prop,fill=Class))+
    #   geom_bar(stat="identity",position="stack",col="black",size=0.15)+
    #   my_theme+
    #   theme(text=element_text(size=8))+
    #   labs(fill="",
    #        x="",
    #        y="Proportion of total PVVs",
    #        title=paste0(sample,": Data v simulation;\n Proportion of PVVs having orientation\n expected for persistent DNA lesion"))
    # 
    # plot_comb1<-arrangeGrob(p1,p5,ncol=2)
    # plot_comb<-arrangeGrob(plot_comb1,p3,nrow=2)
    # 
    #ggsave(filename=paste0("output2/PVV_vs_sim_plots/",sample,"_PVV_vs_sim_plots.pdf"),plot=plot_comb,device="pdf",width=8,height=6)
    
    return(res)
  })
  return(data_set_out)
})

#Unlist the first layer (data set layer) of the list
all_PVV_res_list=unlist(all_PVV_res_list,recursive = F)

#Extract the simulation outcomes and save as .tsv
sim_outcomes_list=lapply(all_PVV_res_list,function(list) return(list[[1]])) 
sim_outcomes_df=dplyr::bind_rows(sim_outcomes_list)
write.table(sim_outcomes_df,file="output2/PVV_sim_res/PVV_sim_outcomes.tsv",quote=F,sep="\t",row.names=F)

#Extract the actual simulation/ data PVVs including nodes/ FAIL-PASS outcomes and save as .tsv
PVV_res_list=lapply(all_PVV_res_list,function(list) return(list[[2]]))
PVV_res_df=dplyr::bind_rows(PVV_res_list)
write.table(PVV_res_df,file="output2/PVV_sim_res/PVV_sim_results.tsv",quote=F,sep="\t",row.names=F)


