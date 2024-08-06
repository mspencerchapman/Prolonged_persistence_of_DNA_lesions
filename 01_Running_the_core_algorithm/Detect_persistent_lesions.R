#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("stringr","ape","seqinr","tidyr","dplyr","ggplot2","phangorn","optparse","parallel","MASS","VGAM","remotes","dichromat")
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
# Parse the options list (using optparse) ####
#========================================#

option_list = list(
  make_option(c("-w", "--working_dir"), action="store", default=NULL, type='character', help="Working directory for analysis - should be set to the cloned github directory. If not set will default to the current directory"),
  make_option(c("-s", "--sample_id"), action="store", default='unspecified', type='character', help="Sample ID to include in output file names"),
  make_option(c("-t", "--tree_path"), action="store", default=NULL, type='character', help="path for the tree file - mandatory"),
  make_option(c("-f", "--filtering_muts_path"), action="store", type='character', help="path for the filtered_muts file - mandatory, and needs to be in a very specific structure"),
  make_option(c("-g", "--genome_file"), action="store", default='/nfs/cancer_ref02/human/GRCh37d5/genome.fa',type='character', help="path to the relevant genome file"),
  make_option(c("-a", "--ancestral"), action="store_true", default=FALSE, type='logical', help="Set this flag if input trees contain an ancestral tip"),
  make_option(c("-d", "--duplicate_remove"), action="store_true", default=FALSE, type='logical', help="Remove any duplicate samples that are present"),
  make_option(c("-r", "--re_run"), action="store_true", default=FALSE, type='logical', help="Force rerun of TREEMUT and the PVV assessment, even if conditions are met for re-importing existing intermediate files"),
  make_option(c("-c", "--cores"), action="store_true", default=1, type='numeric', help="Number of cores to use"),
  make_option(c("-o", "--output_dir"), action="store", default="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/output", type='character', help="output directory for files")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

print(opt)

#========================================#
# Set various variables from the options list ####
#========================================#

#Set file paths and script options from command args
if(!is.null(opt$w)) {root_dir = opt$w} else {root_dir<-getwd()}
sampleID=opt$s
tree_file_path = opt$t
filtered_muts_file = opt$f
genome_file=opt$g
output_dir=opt$o
include_ancestral_tip=opt$a  #Does your tree include an ancestral tip that you want to keep in the analysis?
remove_duplicates=opt$d
MC_CORES=opt$c
re_run=opt$r
print(paste(MC_CORES,"cores available"))

functions_file=source(paste0(root_dir,"Data/Prolonged_persistence_functions.R"))
filter_output_df_file=paste0(output_dir,"/",sampleID,"_filter_df.tsv")

#Load up the objects, and pull out the details data frame and the NV and NR matrices
print("Loading the tree and converting to multifurcating structure")
tree <- di2multi(read.tree(tree_file_path))
print("Loading the mutations matrix")
load(filtered_muts_file)

details=filtered_muts$COMB_mats.tree.build$mat
NV = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
NR = as.matrix(filtered_muts$COMB_mats.tree.build$NR)

#Re-derive the res object & put the pval and node into the table
if(re_run|!all(c("node","pval")%in%colnames(details))|sum(tree$edge.length)!=nrow(details)|any(!details$node%in%tree$edge[,2])) {
  print("Reassigning mutations to branches using the treemut package")
  if(!include_ancestral_tip) {
    p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
  } else if (include_ancestral_tip) { #If the tree includes an ancestral tip, add in a dummy ancestral sample to the NV, NR and p.error objects
    NV<-cbind(NV,"Ancestral"=rep(0,nrow(NV)))
    NR<-cbind(NR,"Ancestral"=rep(10,nrow(NR)))
    p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)),1e-6)
  }
  res = assign_to_tree(tree=tree,mtr=NV[,tree$tip.label], depth=NR[,tree$tip.label], error_rate = p.error) #Get res (results!) object
  details$node <- tree$edge[res$summary$edge_ml,2]
  details$pval <- res$summary$pval
  
  #Update the original objects and re-save to the original file-paths so that can avoid re-running treemut in future
  tree$edge.length<-res$df$df$edge_length
  filtered_muts$COMB_mats.tree.build$mat<-details
  filtered_muts$COMB_mats.tree.build$NV<-NV
  filtered_muts$COMB_mats.tree.build$NR<-NR
  save(filtered_muts,file=filtered_muts_file)
  write.tree(tree,file = tree_file_path)
}

#CORRECTLY CLASSIFY MNVs
filtered_muts$COMB_mats.tree.build=reclassify_MNVs(COMB_mats=filtered_muts$COMB_mats.tree.build,region_size=2,genomeFile = genome_file)

#Add a "Chrom_pos" column for identifying multi-allelic variants
details=filtered_muts$COMB_mats.tree.build$mat
details$Chrom=as.character(details$Chrom)
details$Ref=as.character(details$Ref)
details$Alt=as.character(details$Alt)
NV = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
NR = as.matrix(filtered_muts$COMB_mats.tree.build$NR)

details$Chrom_pos=paste(details$Chrom,details$Pos,sep = "-") #Create a "chrom_pos" column

#If have multiple samples from the same colony the actual "private branches" become apparently internal.
#Need to handle these as they mess up the filtering
if(remove_duplicates) {
  pseudo_terminal_nodes=sapply(tree$edge[,2][!tree$edge[,2]%in%1:length(tree$tip.label)],function(node) {
    node_height=nodeheight(tree = tree,node=node)
    samples=get_edge_info(tree,details,node)$samples
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

#Now pull out the mutations of interest and plot
treefit_pval_cutoff = 0.5 #Decide how bad the binomial fit should be to review the mutation
poor_fit = details$pval < treefit_pval_cutoff  #See how many mutations have read counts that don't fit the tree very well
sum(poor_fit) #How many poor fit mutations are there?  Perhaps there is a problem with tree building?
#Get the mut refs of those mutations that are "poor fit" and not assigned to terminal branches (these are almost always artefact)
if(remove_duplicates) {
  mutations_to_test=details$mut_ref[!details$node%in%c(1:length(tree$tip.label),pseudo_terminal_nodes)]
  #poor_fit_mutations = details$mut_ref[poor_fit & !details$node%in%c(1:length(tree$tip.label),pseudo_terminal_nodes)]
} else {
  mutations_to_test=details$mut_ref[!details$node%in%1:length(tree$tip.label)]
  #poor_fit_mutations = details$mut_ref[poor_fit & !details$node%in%1:length(tree$tip.label)]
}

print(paste0("There are ",length(mutations_to_test)," mutations on internal branches to test"))
#print(paste0("There are ",length(poor_fit_mutations)," mutations with a pval <",treefit_pval_cutoff," on internal branches"))

#Apply a filter to try and pull out those mutations that are truly "phylogeny breaking"
#Notes: a "phylogeny breaking" mutation will have either:
#(1) a negative subclade within the allocated clade
#(2) a positive clade outside (though close to) the allocated clade
#This represents the same process but just depends on which node the "treemut" package has decided to place the mutation when presented with one that breaks the phylogeny.
#In the case of (1) there will be overdispersion of the read counts within the allocated clade
#In the case of (2) there will be overdispersion of the read counts outside the allocated clade
#This filter assesses the overdispersion within ("pos_rho") and outside ("neg_rho) the allocated clade.  It does not just assess individual samples, but all the clades (i.e. phylogenetically
#grouped samples) that are anticipated to be entirely positive or negative, by summing the read counts in each clade.
if(length(mutations_to_test)>0) {
  if(file.exists(filter_output_df_file)&!re_run) {
    print('Existing filter output file found in the output directory and will be imported')
    filter_output_df=read.delim(filter_output_df_file,header=T)
  } else {
    filter_output_df=create_PVV_filter_table(mutations_to_test,details,tree,matrices=list(NV=NV,NR=NR),look_back="all",remove_duplicates = remove_duplicates,duplicate_samples=duplicate_samples,MC_CORES=MC_CORES)
    write.table(filter_output_df,file=filter_output_df_file,quote=F,sep="\t",row.names = F)
  }
  
  #Then apply cutoffs to these parameters to decide which of the poor fit mutations are "TRUE" phylogeny breakers
  filter_res=sapply(1:nrow(filter_output_df),
                    function(i) (filter_output_df[i,"pos_rho"]>0.1 & filter_output_df[i,"pos_test"])|
                      (filter_output_df[i,"neg_rho"]>0.1 & filter_output_df[i,"neg_test"]))
  
  print(paste0(sum(filter_res)," mutations pass the filtering to be phylogeny breaking mutations"))
} else {
  filter_res<-F
}

#Create reference set of sample sets that form clades - used in assessing PVVs and MAVs
all_clades=unique(tree$edge[,2])
all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))

if(sum(filter_res)>0) {
  #Create a summary table of these mutations
  poor_fit_list=lapply(mutations_to_test[filter_res],function(mut_ref) {
    node=details$node[details$mut_ref==mut_ref]
    alts=details$Alt[details$mut_ref==mut_ref]
    df=data.frame(Chrom_pos=details$Chrom_pos[details$mut_ref==mut_ref],
                  mut_ref1=mut_ref,
                  mut_ref2=NA,
                  Ref=details$Ref[details$mut_ref==mut_ref][1],
                  Alt1=alts[1],
                  Alt2=details$Ref[details$mut_ref==mut_ref],
                  Mut_type1=details$Mut_type[details$mut_ref==mut_ref],
                  Mut_type2="No_mut",
                  Node1=node,
                  Node2=NA,
                  Type="PVV",
                  Class=NA,
                  Sample_ID=sampleID,
                  lesion_node=NA,
                  lesion_timing=NA,
                  lesion_repair_node=NA,
                  lesion_repair_timing=NA,
                  lesion_duration=NA,
                  no_of_cell_divisions=NA,
                  base_order=NA,
                  MAV_is_PVV=NA,
                  colonies_per_negative_subclade=NA,
                  depth_per_negative_subclade=NA
    )
    return(df)
  })
  poor_fit_df=dplyr::bind_rows(poor_fit_list)
  
  #Work out the "initial lesion" node, and the "latest normal clade" nodes
  print(paste0("Starting assessment of the 'lesion node' and 'lesion repair node' for the PVVs"))
  for(i in 1:nrow(poor_fit_df)) {
    print(i)
    allocated_node=poor_fit_df$Node1[i]
    mut=poor_fit_df$mut_ref1[i]
    print(mut)
    
    initial_lesion_node=find_PVV_lesion_node(mut=mut,allocated_node=poor_fit_df$Node1[i],pos_test = filter_output_df$pos_test[filter_output_df$mut==mut],filter_output_df$neg_test[filter_output_df$mut==mut],tree,matrices=list(NV=NV,NR=NR))
    initial_lesion_timing=nodeheight(tree,initial_lesion_node)
    pure_subclades=get_pure_subclades(mut1 = mut, lesion_node = initial_lesion_node,tree=tree,matrices=list(NV=NV,NR=NR))
    if(any(pure_subclades=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")|
       any(pure_subclades=="Not PVV")){next}
    #Now get the earliest time that the lesion may have been repaired
    #1. get daughter nodes of lesion node
    if(poor_fit_df$Mut_type1[i]!="SNV") {
      bases=c("Ref","Mut"); names(bases)=c("pure_negative","pure_positive")
    } else {
      bases=str_split(mut,pattern="-",simplify = TRUE)[,3:4];names(bases)=c("pure_negative","pure_positive")
    }
    base_order=bases[names(pure_subclades)] #PVVs always start with an "Alt" base
    
    #Determine how strong evidence is for each negative subclade by recording number of colonies & total depth in each
    n_negative_colonies=sapply(pure_subclades[which(names(pure_subclades)=="pure_negative")],function(node) length(getTips(tree,node))) #How many colonies per negative subclade
    depth_in_negative=sapply(pure_subclades[which(names(pure_subclades)=="pure_negative")],function(node) {sum(NR[mut,getTips(tree,node)])}) #Total depth at mutation site for each
    
    #Update the df
    poor_fit_df$lesion_node[i]<-initial_lesion_node
    poor_fit_df$lesion_timing[i]<-nodeheight(tree,initial_lesion_node)
    poor_fit_df$lesion_repair_node[i]<-get_ancestor_node(node = tail(pure_subclades,n=1),tree)
    poor_fit_df$lesion_repair_timing[i]<-nodeheight(tree,poor_fit_df$lesion_repair_node[i])
    poor_fit_df$lesion_duration[i]<-poor_fit_df$lesion_repair_timing[i]-poor_fit_df$lesion_timing[i]
    poor_fit_df$no_of_cell_divisions[i]<-length(pure_subclades)
    poor_fit_df$base_order[i]<-paste(base_order,collapse = "-")
    poor_fit_df$depth_per_negative_subclade[i]<-paste(depth_in_negative,collapse=",")
    poor_fit_df$colonies_per_negative_subclade[i]<-paste(n_negative_colonies,collapse=",")
  }
}

#Now look for the different mutations at the same site
print("Assessing mutations for duplicates at the same site")
MAV_list=get_multi_allelic_variant_list(details)
print(paste0("There are ",length(MAV_list)," mutations at overlapping sites"))

if(any(sapply(MAV_list,length)>2)) {
  print("Evidence of triallelic variants that cannot be handled automatically. These are being removed and should be analysed manually.")
  print(MAV_list[sapply(MAV_list,length)>2])
  MAV_list[sapply(MAV_list,length)>2]<-NULL
}

if(length(MAV_list)!=0) {
  multi_allelic_list=lapply(MAV_list,function(muts) {
    Ref1=details$Ref[details$mut_ref==muts[1]]; Ref2=details$Ref[details$mut_ref==muts[2]]
    Alt1<-details$Alt[details$mut_ref==muts[1]]; Alt2<-details$Alt[details$mut_ref==muts[2]]
    Pos1<-as.numeric(details$Pos[details$mut_ref==muts[1]]); Pos2<-as.numeric(details$Pos[details$mut_ref==muts[2]])
    
    new_refs=establish_ref_and_alt(Ref1,Ref2,Alt1,Alt2,Pos1,Pos2)
    
    df=data.frame(Chrom_pos=details$Chrom_pos[details$mut_ref==muts[1]],
                  mut_ref1=muts[1],
                  mut_ref2=muts[2],
                  Ref=new_refs$Ref,
                  Alt1=new_refs$Alt1,
                  Alt2=new_refs$Alt2,
                  Mut_type1=details$Mut_type[details$mut_ref==muts[1]],
                  Mut_type2=details$Mut_type[details$mut_ref==muts[2]],
                  Node1=details$node[details$mut_ref==muts[1]],
                  Node2=details$node[details$mut_ref==muts[2]],
                  Type="MAV",
                  Class=NA,
                  Sample_ID=sampleID,
                  lesion_node=NA,
                  lesion_timing=NA,
                  lesion_repair_node=NA,
                  lesion_repair_timing=NA,
                  lesion_duration=NA,
                  no_of_cell_divisions=NA,
                  base_order=NA,
                  MAV_is_PVV=NA,
                  colonies_per_negative_subclade=NA,
                  depth_per_negative_subclade=NA
    )
    return(df)
  })
  multi_allelic_df=dplyr::bind_rows(multi_allelic_list)
  
  #Apply a filtering step to only include multi-allelic mutations that are adjacent to each other, or one is "within" the other
  print("Filtering MAVs and allocating lesion node/ lesion repair nodes")
  for(i in 1:nrow(multi_allelic_df)){
    print(paste(i,multi_allelic_df$Chrom_pos[i]))
    
    lesion_node_res=find_MAV_lesion_node(node1=multi_allelic_df$Node1[i],node2=multi_allelic_df$Node2[i],tree,Chrom=str_split(multi_allelic_df$Chrom_pos[i],pattern = "-",simplify=T)[,1])
    
    #multi_allelic_df$Filter[i]<-lesion_node_res$Filter
    multi_allelic_df$Class[i]<-lesion_node_res$Class
    multi_allelic_df$lesion_node[i]<-lesion_node_res$initial_lesion_node
    
    #select the duplicates that may represent lesion persistence events & the relationship between the branches
    if(lesion_node_res$Class=="simple"){
      multi_allelic_df$lesion_timing[i]<-nodeheight(tree,lesion_node_res$initial_lesion_node)
      multi_allelic_df$lesion_repair_node[i]<-lesion_node_res$initial_lesion_node
      multi_allelic_df$no_of_cell_divisions[i]=2
      multi_allelic_df$lesion_repair_timing[i]<-nodeheight(tree,lesion_node_res$initial_lesion_node)
      multi_allelic_df$lesion_duration[i]=0
      multi_allelic_df$base_order[i]=paste(c(multi_allelic_df$Alt1[i],multi_allelic_df$Alt2[i]),collapse="-")
      multi_allelic_df$MAV_is_PVV[i]="No"
      
    }else if(lesion_node_res$Class=="removed") {
      
      pure_subclades=get_pure_subclades(mut1 = multi_allelic_df$mut_ref1[i], mut2 = multi_allelic_df$mut_ref2[i],lesion_node = lesion_node_res$initial_lesion_node,tree=tree,matrices=list(NV=NV,NR=NR))
      if(any(pure_subclades=="More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion")|
         any(pure_subclades=="Not MAV")){
        multi_allelic_df$Class[i]<-"FAIL"
        next
        }
      
      #3. Define the "bases" vector for this mutation - used as a reference for the base order
      if(all(multi_allelic_df[i,c("Mut_type1","Mut_type2")]=="SNV")) {
        bases=c(multi_allelic_df$Ref[i],multi_allelic_df$Alt1[i],multi_allelic_df$Alt2[i])
        names(bases)=c("pure_negative","pure_mut1","pure_mut2")
      } else {
        bases=c("Ref","Mut1","Mut2") #Use a generic code if they are not both SNVs
        names(bases)=c("pure_negative","pure_mut1","pure_mut2")
      }
      base_order=bases[names(pure_subclades)]
      
      #Determine how strong evidence there is for each pure subclade by recording number of colonies & total depth in each
      #This is particularly important for the negative subclades
      if(any(names(pure_subclades)=="pure_negative")) {
        colonies_per_fixed_clade=sapply(pure_subclades[which(names(pure_subclades)=="pure_negative")],function(node) length(getTips(tree,node))) #How many colonies per negative subclade
        depth_per_fixed_clade=sapply(pure_subclades[which(names(pure_subclades)=="pure_negative")],function(node) {sum(NR[multi_allelic_df$mut_ref1[i],getTips(tree,node)])}) #Total depth at mutation site for each
      } else {
        colonies_per_fixed_clade<-depth_per_fixed_clade<-NA
      }
      
      #Fill out the df with this acquired info
      multi_allelic_df$lesion_repair_node[i]=get_ancestor_node(tail(pure_subclades,n=1),tree)
      multi_allelic_df$no_of_cell_divisions[i]=length(pure_subclades)
      multi_allelic_df$lesion_timing[i]<-nodeheight(tree,lesion_node_res$initial_lesion_node)
      multi_allelic_df$lesion_repair_timing[i]<-nodeheight(tree,multi_allelic_df$lesion_repair_node[i])
      multi_allelic_df$lesion_duration[i]<-nodeheight(tree,multi_allelic_df$lesion_repair_node[i])-nodeheight(tree,multi_allelic_df$lesion_node[i])
      multi_allelic_df$base_order[i]=paste(base_order,collapse="-")
      multi_allelic_df$MAV_is_PVV[i]=ifelse(sum(base_order==bases["pure_mut1"])>1|sum(base_order==bases["pure_mut2"])>1,"Yes","No")
      multi_allelic_df$colonies_per_negative_subclade[i]=paste(colonies_per_fixed_clade,collapse=",")
      multi_allelic_df$depth_per_negative_subclade[i]=paste(depth_per_fixed_clade,collapse=",")
    } else {
      #multi_allelic_df$Filter[i]<-lesion_node_res$Filter
    }
  }
    
  #Keep only those that pass the above filter, and remove the "Filter" variable
  print(paste0("There are ",sum(multi_allelic_df$Class!="FAIL")," mutations passing filtering to be true multi-allelic variant"))
  if(exists("poor_fit_df")) {
    poor_fit_df <- poor_fit_df[!poor_fit_df$Chrom_pos%in%multi_allelic_df$Chrom_pos,]
    combined_df=rbind(poor_fit_df,multi_allelic_df)
  } else {
     combined_df<-multi_allelic_df
  }
} else if(exists("poor_fit_df")){
   combined_df<-poor_fit_df
}

#SAVE PLOTS OF THE MUTATIONS
if(exists("poor_fit_df")){
  #Plot the poor fit mutations
  print("Saving plots of the PVVs")
  pdf(paste0(output_dir,"/",sampleID,"_PVVs_plots.pdf"),width=20,height=6)
  for(i in 1:nrow(poor_fit_df)) {
    tree=plot_tree(tree=tree,cex.label=0)
    add_annotation(tree,
                   details=details,
                   matrices=list(NV=NV, NR=NR),
                   annot_function = plot_MAV_mut,
                   mut1=paste(poor_fit_df$Chrom_pos[i],poor_fit_df$Ref[i],poor_fit_df$Alt1[i],sep="-"),
                   lesion_node=poor_fit_df$lesion_node[i],
                   cex=0.4)
    if(!is.na(poor_fit_df$lesion_node[i]) & poor_fit_df$lesion_node[i]!=(Ntip(tree) + 1)) {
      lesion_node_info=get_edge_info(tree,details,poor_fit_df$lesion_node[i])
      points(x=lesion_node_info$x,y=lesion_node_info$yb,type ="p" ,pch=19,bg="black",col="black")
      lesion_repair_node_info=get_edge_info(tree,details,poor_fit_df$lesion_repair_node[i])
      points(x=lesion_repair_node_info$x,y=lesion_repair_node_info$yb,type ="p" ,pch=15,bg="black",col="black")
    }
  }
  dev.off()
}

if(exists("multi_allelic_df")){
  #Plot the "PASS" MAVs on the phylogeny for manual review
  multi_allelic_df_pass<-multi_allelic_df[multi_allelic_df$Class!="FAIL",]
  if(nrow(multi_allelic_df_pass)>0){
    print("Saving plots of the multi-allelic variants")
    pdf(paste0(output_dir,"/",sampleID,"_MAVs_single_plots.pdf"),width=20,height=6)
    for(i in 1:nrow(multi_allelic_df_pass)) {
      if(multi_allelic_df_pass$Mut_type1[i]=="SNV"){
        mut1=multi_allelic_df_pass$mut_ref1[i]
      } else {
        mut1=multi_allelic_df_pass$mut_ref1[i]
      }
      if(multi_allelic_df_pass$Mut_type2[i]=="SNV"){
        mut2=multi_allelic_df_pass$mut_ref2[i]
      } else {
        mut2=multi_allelic_df_pass$mut_ref2[i]
      }
      
      tree=plot_tree(tree=tree,cex.label=0)
      add_annotation(tree,
                     details=details,
                     matrices=list(NV=NV,NR=NR),
                     annot_function = plot_MAV_mut,
                     mut1=mut1,
                     mut2=mut2,
                     lwd=2,
                     cex=0.4)
      if(!is.na(multi_allelic_df_pass$lesion_node[i]) & multi_allelic_df_pass$lesion_node[i]!=(Ntip(tree) + 1)) {
        lesion_node_info=get_edge_info(tree,details,multi_allelic_df_pass$lesion_node[i])
        lesion_repair_node_info=get_edge_info(tree,details,multi_allelic_df_pass$lesion_repair_node[i])
        points(x=lesion_node_info$x,y=lesion_node_info$yb,type ="p" ,pch=19,bg="black",col="black")
        points(x=lesion_repair_node_info$x,y=lesion_repair_node_info$yb,type ="p" ,pch=15,bg="black",col="black")
      }
    }
    dev.off()
  }
}

#Save a table of these mutations
if(exists("combined_df")) {
  print("Saving table of all variants")
  write.table(combined_df,file = paste0(output_dir,"/",sampleID,"_mut_table.tsv"),quote = FALSE,sep = "\t",row.names = FALSE)
} else {
  print("No PVV or MAV mutations found")
}
