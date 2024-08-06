#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("stringr","ape","remotes")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}
options(stringsAsFactors = FALSE)

#========================================#
# Set file paths ####
#========================================#
tree_file_path=""
filtered_muts_file=""

#========================================#
# Load data, reassign mutations based on tree provided, and resave ####
#========================================#

#Load up the objects, and pull out the details data frame and the NV and NR matrices
print("Loading the tree and converting to multifurcating structure")
tree <- di2multi(read.tree(tree_file_path))
print("Loading the mutations matrix")
load(filtered_muts_file)

details=filtered_muts$COMB_mats.tree.build$mat
NV = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
NR = as.matrix(filtered_muts$COMB_mats.tree.build$NR)

re_run=T #Forces reassignment even if everything seems to match up

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
