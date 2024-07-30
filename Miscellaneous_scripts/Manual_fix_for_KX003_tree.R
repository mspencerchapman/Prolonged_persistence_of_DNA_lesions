#CODE to manually fix the dodgy clade in the KX003 tree

library(ape)
tree=read.tree("~/R_work/Prolonged_persistence_of_DNA_lesions/Data/input_data/EM/tree_KX003_5_01_standard_rho01.tree")
lesion_node=472
sub_tree=extract.clade(phy=tree,node=lesion_node)
sub_tree$edge.length<-NULL
plot(sub_tree,direction="downwards")
correct_subtree=read.tree(text = "((PD43974io_new:1,PD43974n2_new:1):1,(PD43974fs2_new:1,PD43974mq_new:1):1);")
plot(correct_subtree,direction="downwards")

new_tree<-ape::bind.tree(tree,correct_subtree,where=lesion_node)
new_tree<-drop.tip(new_tree,tip=sub_tree$tip.label)
new_tree$tip.label<-gsub("_new","",new_tree$tip.label)
plot(new_tree,direction="downwards")
write.tree(new_tree,paste0(lesion_seg_input_dir,"/EM/tree_KX003_5_01_standard_rho01.tree"))

##CAN NOW RERUN THE LESION SEGREGRATION SCRIPT USING THE NEW TREE##