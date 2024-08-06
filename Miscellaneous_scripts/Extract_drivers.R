#========================================#
# Define function to extract drivers ####
#========================================#

extract_driver_list=function(details,
                             driver_list #Vector of driver genes
                             ) {
  library(dplyr)
  details$coding_change <- ifelse(details$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                      "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                      "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                      "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                      "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                      "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss")|
                                    grepl("splice_site_variant",details$Type),
                                  "Coding change",
                                  "no")
  
  details$coding_change_driver<-ifelse(details$coding_change=="Coding change" & details$Gene%in%driver_list,"yes","no")
  details$variant_ID=paste(details$Gene, details$Protein, sep = " ")
  
  return(details%>%filter(coding_change_driver=="yes")%>%dplyr::select(mut_ref,Chrom,Pos,Ref,Alt,Mut_type,node,Gene,Transcript,RNA,Protein,variant_ID))
}

#========================================#
# Parse the options list (using optparse) ####
#========================================#
#Load up the filtered mutations file and a list of driver genes to pull out
load("~/R_work/Prolonged_persistence_of_DNA_lesions/Data/input_data/EM/annotated_mut_set_PX001_2_01_standard_rho01")
driver_list=readLines("~/R_work/reference_files/blood_driver_list.txt")

#Run the function to get the list of mutations in genes of interest
extract_driver_list(details=filtered_muts$COMB_mats.tree.build$mat,driver_list = driver_list)
