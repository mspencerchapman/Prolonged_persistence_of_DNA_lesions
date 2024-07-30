#----------------Lesion type simulations----------------

library(dplyr)

##Define the key function for calculating outcome probabilites
#Easily editable to get different info from this
#e.g. the proportion of persistent lesions that will be recognized as such by detection of PVVs/ MAVs

calculated_PVV_to_MAV_ratio=function(Pairing_probabilities,
                                     Two_to_three_replication_ratio) {
  
  #Define the different variables
  Pairing_possibilities=c("G*/C","G*/T","G*/A","G*/G")
  Correct_pairing="G*/C"
  Possible_outcomes=c("No mutation","Mutation","Simple MAV","Separated MAV","Separated MAV (triallelic)","PVV")
  Recognized_persistent_lesion_outcomes=c("Simple MAV","Separated MAV","Separated MAV (triallelic)","PVV")
  
  ## CREATE THE 'TWO REPLICATION' DATA FRAME
  Two_replications_df=expand.grid(Pairing_possibilities,Pairing_possibilities)
  colnames(Two_replications_df)=c("Replication_1","Replication_2")
  Two_replications_df$Result<-apply(Two_replications_df,1,function(vec) {
    if(vec[1]==vec[2]&vec[1]==Correct_pairing) {
      return("No mutation")
    } else if(vec[1]==vec[2]&vec[1]!=Correct_pairing) {
      return("Mutation")
    } else if(vec[1]==Correct_pairing|vec[2]==Correct_pairing) {
      return("Mutation")
    } else if(vec[1]!=Correct_pairing & vec[2]!=Correct_pairing & vec[1]!=vec[2]) {
      return("Simple MAV")
    }
  })
  
  Two_replications_df$Probability<-apply(Two_replications_df,1,function(vec) {
    prod(Pairing_probabilities[vec[1:2]])
  })
  
  ## CREATE THE 'THREE REPLICATION' DATA FRAME
  Three_replications_df=data.frame("Replication_1"=rep(Pairing_possibilities,each=16),
                                   "Replication_2"=rep(Two_replications_df$Replication_1,times=4),
                                   "Replication_3"=rep(Two_replications_df$Replication_2,times=4))
  
  Three_replications_df$Result<-apply(Three_replications_df,1,function(vec) {
    if(sum(vec==Correct_pairing)==3) {
      return("No mutation")
    } else if(vec[1]==Correct_pairing & vec[2]==vec[3] & vec[2]!=Correct_pairing) {
      return("Mutation")
    } else if(sum(vec==Correct_pairing)==2) {
      return("Mutation")
    } else if(vec[1]==Correct_pairing & vec[2]!=vec[3]) {
      return("Simple MAV")
    } else if(sum(vec==Correct_pairing)==1 & length(unique(vec))==2 & vec[2]!=vec[3]){
      return("PVV")
    } else if(sum(vec==Correct_pairing)==1 & length(unique(vec))==3 & vec[2]!=vec[3]){
      return("Separated MAV")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==3){
      return("Separated MAV (triallelic)")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==1){
      return("Mutation")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==2 & vec[2]!=vec[3]){
      return("Separated MAV")
    } else if(sum(vec==Correct_pairing)==0 & length(unique(vec))==2 & vec[2]==vec[3]){
      return("Simple MAV")
    } else {
      return("Undefined")
    }
  })
  
  Three_replications_df$Probability<-apply(Three_replications_df,1,function(vec) {
    prod(Pairing_probabilities[vec[1:3]])
  })
  
  ##Create outcome summaries
  Two_replications_summary<-sapply(Possible_outcomes,function(outcome) {Two_replications_df%>%filter(Result==outcome)%>%pull(Probability)%>%sum()})
  Three_replications_summary<-sapply(Possible_outcomes,function(outcome) {Three_replications_df%>%filter(Result==outcome)%>%pull(Probability)%>%sum()})
  
  outcomes_summary=data.frame(Outcome=Possible_outcomes,
                              Two_replications=Two_replications_summary,
                              Three_replications=Three_replications_summary,
                              Totals=sapply(Possible_outcomes,function(outcome) {Two_to_three_replication_ratio*Two_replications_summary[outcome] + Three_replications_summary[outcome]})
  )
  
  outcomes_summary$Totals_normalized=outcomes_summary$Totals/sum(outcomes_summary$Totals)
  
  outcomes_summary$Recognized_as_persistent_lesions=sapply(Possible_outcomes,function(outcome) {
    if(outcome%in%Recognized_persistent_lesion_outcomes) {
      return(outcomes_summary%>%filter(Outcome==outcome)%>%pull(Totals))
    } else {
      return(0)
    }
  })
  
  outcomes_summary$Recognized_as_persistent_lesions_normalized=outcomes_summary$Recognized_as_persistent_lesions/sum(outcomes_summary$Recognized_as_persistent_lesions)
  
  MAV_total=sum(outcomes_summary$Recognized_as_persistent_lesions_normalized[grepl("MAV",outcomes_summary$Outcome)])
  PVV_total=sum(outcomes_summary$Recognized_as_persistent_lesions_normalized[grepl("PVV",outcomes_summary$Outcome)])
  
  PVV_to_MAV_ratio=PVV_total/MAV_total
  return(PVV_to_MAV_ratio)
}

##
Pairing_probabilities=c("G*/C"=0.49,
                        "G*/T"=0.49,
                        "G*/A"=0.005,
                        "G*/G"=0.005)


calculated_PVV_to_MAV_ratio(Pairing_probabilities=Pairing_probabilities,
                            Two_to_three_replication_ratio=3)




Transition_to_transversion_ratio=100
Pairing_probabilities_unnormalized=c("G*/C"=Transition_to_transversion_ratio,
                                     "G*/T"=Transition_to_transversion_ratio,
                                     "G*/A"=1,
                                     "G*/G"=1)

Pairing_probabilities=Pairing_probabilities_unnormalized/sum(Pairing_probabilities_unnormalized)

outcomes_df=data.frame(Transition_to_transversion_ratio=10^seq(-3,3,0.1))

all_outcomes_df<-lapply(c(1,10,50),function(Two_to_three_replication_ratio) {
  outcomes_df=data.frame(Transition_to_transversion_ratio=10^seq(-3,3,0.1))
  outcomes_df$PVV_to_MAV_ratio=sapply(outcomes_df$Transition_to_transversion_ratio,function(Transition_to_transversion_ratio) {
    Pairing_probabilities_unnormalized=c("G*/C"=Transition_to_transversion_ratio,
                                         "G*/T"=Transition_to_transversion_ratio,
                                         "G*/A"=1,
                                         "G*/G"=1)
    
    Pairing_probabilities=Pairing_probabilities_unnormalized/sum(Pairing_probabilities_unnormalized)
    
    calculated_PVV_to_MAV_ratio(Pairing_probabilities=Pairing_probabilities,
                                Two_to_three_replication_ratio=Two_to_three_replication_ratio)
  })
  return(outcomes_df%>%mutate(Two_to_three_replication_ratio=Two_to_three_replication_ratio))
})%>%dplyr::bind_rows()

breaks=c(0.01,0.1,1,10,100,1000)
PVV_to_MAV_ratio<-all_outcomes_df%>%
  ggplot(aes(x=Transition_to_transversion_ratio,y=PVV_to_MAV_ratio,col=factor(Two_to_three_replication_ratio)))+
  geom_line()+
  theme_classic()+
  scale_x_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  scale_y_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  geom_hline(yintercept=5,linetype=2)+
  geom_vline(xintercept=240,linetype=2)+
  labs(x="Transition to transversion ratio of TLS",
       y="PVV:MAV ratio of resulting\n(recognized) persistent lesions",
       col="Ratio of lesions persisting across only 1 node,\ncompared to those persisting across 2 or more")+
  my_theme

ggsave(filename=paste0(root_dir,"Rebuttal_plots/PVV_to_MAV_ratio_simulated.pdf"),plot=PVV_to_MAV_ratio,width = 5,height=3)
