library(cpgAccessData)
library(data.table)
library(GenomicRanges)
library(ggplot2)

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

# read in CpG methylation data from cord blood donor
CpgConfig$DATA="/lustre/scratch126/casm/team273jn/methylation/cpgdata/post_qc/latest"
CpgConfig$METADATA="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/methylation/metadata"
cb_meth=cpg_get_cpgdata("PD45517", chrom=c(1:22,"X"))

cb_meth_pseudo_vct <- rowSums(cb_meth$MF + cb_meth$MR)/rowSums(cb_meth$MF + cb_meth$MR + cb_meth$UF + cb_meth$UR)
cb_meth_pseudo_dt <- data.table(Chrom = cb_meth$details$Chrom,
                                Pos = cb_meth$details$Pos,
                                mm = cb_meth_pseudo_vct)
# read in CpG density data
CpgConfig$DATA="/lustre/scratch126/casm/team273jn/methylation/build_example/post_qc"
CpgConfig$METADATA="/lustre/scratch126/casm/team273jn/methylation/build_example/metadata/"
cpg_dens=cpg_get_cpgdata("PD5163",chrom=c(1:22, "X"))

nrow(cpg_dens$details)

# read in PVVs
root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
mutations<-read.delim(paste0(root_dir,"/lesion_segregation/mutations_filtered.tsv"))
ABC_Sample_IDs=c("KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01")
PVV_blood_table<-mutations%>%
  filter(Type=="PVV" & cat%in%c("Adult_HSPC","Chemo_HSPC","Foetal_HSPC") & Class=="PASS" & Sub1=="C>T" & lesion_duration<200 & grepl("\\[C>T\\]T",mut_profile_1))

pvv_dt <- data.table(PVV_blood_table[,1])
setnames(pvv_dt, new = colnames(PVV_blood_table)[1])
pvv_dt[, c("Chrom", "Pos") := tstrsplit(Chrom_pos, "-")]
pvv_dt[, Pos := as.numeric(Pos)]
pvv_dt[, chr_pos := paste0(Chrom, "_", Pos)]
pvv_dt[, start := Pos-1]
pvv_dt[, end := Pos]
setkey(pvv_dt, Chrom, start, end)


# get 1000bp bins on genome
maxpos_dt <- data.table(cpg_dens$details)[, (maxpos=max(Pos)),by=Chrom]
roundUp <- function(x,to=1000){ to*(x%/%to + as.logical(x%%to))}
bins_list_dt <- mapply(FUN=function(chrom, maxpos){ bins_dt <- data.table(Chrom = chrom,
                                                                          start = seq(1, maxpos, 1000),
                                                                          end = c(seq(1000, maxpos, 1000), roundUp(maxpos)))
return(bins_dt)}, chrom = maxpos_dt[["Chrom"]], maxpos = maxpos_dt[["V1"]], SIMPLIFY = FALSE)
bins_dt <- rbindlist(bins_list_dt)
setkey(bins_dt, Chrom, start, end)

# preprocess cpg sites
dt_cpg <- data.table(cpg_dens$details)
setkey(dt_cpg, Chrom, Pos)
dt_cpg[, start := Pos-1]
dt_cpg[, end := Pos]
setkey(dt_cpg, Chrom, start, end)

# overlaps between CpG sites and bins
dt_binned <- foverlaps(dt_cpg, bins_dt, mult = "first")

# add per CpG pseudo-bulk mean methylation
dt_binned2 <- merge(dt_binned, cb_meth_pseudo_dt, by = c("Chrom", "Pos"))

# group by chrom, bin start, bin end and cpg_context and get mean methylation and cpg_count
dt_dens <- dt_binned2[, list(mean_meth = mean(mm),
                             cpg_count = .N),
                      by=.(Chrom, start, end, cpg_context)]
# add PVV count
# overlap PVV with 1K bp bins
pvv_bins <- foverlaps(pvv_dt, bins_dt, mult = "first")
# number of PVVs per 1K bp bin
pvv_count_bin <- pvv_bins[, list(pvv_count=.N), by=.(Chrom, start, end)]
table(pvv_count_bin$pvv_count)


# merge PVV with CpG density by bin
dt_full <- merge(dt_dens, pvv_count_bin, by=c("Chrom", "start", "end"), all=TRUE)
setnafill(dt_full, fill=0, cols=c("pvv_count", "cpg_count","mean_meth"))
# split CpG density into quantiles
dt_full[, quantile:= cut(cpg_count, quantile(cpg_count, probs = 0:5/5),
                         labels = FALSE, include.lowest = TRUE)]


# plot
plot1 <- ggplot(dt_full[!is.na(mean_meth),], aes(x=mean_meth)) +
  facet_grid(pvv_count~cpg_context, scales="free_y") +
  geom_histogram()+
  my_theme

ggsave(filename="temp.pdf",plot1,width=5,height=3)

#
feature_proportions<-as.data.frame(table(dt_full$cpg_context))%>%dplyr::rename("Feature"=Var1,"N"=Freq)%>%mutate(Prop=N/sum(N))
pvv_feature_proportions<-as.data.frame(table(dt_full[pvv_count==1,]$cpg_context))%>%dplyr::rename("Feature"=Var1,"N"=Freq)%>%mutate(Prop=N/sum(N))

prop.test(x=pvv_feature_proportions$N,n=feature_proportions$N)

do_chi_sq_by_feature=function(df,feature,bins=5) {
  df_ctrl<-df%>%
    filter(cpg_context==feature)%>%
    mutate(quant = cut(mean_meth, unique(quantile(mean_meth, seq(0, 1, 1/bins))), labels = FALSE))
  
  df_test<-df_ctrl%>%
    filter(cpg_context==feature & pvv_count>0)
  if(nrow(df_test)<2) {stop(return("Insufficient PVVs"))}
  test_set=table(df_test$quant)[as.character(1:bins)]
  test_set[is.na(test_set)]<-0
  return(chisq.test(test_set))
}

lapply(c("cpg_islands","cpg_shores","cpg_inter","cpg_shelves"), function(feature) {
  do_chi_sq_by_feature(dt_full,feature=feature,bins=10)
})


bins=10
dt_full<-dt_full%>%
  mutate(mean_meth_cent = cut(mean_meth, unique(quantile(mean_meth, seq(0, 1, 1/bins))), labels = FALSE))%>%
  mutate(cpg_count_cent = cut(cpg_count, unique(quantile(cpg_count, seq(0, 1, 1/bins))), labels = FALSE))


PVV_CpG_feature_2D<-dt_full%>%
  filter(pvv_count>0)%>%
  ggplot(aes(x=mean_meth_cent,y=cpg_count_cent))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white",linewidth=0.2)+
  scale_x_continuous(limits=c(0,bins+1),breaks=seq(0,10,1))+
  scale_y_continuous(limits=c(0,bins+1),breaks=seq(0,10,1))+
  theme_classic()+
  my_theme+
  labs(x="Mean methylation quantiles",y="CpG count quantiles")

ggsave(filename = "PVV_CpG_features_2D.pdf",PVV_CpG_feature_2D,width=3,height=3)

