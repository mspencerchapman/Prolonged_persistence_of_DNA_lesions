##Using Nicola Roberts' code on assessing genomic positions

#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("dplyr","tidyr","ggplot2","stringr","readr","tidyr")
bioconductor_packages=c("MutationalPatterns","GenomicFeatures","BSgenome.Hsapiens.UCSC.hg19")


for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(package)
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("nrmisc", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/nrmisc")
  library("nrmisc",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Now define the key function for getting the properties (and some plotting functions)
#----------------------------------

#If no 'control positions' are used, then the whole callable genome is the default control
get_hg19_properties <- function(gpos, location='farm',control_gpos=NULL){
  require(GenomicFeatures)
  require(stringr)
  require("BSgenome.Hsapiens.UCSC.hg19")
  
  # all gpos elements must have width 1.
  if(any(width(gpos)!=1)){stop('All gpos elements must have width 1.')}
  # gpos must have mcol called "hist"
  if(!"hist" %in% colnames(mcols(gpos))) {stop("gpos must have column 'hist'")}
  
  gpos <- sort(gpos, ignore.strand=TRUE)
  
  # check location value, and set root accordingly
  location <- tolower(location)
  if (!location %in% c('farm', 'local')){stop("location must be one of 'farm' or 'local'")}
  if (location=='farm'){
    root <- '/lustre/scratch126/casm/team154pc/ms56/reference_files/nr3/'
  } else {
    root <- '/Volumes/nfs_home/data/'
  }
  
  # genome property files - list by property type
  in_dir <- file.path(root, 'results/GRanges/')
  f_in <- list.files(in_dir)
  f_in <- split(f_in, sapply(strsplit(f_in, "_"), `[`, 1))
  
  # don't include plain g_quadruplex_GR.RData, just use g4 subset with loops <=4
  f_in <- f_in[-which(names(f_in)=='g')]
  
  # load callable genome for ECDF calc
  load(file.path(root, 'results/callable_genome.RData'))
  gpos<-subsetByOverlaps(gpos,callable_genome)
  
  for (i in 1:length(f_in)){
    
    if (length(f_in[[i]])==1) {
      load(file.path(in_dir, f_in[[i]]))
      pname <- gsub(".RData", "", f_in[[i]])
      cat(pname,sep="\n")
      prop <- get(pname)
      rm(list=pname)
      
      suppressWarnings(seqinfo(prop) <- seqinfo(gpos))
      prop <- trim(prop)
      
      ovl <- findOverlaps(gpos, prop)
      
      # only use one metric
      if (ncol(mcols(prop))==1) j <- 1 else {
        if (pname %in% c("LAD_GR", "gene_GR")) {
          j <- which(names(mcols(prop))=="dens_1e6")
        } else if (pname %in% c("cruciform_inverted_rep_GR", "short_tandem_rep_GR")){
          j <- which(names(mcols(prop))=="dens_3e3")
        } else  {
          j <- which(names(mcols(prop))=="dist_log10")
        }
      }
      
      cname <- gsub(" ", "_", paste(gsub("_GR", "", pname), names(mcols(prop))[j]))
      mcols(gpos)[,cname] <- mcols(prop)[subjectHits(ovl),j]
      
      # get quantile value - callable overlap (force 1kb tiles)
      kbs <- unlist(tile(prop, width=1e3))
      ovl <- findOverlaps(kbs, prop)
      kbs$value <- mcols(prop)[subjectHits(ovl),j]
      
      if(is.null(control_gpos)) {
        prop <- subsetByOverlaps(kbs, callable_genome)
        p_ecdf <- ecdf(jitter(prop$value))
        mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
        rm(prop, kbs, j, p_ecdf, ovl, cname, pname)
      } else {
        #New code to use the full mutation set as 'control' to calculate background quantiles
        overlaps <- findOverlaps(kbs, control_gpos)
        temp<-control_gpos
        temp$value[overlaps@to]<-kbs$value[overlaps@from]
        p_ecdf <- ecdf(jitter(temp$value))
        mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
        rm(kbs, j, p_ecdf, ovl, overlaps, cname, pname,temp)
      }
      
      
    } else {
      for(k in seq_along(f_in[[i]])){
        fname <- f_in[[i]][k]
        pname <- gsub(".RData", "", fname)
        tname <- strsplit(pname, "_")[[1]][2]
        
        if (!tname %in% gpos$hist) next
        
        cat(pname,sep="\n")
        load(file.path(in_dir, fname))
        prop <- get(pname)
        suppressWarnings(seqinfo(prop) <- seqinfo(gpos))
        prop <- trim(prop)
        rm(list=pname)
        
        ovl <- findOverlaps(gpos[gpos$hist==tname], prop)
        cname <- names(f_in)[i]
        
        mcols(gpos)[gpos$hist==tname,cname] <- prop$value[subjectHits(ovl)]
        
        # get quantile value - callable overlap (force 1kb tiles)
        kbs <- unlist(tile(prop, width=1e3))
        ovl <- findOverlaps(kbs, prop)
        kbs$value <- mcols(prop)[subjectHits(ovl),1] #edited this to 1 as there is no 'j' in this section
        
        if(is.null(control_gpos)) {
          prop <- subsetByOverlaps(kbs, callable_genome)
          p_ecdf <- ecdf(jitter(prop$value))
          mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
          rm(prop, kbs, p_ecdf, ovl, cname, pname)
        } else {
          #New code to use the full mutation set as 'control' to calculate background quantiles
          overlaps <- findOverlaps(kbs, control_gpos)
          temp<-control_gpos
          temp$value[overlaps@to]<-kbs$value[overlaps@from]
          p_ecdf <- ecdf(jitter(temp$value))
          mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
          rm(kbs, p_ecdf, ovl, overlaps, cname, pname,temp)
        }
      }
    }
  }
  return(gpos)
}

get_quantile_metrics=function(df,
                              signif_levels=c(1e-10,1e-5,1e-2)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  qs <- df[,union(1:6, grep('^q_', colnames(df)))]
  
  ##Calculate significance of each metric
  bin_data=function(num_vec,bins) {
    x=vector(mode="numeric",length=(length(bins)-1))
    for(i in 1:(length(bins)-1)) {
      x[i]<-sum(num_vec>bins[i] & num_vec<=bins[i+1])
    }
    return(x)
  }
  all_metrics=grep("^q_",colnames(qs),value = T)
  metric_res=data.frame(metric=all_metrics)
  metric_res$p.value=sapply(all_metrics,function(metric) {
    #cat(metric,sep="\n")
    nbin=10; bins=seq(0,1,1/nbin)
    dat_binned=bin_data(num_vec = qs[[metric]][!is.na(qs[[metric]])],bins=bins)
    chisq_res=chisq.test(x = dat_binned)
    return(chisq_res$p.value)
  })
  metric_res$q.value=p.adjust(p=metric_res$p.value,method="BH")
  metric_res$signif=ifelse(metric_res$q.value<signif_levels[1],"***",ifelse(metric_res$q.value<signif_levels[2],"**",ifelse(metric_res$q.value<signif_levels[3],"*","n.s.")))
  return(metric_res)
}

plot_quantile_metrics=function(df,
                               metrics_to_plot="signif",
                               signif_levels=c(1e-10,1e-5,1e-2)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  
  my_theme=theme_classic(base_family="Helvetica")+theme(text=element_text(size=7,family="Helvetica"),
                                                        axis.text=element_text(size=5,family="Helvetica"),
                                                        strip.text = element_text(size=6,family="Helvetica"),
                                                        legend.key.height = unit(0.25, 'cm'),
                                                        legend.key.width = unit(0.25, 'cm'),
                                                        legend.title = element_text(size=6),
                                                        legend.text = element_text(size=5,family="Helvetica"))
  
  qs <- df[,union(1:6, grep('^q_', colnames(df)))]
  
  ##Calculate significance of each metric
  bin_data=function(num_vec,bins) {
    x=vector(mode="numeric",length=(length(bins)-1))
    for(i in 1:(length(bins)-1)) {
      x[i]<-sum(num_vec>bins[i] & num_vec<=bins[i+1])
    }
    return(x)
  }
  all_metrics=grep("^q_",colnames(qs),value = T)
  metric_res=data.frame(metric=all_metrics)
  metric_res$p.value=sapply(all_metrics,function(metric) {
    #cat(metric,sep="\n")
    nbin=10; bins=seq(0,1,1/nbin)
    dat_binned=bin_data(num_vec = qs[[metric]][!is.na(qs[[metric]])],bins=bins)
    chisq_res=chisq.test(x = dat_binned)
    return(chisq_res$p.value)
  })
  metric_res$q.value=p.adjust(p=metric_res$p.value,method="BH")
  metric_res$signif=ifelse(metric_res$q.value<signif_levels[1],"***",ifelse(metric_res$q.value<signif_levels[2],"**",ifelse(metric_res$q.value<signif_levels[3],"*","n.s.")))
  
  if(metrics_to_plot[1]=="signif") {
    significant_metrics=metric_res%>%filter(signif!="n.s.")%>%pull(metric)
  } else if(metrics_to_plot[1]=="all"){
    significant_metrics<-all_metrics
  } else {
    significant_metrics<-metrics_to_plot
  }
  p1<-qs[,-c(1:6)]%>%
    gather(key="Metric",value="Quantile")%>%
    filter(Metric%in%significant_metrics)%>%
    mutate(Metric=gsub("^q_","",Metric))%>%
    ggplot(aes(x=Quantile))+
    geom_density(fill="lightblue")+
    facet_wrap(~Metric,ncol=5)+
    theme_classic()+
    my_theme
  
  return(p1)
}

combined_plot=function(df_list,metrics_to_plot=NULL) {
  library(dplyr)
  library(ggplot2)
  my_theme=theme_classic(base_family="Helvetica")+theme(text=element_text(size=7,family="Helvetica"),
                                                        axis.text=element_text(size=5,family="Helvetica"),
                                                        strip.text = element_text(size=6,family="Helvetica"),
                                                        legend.key.height = unit(0.25, 'cm'),
                                                        legend.key.width = unit(0.25, 'cm'),
                                                        legend.title = element_text(size=6),
                                                        legend.text = element_text(size=5,family="Helvetica"))
  
  datasets=names(df_list)
  
  if(is.null(metrics_to_plot)) {
    all_quantile_metrics<-Map(df=df_list,dataset=datasets,function(df,dataset) {
      df<-df[,grep('^q_', colnames(df))]%>%
        gather(key="Metric",value="Quantile")%>%
        mutate(Metric=gsub("^q_","",Metric))%>%
        mutate(dataset=dataset,.before=1)
      return(df)
    })%>%dplyr::bind_rows()
  } else {
    all_quantile_metrics<-Map(df=df_list,dataset=datasets,function(df,dataset) {
      df<-df[,grep('^q_', colnames(df))]%>%
        gather(key="Metric",value="Quantile")%>%
        filter(Metric%in%metrics_to_plot)%>%
        mutate(Metric=gsub("^q_","",Metric))%>%
        mutate(dataset=dataset,.before=1)
      return(df)
    })%>%dplyr::bind_rows()
  }
  
  
  
  plot<-ggplot(all_quantile_metrics,aes(x=Quantile))+
    geom_density(fill="lightblue")+
    facet_grid(Metric~dataset)+
    scale_x_continuous(breaks=c(0,1))+
    theme_classic()+
    my_theme+
    theme(strip.text.y = element_text(angle=0),
          axis.text.x=element_blank())
  
  plot
}

# Set the ggplot2 theme and palettes for plotting-----
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
# Create the GRanges objects of the mutation sets----
#----------------------------------

##Create the GRanges object of all HCT mutations to use as a 'control'
all_HCT_mutations_path<-"/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/HCT_sites_per_branch.txt"
all_HCT_mutations<-read.delim(all_HCT_mutations_path,stringsAsFactors = F)%>%
  filter(Mut_type=="SNV")%>%
  mutate(Chrom=paste0("chr",Chrom))
all_HCT_mutations_GR<-GenomicRanges::makeGRangesFromDataFrame(df=all_HCT_mutations,
                                        seqnames.field="Chrom",
                                        start.field="Pos",
                                        end.field="Pos")
all_HCT_mutations_GR$hist<-"Myeloid"

##Now do the same for the blood PVVs
confirmatory_muts<-read.delim("/lustre/scratch126/casm/team154pc/ms56/confirmatory_mutations.tsv",stringsAsFactors = F)%>%
  mutate(Chrom=paste0("chr",Chrom))
confirmatory_muts_GR<-GenomicRanges::makeGRangesFromDataFrame(df=confirmatory_muts,
                                                              seqnames.field="Chrom",
                                                              start.field="Pos",
                                                              end.field="Pos")
confirmatory_muts_GR$hist<-"Myeloid"

#----------------------------------
# Run function over the gpos object-----
#----------------------------------

confirmatory_muts_with_properties=get_hg19_properties(gpos=confirmatory_muts_GR,control_gpos = all_HCT_mutations_GR)
confirmatory_muts_with_properties_df=as.data.frame(confirmatory_muts_with_properties)

metrics_df=get_quantile_metrics(df=confirmatory_muts_with_properties_df,
                     signif_levels = c(1e-3,1e-2,0.1))%>%arrange(q.value)

readr::write_csv(metrics_df%>%mutate(metric=gsub("^q_","",metric))%>%arrange(p.value),file="PVV_association_metrics.csv")

PVVs_plot<-plot_quantile_metrics(df=PVVs_with_properties_df,
                      metrics_to_plot = "all",
                      signif_levels = c(1e-3,1e-2,0.1))

confirmatory_muts_signif_plot<-plot_quantile_metrics(df=confirmatory_muts_with_properties_df,
                                 metrics_to_plot = "signif",
                                 signif_levels = c(1e-3,1e-2,0.1))

ggsave("Confirmatory_muts_signif_metrics.pdf",plot = confirmatory_muts_signif_plot,width=3,height=2.5)

