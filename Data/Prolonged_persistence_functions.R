#Function to subset the list of matrices
list_subset = function(list, select_vector) {
  for (i in 1:length(list)) {
    if(!is.null(dim(list[[i]]))) {
      list[[i]] = list[[i]][select_vector,]
    }
  }
  return(list)
}

#Function to create pval matrix based on the NV and NR matrix
pval_matrix = function(COMB_mats) {
  cat("Starting pval matrix generation\n")
  COMB_mats$NR[COMB_mats$NR == 0] <- 1 
  pval_mat <- matrix(0, nrow = nrow(COMB_mats$NV), ncol = ncol(COMB_mats$NV))
  if(COMB_mats$gender == "male") {
    for(i in 1:nrow(COMB_mats$NV)) {
      for (j in 1:ncol(COMB_mats$NR)) {
        if (!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {pval_mat[i,j] <- binom.test(COMB_mats$NV[i,j], COMB_mats$NR[i,j], p = 0.5, alternative = "less")$p.value}
        else {pval_mat[i,j] <- binom.test(COMB_mats$NV[i,j], COMB_mats$NR[i,j], p = 0.95, alternative = "less")$p.value}
      }
      if (i %% 1000 == 0) {print(i)}
    }
  } else if(COMB_mats$gender == "female") {
    for(i in 1:nrow(COMB_mats$NV)) {
      for (j in 1:ncol(COMB_mats$NR)) {
        pval_mat[i,j] <- binom.test(COMB_mats$NV[i,j], COMB_mats$NR[i,j], p = 0.5, alternative = "less")$p.value
      }
      if (i %% 1000 == 0) {print(i)}
    }
  }
  cat("Completed pval matrix generation\n")
  return(pval_mat)
}

#mat object needs the Chrom column with chromosome.  Returns the pval vector for each mutation.
germline.binomial.filter = function(COMB_mats){
  cat("Starting the germline binomial filter\n")
  XY_chromosomal = COMB_mats$mat$Chrom %in% c("X","Y")
  autosomal = !XY_chromosomal
  
  if(COMB_mats$gender=="female"){
    NV_vec = rowSums(COMB_mats$NV) #create vector of COMBINED variant reads across all samples
    NR_vec = rowSums(COMB_mats$NR) #create vector of COMBINED depth across all samples
    pval = rep(1,length(NV_vec))
    #For loop to calculate whether each one is likely to have come from a binomial distribution where true probability is 0.5 (i.e. would be expected for autosomal heterozygous germline variant)
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      } 
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # IF MALE - need to separate off the XY chromosomes, as expected probability of germline variant is now close to 1.
  if(COMB_mats$gender=="male"){
    pval=rep(1,nrow(COMB_mats$NV))
    NV_vec = rowSums(COMB_mats$NV)[autosomal]
    NR_vec = rowSums(COMB_mats$NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(COMB_mats$NV)[XY_chromosomal]
    NR_vec = rowSums(COMB_mats$NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  cat("Completed the germline binomial filter\n")
  return(pval)
}

#THE BETA-BINOMIAL FUNCTIONS
require(VGAM)
estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(COMB_mats){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  cat("Starting the beta-binomial filter\n")
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  rho_est = rep(NA,nrow(COMB_mats$NR))
  for (k in 1:nrow(COMB_mats$NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(COMB_mats$NV[k,]),
                                  NR_vec=as.numeric(COMB_mats$NR[k,]))
    if (k%%1000==0){
      print(k)
    }
  }
  cat("Completed the beta-binomial filter\n")
  return(rho_est)
}

#FILTER VARIANTS WITH LOW VAF AMONGST POSITIVE SAMPLES
#Filters mutations present in multiple samples that are (on aggregate) unlikely to have come from true somatic mutation distribution (i.e. 0.5 in auto, 1 in XY)
low_vaf_in_pos_samples_dp2 = function(COMB_mats, define_pos = 2) {
  cat("Starting the low_vaf_in_pos_samples_dp2 filter\n")
  if(all(c(nrow(COMB_mats$mat) == nrow(COMB_mats$NV),dim(COMB_mats$NV) == dim(COMB_mats$NR)))) {print("Input matrices are of correct dimensions")}
  pval=rep(0,nrow(COMB_mats$mat))
  if(COMB_mats$gender == "male") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        if (COMB_mats$mat$Chrom[k] %in% c("X","Y")) {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        }
        else {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
        if(k %% 1000 == 0) {print(k)}
      }
    }
  } else if(COMB_mats$gender == "female") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        if(k %% 1000 == 0) {print(k)}
      }
    }
  }
  cat("Completed the low_vaf_in_pos_samples_dp2 filter\n")
  return(pval)
}


low_vaf_in_pos_samples_dp3 = function(COMB_mats, define_pos = 3) {
  cat("Starting the low_vaf_in_pos_samples_dp3 filter\n")
  if(all(c(nrow(COMB_mats$mat) == nrow(COMB_mats$NV),dim(COMB_mats$NV) == dim(COMB_mats$NR)))) {print("Input matrices are of correct dimensions")}
  pval=rep(0,nrow(COMB_mats$mat))
  if(COMB_mats$gender == "male") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        if (COMB_mats$mat$Chrom[k] %in% c("X","Y")) {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        }
        else {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
        if(k %% 1000 == 0) {print(k)}
      }
    }
  } else if(COMB_mats$gender == "female") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        if(k %% 1000 == 0) {print(k)}
      }
    }
  }
  cat("Completed the low_vaf_in_pos_samples_dp3 filter\n")
  return(pval)
}

get_mean_depth = function(COMB_mats) {
  mean_depth = rowMeans(COMB_mats$NR)
  cat("Completed the mean depth filter\n")
  return(mean_depth)
}

get_max_depth_in_pos = function(COMB_mats) {
  cat("Starting the max_depth_in_pos filter\n")
  apply_max_depth_in_pos = function(i, mat_list) {
    if(!mat_list$mat$Chrom[i] %in% c("X","Y") | COMB_mats$gender == "female") {
      max_depth_in_pos_samples <- max(mat_list$NR[i,which(mat_list$NV[i,]>=min_variant_reads_auto)])
    } else {
      max_depth_in_pos_samples <- max(mat_list$NR[i,which(mat_list$NV[i,]>=min_variant_reads_xy)])
    }
    return(max_depth_in_pos_samples)
  }
  max_depth_in_pos_vec = sapply(1:nrow(COMB_mats$NV), apply_max_depth_in_pos, mat_list = COMB_mats)
  cat("Completed the max_depth_in_pos filter\n")
  return(max_depth_in_pos_vec)
}

get_max_pval_in_pos = function(COMB_mats) {
  apply_max_pval_in_pos = function(i, COMB_mats) {
    if(COMB_mats$gender == "male") {
      if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
        max_pval_in_pos_samples <- max(COMB_mats$PVal[i,which(COMB_mats$NV[i,]>=min_variant_reads_auto)])
      } else {
        max_pval_in_pos_samples <- max(COMB_mats$PVal[i,which(COMB_mats$NV[i,]>=min_variant_reads_xy)])
      }
    } else if(COMB_mats$gender == "female") {
      max_pval_in_pos_samples <- max(COMB_mats$PVal[i,which(COMB_mats$NV[i,]>=min_variant_reads_auto)])
    }
    return(max_pval_in_pos_samples)  
  }
  max_pval_in_pos_vec = sapply(1:nrow(COMB_mats$NV), apply_max_pval_in_pos, COMB_mats = COMB_mats)
  return(max_pval_in_pos_vec)
}

get_max_vaf = function(COMB_mats) {
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  vaf_mat = COMB_mats$NV/COMB_mats$NR
  max_vaf = apply(vaf_mat,1,max)
  return(max_vaf)
}

#Removing columns with low coverage, or that are otherwise unwanted
remove_low_coverage_samples = function(COMB_mats,
                                       filter_params=NULL,
                                       min_sample_mean_cov,
                                       other_samples_to_remove = NULL,
                                       min_variant_reads_auto = 3,
                                       min_variant_reads_xy = 2) {
  mean_sample_cov = colMeans(COMB_mats$NR)
  low_cov_names = gsub(pattern = "_DEP", replacement = "", x = names(mean_sample_cov)[mean_sample_cov < min_sample_mean_cov])
  remove_cols <- which(gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats$NV)) %in% c(low_cov_names, other_samples_to_remove))
  if(length(remove_cols) > 0) {
    cat("Removing",length(remove_cols),"samples\n")
    COMB_mats$NV <- COMB_mats$NV[,-remove_cols]
    COMB_mats$NR <- COMB_mats$NR[,-remove_cols]
    COMB_mats$PVal <- COMB_mats$PVal[,-remove_cols]
    
    null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto|(COMB_mats$NV >= min_variant_reads_xy & COMB_mats$mat$Chrom %in% c("X","Y") & COMB_mats$gender == "male")) == 0
    cat(sum(null_remove),"mutations removed as no positives in any remaining samples.\n")
    COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
    if(!is.null(filter_params)) {
      filter_params = filter_params[!null_remove,]
      output = list(COMB_mats=COMB_mats,filter_params=filter_params)
    } else {
      output=COMB_mats
    }
  } else {
    cat("No samples removed\n")
    if(!is.null(filter_params)) {
      output = list(COMB_mats=COMB_mats,filter_params=filter_params)
    } else {
      output=COMB_mats
    }
  }
  return(output)
}

#Functions for filtering from the filter_params and COMB_mats object, setting the desired cut-offs
assess_mean_depth = function(i, COMB_mats, AUTO_low_depth_cutoff, AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff) {
  if(COMB_mats$gender == "male") {
    if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
      ifelse((filter_params$mean_depth[i] >= AUTO_low_depth_cutoff & filter_params$mean_depth[i] <= AUTO_high_depth_cutoff), 1,0)
    } else {
      ifelse(filter_params$mean_depth[i] >= XY_low_depth_cutoff & filter_params$mean_depth[i] <= XY_high_depth_cutoff, 1,0)
    }
  } else if(COMB_mats$gender == "female") {
    ifelse((filter_params$mean_depth[i] >= AUTO_low_depth_cutoff & filter_params$mean_depth[i] <= AUTO_high_depth_cutoff), 1,0)
  }
}

assess_max_depth_in_pos = function(i, COMB_mats, min_depth_auto, min_depth_xy) {
  if(COMB_mats$gender == "male") {
    if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
      ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_auto, 1,0)
    } else {
      ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_xy, 1,0)
    }
  } else if(COMB_mats$gender == "female") {
    ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_auto, 1,0)
  }
}

assess_max_vaf = function(i, COMB_mats, min_vaf_auto, min_vaf_xy) {
  if(COMB_mats$gender == "male") {
    if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
      ifelse(filter_params$max_mut_vaf[i] >= min_vaf_auto, 1,0)
    } else {
      ifelse(filter_params$max_mut_vaf[i] >= min_vaf_xy, 1,0)
    }
  } else if(COMB_mats$gender == "female") {
    ifelse(filter_params$max_mut_vaf[i] >= min_vaf_auto, 1,0)
  }
}

get_filtered_mut_set = function(input_set_ID,
                                COMB_mats,
                                filter_params,
                                gender,
                                
                                ##PARAMETERS FOR FILTERING THE MUTATION SET. These arguments set the cut-offs for selecting the set of "true somatic mutations"
                                #NB. If parameter is set to NA, filter will not be applied.
                                retain_muts = NA, #vector of "mut_refs" (in format Chrom-Pos-Ref-Alt) to retain even if fail the filtering
                                exclude_muts = NA, #vector of "mut_refs" (in format Chrom-Pos-Ref-Alt) to exclude even if pass the filtering
                                germline_pval = -10,  
                                rho = 0.1, #Beta-binomial filter, rho value cut-off (i.e. rho must be > cut-off)
                                mean_depth=NA, #Numeric vector of length 4 in the order (1) AUTO min depth cutoff, (2) AUTO max depth cutoff, (3) XY min depth cutoff, (4) XY max depth cutoff
                                pval_dp2=NA, #Numeric vector of length 1. 
                                pval_dp3=0.01, #Numeric vector of length 1. 
                                min_depth = c(6,4), #Numeric vector of length 2 for (1) AUTO and (2) XY muts. Minimum depth that at least one positive sample must have.
                                min_pval_for_true_somatic = 0.1, #Numeric vector of length 1. At least one sample must have a p-value of >this cutoff for coming from a true somatic distribution.
                                min_vaf = c(0.2,0.8), #Numeric vector of length 2 for (1) AUTO and (2) XY muts. At least one sample must meet the minimum vaf for the mutation.
                                
                                #These arguments decide the thresholds for genotyping each sample. This may be less stringent.
                                min_variant_reads_SHARED = 2, #PARAMETERS FOR DECIDING GENOTYPE FOR EACH SAMPLE FOR EACH RETAINED MUTATION - should be equal to, or more relaxed than the above. At least one parameter must be used.
                                min_pval_for_true_somatic_SHARED = 0.05,
                                min_vaf_SHARED = c(0.2,0.7) #Numeric vector of length 2 for AUTO and XY muts
) {
  #List (1) the names of each of the initial filters that can be applied,
  #(2) the names of their input parameters,
  #(3) the variable name that must be present in the filter_params dataframe if filter is being used.
  #The elements of each must match
  filter_name=c("Germline filter","Beta-binomial filter","Mean depth filter","P-value within positive samples [defining positive as ≥2 variant reads] filter","P-value within positive samples [defining positive as ≥3 variant reads] filter","Minimum depth filter","Minimum p-value for true somatic mutation filter","Minimum VAF filter")
  set_parameter=list(germline_pval,rho,mean_depth,pval_dp2,pval_dp3,min_depth,min_pval_for_true_somatic,min_vaf)
  var_name = c("germline_pval","bb_rhoval","mean_depth","pval_within_pos_dp2","pval_within_pos_dp3","max_depth_in_pos_samples","max_pval_in_pos_sample","max_mut_vaf")
  
  #List the functions that operate on the input parameters to decide if the mutation is a pass (1) or fail (0) for each filter
  test_function=list(function(x) {ifelse(log10(filter_params$germline_pval)<x,1,0)},
                     function(x) {ifelse(filter_params$bb_rhoval > x, 1, 0)},
                     function(x) {sapply(1:nrow(COMB_mats$NV), assess_mean_depth, COMB_mats = COMB_mats, AUTO_low_depth_cutoff = x[1], AUTO_high_depth_cutoff = x[2], XY_low_depth_cutoff = x[3], XY_high_depth_cutoff=x[4])},
                     function(x) {ifelse(filter_params$pval_within_pos_dp2 > x, 1, 0)},
                     function(x) {ifelse(filter_params$pval_within_pos_dp3 > x, 1, 0)},
                     function(x) {sapply(1:nrow(COMB_mats$NV), assess_max_depth_in_pos, COMB_mats = COMB_mats, min_depth_auto = x[1], min_depth_xy = x[2])},
                     function(x) {ifelse(filter_params$max_pval_in_pos_sample > x, 1, 0)},
                     function(x) {sapply(1:nrow(COMB_mats$NV), assess_max_vaf, COMB_mats = COMB_mats, min_vaf_auto = x[1], min_vaf_xy = x[2])}
  )
  
  #CHECK THE INPUT DATA is all present for the set filters
  if(any(!c("NR","NV")%in%names(COMB_mats))|!all.equal(dim(COMB_mats$NR),dim(COMB_mats$NV))) {
    stop(return("COMB_mats object must contain matched NR and NV matrices of equal dimensions"))
  }
  
  if(dim(filter_params)[1]!=dim(COMB_mats$NV)[1]) {
    stop(return("The filter_params data frame and the COMB_mats$NV and NR matrices must contain the same number of mutations"))
  }
  #For any filter that has a set parameter, check that the corresponding variable is included in the filter_params data frame. If not, stop & return an error.
  for(i in 1:length(filter_name)) {
    if(!is.na(set_parameter[[i]][1]) & !any(names(filter_params)==var_name[i])) {
      stop(print(paste(filter_name[i],"needs a variable named",var_name[i],"in the filter_params data frame. Update the filter_params data frame or set the appropriate argument to NA")))
    }
  }
  #If using either 'pval for true somatic' parameter, need to include the PVal matrix in COMB_mats
  if((!is.na(min_pval_for_true_somatic_SHARED)|!is.na(min_pval_for_true_somatic))&!any(names(COMB_mats)=="PVal")) {
    stop(return("If using the 'min_pval_for_true_somatic' or 'min_pval_for_true_somatic_SHARED' parameters, the COMB_mats list must contain a matrix called PVal"))
  }
  
  #Apply each of the filter functions, for those with NULL parameters, a vector of 1's is returned (i.e. all mutations "pass" the filter)
  out=mapply(FUN=function(param,test_function,var_name,filter_name) {
    if(is.na(param[1])) {
      return(rep(1,nrow(filter_params)))
    } else if(!any(names(filter_params)==var_name)){
      stop(return(paste(filter_name,"needs a variable named:",var_name,"in the filter_params data frame")))
    } else {
      return(test_function(param))
    }
  },
  param=set_parameter,
  test_function=test_function,
  var_name=var_name,
  filter_name=filter_name,
  SIMPLIFY = F)
  
  filter_pass=Reduce(cbind,out);rownames(filter_pass)<-rownames(filter_params);colnames(filter_pass)<-var_name #Combine the output & name the rows
  select_muts = (apply(filter_pass,1, function(x) all(x == 1))&!rownames(filter_pass)%in%exclude_muts)|rownames(filter_pass) %in% retain_muts #Select the mutations for output. These must pass ALL the applied filters.
  filter_code = apply(filter_pass, 1, paste, collapse = "-") #Save a "filter_code" vector. This can be used as a quick test for which filter is removing most mutations.
  COMB_mats.tree.build = list_subset(COMB_mats, select_vector = select_muts) #Subset matrices to include only the PASS mutations
  
  #Set the rownames to the mut_ref, and colnames to the sample names (without MTR or DEP)
  rownames(COMB_mats.tree.build$NV) = rownames(COMB_mats.tree.build$NR) = rownames(COMB_mats.tree.build$PVal) <- COMB_mats.tree.build$mat$mut_ref
  colnames = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats.tree.build$NV))
  colnames(COMB_mats.tree.build$NV) = colnames(COMB_mats.tree.build$NR) = colnames(COMB_mats.tree.build$PVal) <- colnames
  
  #BUILD THE GENOTYPE MATRIX
  #Select the "definite positive" samples by setting genotype to 1
  #First build individual matrices, the same dimensions as the NV matrix, where each cell is 1 if it passes that criteria, or 0 if not. If criteria is NULL, set to 1.
  if(!is.na(min_variant_reads_SHARED)) {min_variant_reads_mat <- COMB_mats.tree.build$NV >= min_variant_reads_SHARED} else {min_variant_reads_mat <- 1}
  if(!is.na(min_pval_for_true_somatic_SHARED)) {min_pval_for_true_somatic_mat <- COMB_mats.tree.build$PVal > min_pval_for_true_somatic_SHARED} else {min_pval_for_true_somatic_mat <- 1}
  if(!is.na(min_vaf_SHARED[1]) & gender == "female") {
    depth_no_zero = COMB_mats.tree.build$NR
    depth_no_zero[depth_no_zero == 0] <- 1
    min_vaf_mat <- COMB_mats.tree.build$NV/depth_no_zero > min_vaf_SHARED[1]
  } else if(!is.na(min_vaf_SHARED) & gender == "male") {
    min_vaf_mat = matrix(0, ncol = ncol(COMB_mats.tree.build$NV), nrow = nrow(COMB_mats.tree.build$NV))
    xy_muts = COMB_mats.tree.build$mat$Chrom %in% c("X","Y")
    depth_no_zero = COMB_mats.tree.build$NR
    depth_no_zero[depth_no_zero == 0] <- 1
    min_vaf_mat[xy_muts,] <- COMB_mats.tree.build$NV[xy_muts,]/depth_no_zero[xy_muts,] > min_vaf_SHARED[2]
    min_vaf_mat[!xy_muts,] <- COMB_mats.tree.build$NV[!xy_muts,]/depth_no_zero[!xy_muts,] > min_vaf_SHARED[1]
  } else {min_vaf_mat <- 1}
  
  COMB_mats.tree.build$Genotype_bin = min_variant_reads_mat * min_pval_for_true_somatic_mat * min_vaf_mat
  
  #Select the "not sure" samples by setting genotype to 0.5.  THIS IS THE ONLY SLIGHTLY OPAQUE BIT OF THIS FUNCTION - SET EMPIRICALLY FROM EXPERIMENTATION.
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV > 0 & COMB_mats.tree.build$PVal > 0.01 & COMB_mats.tree.build$Genotype_bin != 1] <- 0.5 #If have any mutant reads, set as "?" as long as p-value > 0.01
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV >= 3 & COMB_mats.tree.build$PVal > 0.001 & COMB_mats.tree.build$Genotype_bin != 1] <- 0.5 #If have high numbers of mutant reads, should set as "?" even if incompatible p-value (may be biased sequencing)
  COMB_mats.tree.build$Genotype_bin[(COMB_mats.tree.build$NV == 0) & (COMB_mats.tree.build$PVal > 0.05)] <- 0.5 #Essentially if inadequate depth to exclude mutation, even if no variant reads
  
  Genotype_shared_bin = COMB_mats.tree.build$Genotype_bin[rowSums(COMB_mats.tree.build$Genotype_bin == 1) > 1,]
  
  #Save the input parameters to a list
  params = list(input_set_ID = input_set_ID,
                gender = gender,
                retain_muts = retain_muts,
                exclude_muts = exclude_muts,
                germline_pval = germline_pval,
                rho = rho,
                mean_depth = mean_depth,
                pval_dp2= pval_dp2,
                pval_dp3= pval_dp3,
                min_depth = min_depth,
                min_pval_for_true_somatic = min_pval_for_true_somatic,
                min_vaf = min_vaf,
                min_variant_reads_SHARED = min_variant_reads_SHARED,
                min_pval_for_true_somatic_SHARED = min_pval_for_true_somatic_SHARED,
                min_vaf_SHARED=min_vaf_SHARED)
  
  #Save the summary stats of the run
  summary = data.frame(total_mutations = sum(select_muts),
                       total_SNVs = sum(COMB_mats.tree.build$mat$Mut_type == "SNV"),
                       total_INDELs = sum(COMB_mats.tree.build$mat$Mut_type == "INDEL"),
                       shared_mutations = nrow(Genotype_shared_bin),
                       shared_SNVs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "SNV"),
                       shared_INDELs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "INDEL"))
  #Extract the dummy dna_strings for tree building with MPBoot
  dna_strings = dna_strings_from_genotype(Genotype_shared_bin)
  
  #Print the run stats to the screen
  cat(summary$total_mutations,"total mutations\n")
  cat(summary$total_SNVs,"total SNVs\n")
  cat(summary$total_INDELs,"total INDELs\n")
  cat(summary$shared_mutations,"shared mutations\n")
  cat(summary$shared_SNVs,"shared SNVs\n")
  cat(summary$shared_INDELs,"shared INDELs\n")
  
  output = list(COMB_mats.tree.build = COMB_mats.tree.build, Genotype_shared_bin= Genotype_shared_bin, filter_code=filter_code, params = params, summary = summary, dna_strings = dna_strings)
  return(output)
}


#Function to create the dummy dna strings from the shared genotype matrix
dna_strings_from_genotype = function(genotype_mat) {
  Muts = rownames(genotype_mat)
  Ref = rep("W", length(Muts)) #W = Wild-type
  Alt = rep("V", length(Muts)) #V = Variant
  
  dna_strings = list()
  dna_strings[1] = paste(Ref, sep = "", collapse = "")
  
  for (k in 1:ncol(genotype_mat)){
    Mutations = Ref
    Mutations[genotype_mat[,k]==0.5] = '?'
    Mutations[genotype_mat[,k]==1] = Alt[genotype_mat[,k]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[k+1]=dna_string
  }
  names(dna_strings)=c("Ancestral",colnames(genotype_mat))
  return(dna_strings)
}

#Function to create vcf files from the "mat" object
create_vcf_files = function(mat, select_vector = NULL) {
  if(is.null(select_vector)) {vcf_file = mat[,c("Chrom","Pos","Ref","Alt")]} else {vcf_file = mat[select_vector,c("Chrom","Pos","Ref","Alt")]}
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  vcf_file = vcf_file[,c(1,2,8,3,4,7,6,5)]
  return(vcf_file)
}

#Function to extract the nodes from the first few divisions of the tree (set divisions - default = 2)
get_early_nodes = function(tree,divisions=2) {
  root_no = length(tree$tip.label)+1
  nodes=root_no
  early_nodes = NULL
  for(i in 1:divisions) {
    daughter_nodes = tree$edge[tree$edge[,1] %in% nodes,2]
    early_nodes=c(early_nodes,daughter_nodes)
    nodes<-daughter_nodes
  }
  return(early_nodes)
}

#Function to calculate peak VAFs from sample mutations (after applying germline and beta-binomial filters) to screen for mixed colonies
check_peak_vaf = function(sample, COMB_mats, filter_params,rho_cutoff=0.3) {
  colnames(COMB_mats$NV) <- gsub(pattern = "_MTR", replacement = "",x = colnames(COMB_mats$NV))
  dens <- density((COMB_mats$NV/COMB_mats$NR)[COMB_mats$NV[,sample] >=2 & !COMB_mats$mat$Chrom %in% c("X","Y") & log10(filter_params$germline_pval) <(-10) & filter_params$bb_rhoval > rho_cutoff ,sample])
  return(dens$x[which.max(dens$y)])
}

#Function to visualize VAF plots in similar way to the above
vaf_density_plot = function(sample, COMB_mats, filter_params,rho_cutoff=0.3) {
  colnames(COMB_mats$NV)=colnames(COMB_mats$NR)=gsub(pattern = "_MTR", replacement = "",x = colnames(COMB_mats$NV))
  sample_mean_depth = mean(COMB_mats$NR[,sample])
  dens <- density((COMB_mats$NV/COMB_mats$NR)[COMB_mats$NV[,sample] >=2 & !COMB_mats$mat$Chrom %in% c("X","Y") & log10(filter_params$germline_pval) <(-10) & filter_params$bb_rhoval > rho_cutoff ,sample])
  plot(dens, main = sample,xlim=c(-0.05,1.05))
  abline(v = dens$x[which.max(dens$y)])
  text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=",round(dens$x[which.max(dens$y)], digits = 2),"\nMean coverage=",round(sample_mean_depth,digits =2)), col = "red", cex = 0.7)
}

vaf_density_plot_final=function(sample,tree,COMB_mats){
  node <- which(tree$tip.label==sample)
  sample_muts <- COMB_mats$mat$mut_ref[COMB_mats$mat$node %in% get_ancestral_nodes(node,tree$edge)]
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  dens <- density((COMB_mats$NV/COMB_mats$NR)[sample_muts,sample])
  plot(dens,xlim = c(0,1),main=sample)
  abline(v = dens$x[which.max(dens$y)])
  text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=",round(dens$x[which.max(dens$y)], digits = 2)),col="red",cex = 0.7)
}



#Functions to import the cgpVAF output matrices neatly into R in the format that is used in my scripts
import_cgpvaf_output=function(cgpvaf_output_file,ref_ID="PDv37is") {
  mat<-read.delim(cgpvaf_output_file,stringsAsFactors = FALSE)
  mat<- mat[,!grepl(ref_ID,colnames(mat))] #Remove the reference sample columns
  NV<- mat[,grepl("_MTR", colnames(mat))]; NR <- mat[,grepl("_DEP",colnames(mat))]
  colnames(NR)=colnames(NV)=gsub(pattern = "_MTR",replacement = "", colnames(NV))
  mat$mut_ref=paste(mat$Chrom,mat$Pos,mat$Ref,mat$Alt,sep = "-")
  mat<-mat[,c("Chrom","Pos","Ref","Alt","mut_ref")]
  rownames(NV)=rownames(NR)=mat$mut_ref
  combined_mats=list(mat,NV,NR); names(combined_mats) <-c("mat","NV","NR")
  return(combined_mats)
}

import_cgpvaf_SNV_and_INDEL = function(SNV_output_file,INDEL_output_file=NULL) {
  #Import cgpVAF snp output file for the single-cell colonies, create the mut_ref column, and extract the mut and dep cols
  SNV_mats=import_cgpvaf_output(SNV_output_file)
  SNV_mats$mat$Mut_type="SNV"
  if(!is.null(INDEL_output_file)) {
    INDEL_mats = import_cgpvaf_output(INDEL_output_file)
    INDEL_mats$mat$Mut_type="INDEL"
    
    #Only include samples that are included in both SNV and INDEL cgpVAF output
    samples_in_both=intersect(colnames(SNV_mats$NV),colnames(INDEL_mats$NV))
    print(paste(length(samples_in_both),"samples in both SNV and INDEL cgpVAF output matrices, and will be combined"))
    SNV_mats$NV<-SNV_mats$NV[,samples_in_both]
    SNV_mats$NR<-SNV_mats$NR[,samples_in_both]
    INDEL_mats$NV<-INDEL_mats$NV[,samples_in_both]
    INDEL_mats$NR<-INDEL_mats$NR[,samples_in_both]
    
    combined_mats=mapply(SNV_mats,INDEL_mats,FUN = rbind) #Bind the indel and snp matrices together
    return(combined_mats)
  } else {
    return(SNV_mats)
  }
}

#Function to split up the imported output from VAGRENT in a neat way
split_vagrent_output = function(df,split_col,col_IDs = c("Gene","Transcript","RNA","CDS","Protein","Type","SO_codes")) {
  col = df[[split_col]]
  output = matrix(nrow = nrow(df), ncol = length(col_IDs))
  for(i in 1:length(col_IDs)) {
    output[,i] = str_split(col, pattern = "\\|", simplify = TRUE)[,i]
  }
  colnames(output) = col_IDs
  output<-as.data.frame(output,stringsAsFactors=F)
  return(output)
}

#Function to check for mutations that have been called as germline that are in fact absent in a clade
#Run using: res=check_for_false_germline_calls(tree,COMB_mats = COMB_mats, filter_params=filter_params)
is.snv=function(mut_ref) {
  sub=stringr::str_split(mut_ref,pattern = "-",simplify=T)[3:4]
  res=ifelse(nchar(sub[1])==1 & nchar(sub[2])==1,T,F)
  return(res)
}

check_for_false_germline_calls = function(tree,
                                          COMB_mats,
                                          filter_params,
                                          max_clade_prop=0.1, #the cutoff size (proportion of total samples included in clade) to test the clade for absent germline mutations. At >10% the germline filter is unlikely to wrongly remove mutations.
                                          SNVs_only=T #Only re-add SNVs (indels are likely to be high frequency artefacts)
) {
  #Pull out the root clades
  get_root_clades=function(tree) {
    tree=di2multi(tree)
    ROOT=tree$edge[1,1]
    clades=tree$edge[tree$edge[,1]==ROOT,2]
    root_clade_samples=lapply(clades,function(node) return(getTips(tree,node)))
    return(root_clade_samples)
  }
  
  #Drop the ancestral tip if present, as this messes up the "get_root_clades" function
  tree.noancestral<-drop.tip(tree,"Ancestral")
  
  root_clades=get_root_clades(tree.noancestral)
  nsamp=length(tree.noancestral$tip.label)
  which_small=sapply(root_clades,length)<max_clade_prop*nsamp
  root_clades[!which_small]<-NULL
  
  if(length(root_clades)>0) {
    res<-lapply(root_clades,function(outlier_sample_group) {
      print(paste("Testing outlier group:",paste(outlier_sample_group,collapse=" ")))
      #Ensure all names are consistent
      rownames(COMB_mats$NR)=rownames(COMB_mats$NV)=rownames(COMB_mats$PVal)<-COMB_mats$mat$mut_ref
      colnames(COMB_mats$NR)=colnames(COMB_mats$PVal)=colnames(COMB_mats$NV)<-gsub("_MTR","",colnames(COMB_mats$NV))
      
      #Select mutations that were filtered by the germline filter
      germline_filtered=rownames(filter_params)[log10(filter_params$germline_pval)>(-10)]
      
      #Aggregate counts across an individual outlier sample/ sample group
      NR_outlier=apply(COMB_mats$NR[germline_filtered,outlier_sample_group,drop=F],1,sum)
      NV_outlier=apply(COMB_mats$NV[germline_filtered,outlier_sample_group,drop=F],1,sum)
      
      outlier_pvals=mapply(FUN=function(NV,NR) {if(NR==0){return(1)}else{binom.test(NV,NR,alternative="less")$p.value}},NV=NV_outlier,NR=NR_outlier)
      outlier_pval.adj=p.adjust(outlier_pvals,method = "BH")
      #hist(log10(outlier_pvals),breaks=50,main="Unadjusted p-values for mutations being present in outlier group") #Review the p-value histogram - any clear low outliers?
      
      #Test for germline filtered mutations that are likely to be absent (with Bon-Ferroni correction for multiple testing)
      mut_refs<-germline_filtered[outlier_pval.adj<0.05 & NV_outlier==0]
      if(SNVs_only & length(mut_refs)>0) {
        mut_refs<-mut_refs[sapply(mut_refs,is.snv)]
      }
      if(length(mut_refs)>0) {
        print(paste(mut_refs,"is convincingly absent in this group"))
        return(mut_refs)
      } else {
        print("There are no mutations called as germline that are robustly absent in this outlier group, though this would relies on adequate coverage")
        return(NULL)
      }
    })
    return(res)
  } else {
    print(paste0("There are no clades from the root that include <",max_clade_prop*100,"% of samples"))
    return(NULL)
  }
}

# #This is Nick's version of the function: adds the ancestral tip at the beginning (i.e. as position 1)
# add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
#   tmp=tree$edge
#   N=length(tree$tip.label)
#   ##Renumber what was root->max+1
#   ##Renumber the node with the sameid as new root as max+2
#   renamedroot=max(tmp+1)
#   tmp=ifelse(tmp==N+1,renamedroot,tmp)
#   ##tmp[which(tmp[,1]==(N+1)),1]=renamedroot
#   tmp=ifelse(tmp==N+2,renamedroot+1,tmp)
#   ##Increment tips by 1
#   tmp[,2]=ifelse(tmp[,2]<=N,tmp[,2]+1,tmp[,2])
#   
#   tree$edge=rbind(matrix(c(N+2,N+2,renamedroot,1),ncol=2,byrow  = FALSE),tmp)
#   tree$edge.length=c(0,0,tree$edge.length)
#   
#   tree$tip.label=c(outgroup_name,tree$tip.label)
#   tree$Nnode=tree$Nnode+1
#   mode(tree$Nnode)="integer"
#   mode(tree$edge)="integer"
#   tree
# }

#This version of the function adds the ancestral tip at the end
add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
  tmp=tree$edge
  N=length(tree$tip.label)
  newroot=N+2
  renamedroot=N+3
  ancestral_tip=N+1
  tmp=ifelse(tmp>N,tmp+2,tmp)
  
  tree$edge=rbind(c(newroot,renamedroot),tmp,c(newroot,ancestral_tip))
  tree$edge.length=c(0,tree$edge.length,0)
  
  tree$tip.label=c(tree$tip.label,outgroup_name)
  tree$Nnode=tree$Nnode+1
  mode(tree$Nnode)="integer"
  mode(tree$edge)="integer"
  tree
}

assign_mutations_to_branches=function(tree,filtered_muts,keep_ancestral=T,create_multi_tree=T,p.error.value=0.01,treefit_pval_cutoff=1e-3) {
  tree=drop.tip(tree,"Ancestral")
  if(!keep_ancestral) {
    print("Assigning mutation without an ancestral branch")
    tree <- multi2di(tree)
    tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
    
    #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
    df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
    
    #Get matrices in order, and run the main assignment functions
    mtr = filtered_muts$COMB_mats.tree.build$NV; mtr = as.matrix(mtr)
    depth = filtered_muts$COMB_mats.tree.build$NR; depth = as.matrix(depth)
    p.error = sapply(df$samples,function(x) ifelse(x=="Ancestral",1e-6,p.error.value))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
    
  } else {
    print("Assigning mutation with an ancestral branch")
    tree <- multi2di(tree)
    tree <- add_ancestral_outgroup(tree) #Re add the ancestral outgroup after making tree dichotomous - avoids the random way that baseline polytomy is resolved
    tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
    
    #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
    df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
    
    #Get matrices in order, and run the main assignment functions
    mtr = filtered_muts$COMB_mats.tree.build$NV; mtr$Ancestral=0; mtr = as.matrix(mtr)
    depth = filtered_muts$COMB_mats.tree.build$NR; depth$Ancestral=10; depth = as.matrix(depth)
    p.error = sapply(df$samples,function(x) ifelse(x=="Ancestral",1e-6,p.error.value))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  if(create_multi_tree){
    print("Converting to a multi-furcating tree structure")
    tree$edge.length <- res$df$df$edge_length #Assign edge lengths from the initial res object
    #Maintain the dichotomy with the ancestral branch
    if(keep_ancestral) {
      ROOT=tree$edge[1,1]
      current_length<-tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]
      new_length<-ifelse(current_length==0,1,current_length)
      tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]<-new_length
    }
    tree<-di2multi(tree) #Now make tree multifurcating
    df = reconstruct_genotype_summary(tree) #Define df (data frame) for new treeshape
    
    #Re-run the mutation assignment algorithm from the new tree
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  tree$edge.length <- res$df$df$edge_length #Assign edge lengths from the most recent res object
  res$tree<-tree #Add the tree to the res object
  
  #See how many mutations are "poor fit"
  poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
  print(paste(sum(poor_fit),"mutations do not have read counts that fit any tree branch well"))
  return(res)
}


##SPECIFIC SET OF FUNCTIONS DESIGNED FOR THE "Lesion_segregation_mutation_summaries.R" SCRIPT AND THE ANALYSIS

#Write a vcf file for reading into MutationalPatterns
write.vcf=function(details,vcf_path,select_vector=NULL,vcf_header_path="~/Documents/vcfHeader.txt") {
  if(any(class(details)=="character")) {
    mat=as.data.frame(stringr::str_split(details,pattern="-",simplify = T),stringsAsFactors=F)
    colnames(mat)<-c("Chrom","Pos","Ref","Alt")
    vcf=create_vcf_files(mat=mat,select_vector=select_vector)
  } else {
    vcf=create_vcf_files(mat=details,select_vector=select_vector)
  }
  write.table(vcf,sep = "\t", quote = FALSE,file=paste0(vcf_path,".temp"),row.names = F)
  system(paste0("cat ",vcf_header_path," ",vcf_path,".temp > ",vcf_path))
  system(paste0("rm ",vcf_path,".temp"))
}

#This function is required in the filtering function
get_ancestor_node=function(node,tree,degree=1){ #to get the 1st degree ancestor (i.e. the direct parent) use degree=1.  Use higher degrees to go back several generations.
  curr<-node
  for(i in 1:degree){
    curr=tree$edge[which(tree$edge[,2]==curr),1]
    if(curr==(1+length(tree$tip.label))) {stop(return(curr))}
  }
  return(curr)
}

#Function for Nick's alternative rho estimation approach (faster)
loglik=function(par,nmuts,depth){
  idx=which(!is.na(nmuts/depth))
  sum(VGAM::dbetabinom(x =nmuts[idx],size = depth[idx],prob = par[2],rho = par[1],log = T))
}

findrho=function(NV_vec,NR_vec){
  pseudo=1e-6
  optres=optim(par=c(0.1,0.1),loglik, gr = NULL,method="L-BFGS-B",lower=c(pseudo,pseudo),upper=c(0.8,0.9),control=list(fnscale=-1),nmuts=NV_vec,depth=NR_vec)
  #c(rho=optres$par[1],p=optres$par[2])
  optres$par[1]
}

#Define the "get_node_types" function required for following the lesion journey in the case of PVVs
#It uses the mutation dataframe ("mut_df") to work out whether daughter branches of a node are (a) a mutant allele (b) wild-type or (c) mixed
get_node_types=function(lesion_children,mut_df,tree) {
  types=sapply(lesion_children, function(node) {
    nodes=c(node,get_all_node_children(node,tree))
    if(all(!mut_df$neg_test[mut_df$clades%in%nodes])){
      return("pure_positive")
    } else if(all(!mut_df$pos_test[mut_df$clades%in%nodes])) {
      return("pure_negative")
    } else {
      return("mixed")
    }
  })
  return(types)
}

#The function to assess all the mutations for whether they are phylogeny breaking
create_PVV_filter_table=function(mutations_to_test,details,tree,matrices,look_back=3,type="standard",remove_duplicates=F,duplicate_samples=NULL,MC_CORES=1) {
  require(dplyr)
  require(parallel)
  print(paste("Assessing",length(mutations_to_test),"mutations for whether they are truly phylogeny breaking"))
  if(look_back=="all"){
    all_clades=unique(tree$edge[,2])
    all_clade_nodes=lapply(all_clades,function(node) c(node,get_all_node_children(node,tree)))
  }
  filter_output=mclapply(mutations_to_test,function(mut) {
    if(which(mutations_to_test==mut)%%10==0){print(which(mutations_to_test==mut))}
    allocated_node=details$node[details$mut_ref==mut]
    
    #Get counts of all individual clades within the allocated node i.e. those expected to be positive
    if(allocated_node<=length(tree$tip.label)) {
      positive_clades=allocated_node
    } else {
      positive_clades=get_all_node_children(allocated_node,tree)
    }
    
    if(remove_duplicates) {positive_clades=positive_clades[!positive_clades%in%duplicate_samples]} #don't assess duplicate samples as individual samples
    positive_clade_nodes=lapply(positive_clades,function(node) c(node,get_all_node_children(node,tree)))
    positive_clade_samples=lapply(positive_clade_nodes,function(nodes) return(tree$tip.label[nodes[nodes%in%1:length(tree$tip.label)]]))
    NV_pos=unlist(lapply(positive_clade_samples,function(samples) sum(matrices$NV[mut,samples])))
    NR_pos=unlist(lapply(positive_clade_samples,function(samples) sum(matrices$NR[mut,samples])))
    
    #If we previously removed duplicates, need to rederive the "positive_clade_nodes" including duplicates for filtering out ALL positive clades in the next section
    if(remove_duplicates) { 
      positive_clades=get_all_node_children(allocated_node,tree)
      positive_clade_nodes=lapply(positive_clades,function(node) c(node,get_all_node_children(node,tree)))
    }
    
    #Get counts of nearby individual clades that don't include the allocated node, i.e. those expected to be negative
    if(look_back=="all"){
      select=lapply(all_clade_nodes,function(nodes)if(!any(nodes%in%positive_clade_nodes)){TRUE}else{FALSE})
      negative_clades=all_clades[unlist(select)]
      if(remove_duplicates) {negative_clades=negative_clades[!negative_clades%in%duplicate_samples]}
      negative_clade_nodes=lapply(negative_clades,function(node) c(node,get_all_node_children(node,tree)))
      #negative_clade_nodes<-negative_clade_nodes[-which(sapply(negative_clade_nodes,function(x) all(x==which(tree$tip.label=="Ancestral"))))]
      negative_clade_samples=lapply(negative_clade_nodes,function(nodes) return(tree$tip.label[nodes[nodes%in%1:length(tree$tip.label)]]))
      NV_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NV[mut,samples])))
      NR_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NR[mut,samples])))
    } else {
      #Assess only "nearby" to the allocated node (up to the degree specified by look_back argument)
      ancestral_node=get_ancestor_node(allocated_node,tree,degree=look_back)
      all_clades=c(ancestral_node,get_all_node_children(ancestral_node,tree))
      all_clade_nodes=lapply(all_clades,function(node) c(node,get_all_node_children(node,tree)))
      
      select=lapply(all_clade_nodes,function(nodes)if(!any(nodes%in%positive_clade_nodes)){TRUE}else{FALSE})
      negative_clades=all_clades[unlist(select)]
      if(remove_duplicates) {negative_clades=negative_clades[!negative_clades%in%duplicate_samples]}
      negative_clade_nodes=lapply(negative_clades,function(node) c(node,get_all_node_children(node,tree)))
      negative_clade_samples=lapply(negative_clade_nodes,function(nodes) return(tree$tip.label[nodes[nodes%in%1:length(tree$tip.label)]]))
      NV_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NV[mut,samples])))
      NR_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NR[mut,samples])))
    }
    
    #Apply beta-binomial filter to these counts to get sense of overdispersion
    NR_pos[NR_pos==0]<-1; NR_neg[NR_neg==0]<-1 #Set sites with 0 depth to a depth of 1
    pos_rho=estimateRho_gridml(NV_vec=NV_pos,NR_vec=NR_pos)
    neg_rho=estimateRho_gridml(NV_vec=NV_neg,NR_vec=NR_neg)
    
    #Additional filters to check not just for overdispersion, but that at least one clade unexpectedly truly negative or positive
    pos_test=any(NV_pos==0 & NR_pos>=10) #is there at least one "WT" sub-clade with zero variant reads with a depth of ≥ 10
    
    if(type=="complete") {
      neg_test=any((NV_neg/NR_neg)>=0.2 & NR_neg >=6) #is there at least one anticipated "negative" clade, that in fact has a VAF>0.3 with a depth ≥8
    } else {
      neg_test=any((NV_neg/NR_neg)>=0.3 & NR_neg >=8) #is there at least one anticipated "negative" clade, that in fact has a VAF>0.3 with a depth ≥8
    }
    
    
    mut_params=data.frame(mut=mut,
                          node=allocated_node,
                          pos_rho=pos_rho,
                          neg_rho=neg_rho,
                          pos_test=pos_test,
                          neg_test=neg_test,
                          pval=details$pval[which(details$mut_ref==mut)])
    return(mut_params) 
  },mc.cores = MC_CORES)
  
  filter_output_df=dplyr::bind_rows(filter_output)
  return(filter_output_df)
}

#As above, but can incorporate two alternative mutant alleles, therefore suitable for MAVs
get_MAV_node_types=function(lesion_children,mut_df,tree) {
  types=sapply(lesion_children, function(node) {
    nodes=c(node,get_all_node_children(node,tree))
    if(all(!mut_df$neg_test[mut_df$clades%in%nodes])&any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])&!any(mut_df$mut2_pos_test[mut_df$clades%in%nodes])){
      return("pure_mut1")
    } else if((sum(!mut_df$neg_test[mut_df$clades%in%nodes])/length(!mut_df$neg_test[mut_df$clades%in%nodes])>0.98)&any(mut_df$mut2_pos_test[mut_df$clades%in%nodes])&!any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])){
      return("pure_mut2")
    }else if(all(!mut_df$mut1_pos_test[mut_df$clades%in%nodes]) & all(!mut_df$mut2_pos_test[mut_df$clades%in%nodes])) {
      return("pure_negative")
    } else {
      return("mixed")
    }
  })
  return(types)
}

#Function to reclassify nearby mutations that are on the same branch as multi-nucleotide variants
#Note - this function ASSUMES that the mutations are phased, and this should be checked
#Needs a reference genome file
reclassify_MNVs=function(COMB_mats,region_size=2,genomeFile) {
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  library("stringr")
  
  details=COMB_mats$mat
  NV=COMB_mats$NV
  NR=COMB_mats$NR
  
  Chroms=c(1:22,"X","Y")
  out_list_by_chrom=lapply(Chroms,function(Chrom) {
    details_by_chrom=details[details$Chrom==as.character(Chrom) & details$Mut_type=="SNV",]
    out_list=lapply(1:nrow(details_by_chrom),function(i) {
      Pos=as.numeric(details_by_chrom$Pos[i])
      node=details_by_chrom$node[i]
      near_muts=which(as.numeric(details_by_chrom$Pos)>Pos&
                        as.numeric(details_by_chrom$Pos)<=(Pos+region_size)&
                        details_by_chrom$node==node)
      near_muts<-near_muts[near_muts!=i]
      if(length(near_muts)==0) {
        return(NA)
      } else {
        return(c(details_by_chrom$mut_ref[i],details_by_chrom$mut_ref[near_muts]))
      }
    })
    out_list[sapply(out_list,function(x) is.na(x[1]))]<-NULL
    return(out_list)
  })
  
  out_list_by_chrom=lapply(out_list_by_chrom, function(list) {
    i=1
    while(i<length(list)){
      muts<-list[[i]]
      if(any(list[[i+1]]%in%muts)) {list[[i]]<-unique(c(list[[i]],list[[i+1]]));list[[i+1]]<-NULL} else {i<-i+1}
    }
    return(list)
  })
  
  #Combine the chromosomes into one list
  out_list=unlist(out_list_by_chrom,recursive=F)
  
  print(paste("There are",length(out_list),"SNV pairs that will be reclassified as MNVs"))
  
  # #Trinucleotide substitutions will occur more than once - therefore need to replace these (not yet done)
  if(length(out_list)>0) {
    new_df=lapply(out_list,function(x) {
      chrom=str_split(x[1],pattern="-",simplify=T)[,1]
      pos=as.numeric(str_split(x,pattern="-",simplify=T)[,2])
      ref=str_split(x,pattern="-",simplify=T)[,3]
      alt=str_split(x,pattern="-",simplify=T)[,4]
      node=details$node[details$mut_ref==x[1]]
      pval=details$pval[details$mut_ref%in%x]
      if(length(pos)==2 & pos[2]==(pos[1]+1)){
        new_ref=paste(ref,collapse="")
        new_alt=paste(alt,collapse="")
        return(data.frame(mut_ref=paste(chrom,pos[1],new_ref,new_alt,sep="-"),Chrom=chrom,Pos=pos[1],Ref=new_ref,Alt=new_alt,Mut_type="MNV",node=node,pval=mean(pval)))
      } else {
        new_ref=as.vector(scanFa(genomeFile, GRanges(chrom, IRanges(min(pos),max(pos)))))
        new_ref_vec=as.vector(str_split(new_ref,"",simplify=T))
        names(new_ref_vec)=seq(min(pos),max(pos),by=1)
        new_alt_vec=new_ref_vec
        for(i in 1:length(pos)){
          new_alt_vec[as.character(pos[i])]<-alt[i]
        }
        new_alt=paste(new_alt_vec,collapse="")
        return(data.frame(mut_ref=paste(chrom,pos[1],new_ref,new_alt,sep="-"),Chrom=chrom,Pos=pos[1],Ref=new_ref,Alt=new_alt,Mut_type="MNV",node=node,pval=mean(pval)))
      }
    })
    new_df=dplyr::bind_rows(new_df)
    
    replaced_SNVs=unlist(out_list)
    SNVs_for_counts=sapply(out_list,function(x) x[1])
    
    new_NV=NV[SNVs_for_counts,]
    new_NR=NR[SNVs_for_counts,]
    
    rownames(new_NV)=rownames(new_NR)<-new_df$mut_ref
    
    details_new<-details[!details$mut_ref%in%replaced_SNVs,]
    details_new<-rbind(details_new[,colnames(new_df)],new_df)
    
    return(list(mat=details_new,NV=rbind(NV[!details$mut_ref%in%replaced_SNVs,],new_NV),NR=rbind(NR[!details$mut_ref%in%replaced_SNVs,],new_NR)))
  } else {
    return(COMB_mats)
  }
}

#Function to find potential multi-allelic variants that may be caused by persistent DNA lesions
#Assesses for overlapping positions of reference
get_multi_allelic_variant_list=function(details,SNV_only=F) {
  Chroms=c(1:22,"X","Y")
  details$Ref=as.character(details$Ref)
  details$Alt=as.character(details$Alt)
  if(SNV_only){
    if(!"Chrom_pos"%in%colnames(details)){
      details$Chrom_pos<-paste(details$Chrom,details$Pos,sep = "-")
    }
    duplicates=details$Chrom_pos[duplicated(details$Chrom_pos)]
    dups_list<-lapply(duplicates,function(Chrom_pos) {return(details$mut_ref[details$Chrom_pos==Chrom_pos])})
    return(dups_list)
  } else {
    out_list_by_chrom=lapply(Chroms,function(Chrom) {
      if(sum(details$Chrom==Chrom)>1) {
        print(paste("Analysing chromosome",Chrom))
        #Split comparison of mutations by chromosome
        details_by_chrom=details[details$Chrom==as.character(Chrom),]
        details_by_chrom=details_by_chrom[order(details_by_chrom$Pos),]
        dups_list=lapply(1:(nrow(details_by_chrom)-1),function(i) {
          if(i%%1000==0) {print(i)}
          if(nchar(details_by_chrom$Ref[i])==1) {
            Pos=as.numeric(details_by_chrom$Pos[i])
          } else {
            Pos=as.numeric(details_by_chrom$Pos[i]):(as.numeric(details_by_chrom$Pos[i])+nchar(details_by_chrom$Ref[i])-1)
          }
          out=sapply((i+1):min(i+5,nrow(details_by_chrom)),function(j) {
            if(nchar(details_by_chrom$Ref[j])==1) {
              pos=as.numeric(details_by_chrom$Pos[j])
            } else {
              pos=as.numeric(details_by_chrom$Pos[j]):(as.numeric(details_by_chrom$Pos[j])+nchar(details_by_chrom$Ref[j])-1)
            }
            if(length(intersect(pos,Pos))>0) {
              return(T)
            } else {
              return(F)
            }
          })
          if(any(out)) {return(c(details_by_chrom$mut_ref[i],details_by_chrom$mut_ref[(i+1):min(i+5,nrow(details_by_chrom))][out]))} else {return(NA)}
        })
        dups_list[sapply(dups_list,function(x) is.na(x[1]))]<-NULL
        return(dups_list)
      } else {
        return(NULL)
      }
      
    })
    out_list=unlist(out_list_by_chrom,recursive=F)
    return(out_list)
  }
}

#Find the latest possible timing of the acquisition of the lesion
find_PVV_lesion_node=function(mut,allocated_node,pos_test,neg_test,tree,matrices) {
  
  #Define the get_ancestor_node function
  get_ancestor_node=function(node,tree,degree=1){ #to get the 1st degree ancestor (i.e. the direct parent) use degree=1.  Use higher degrees to go back several generations.
    curr<-node
    for(i in 1:degree){
      curr=tree$edge[which(tree$edge[,2]==curr),1]
      if(curr==(1+length(tree$tip.label))) {stop(return(curr))}
    }
    return(curr)
  }
  
  mut_df=create_mut_df(mut=mut,tree=tree,matrices=matrices)
  
  if(pos_test) {
    #if the negative sub-clade is within the allocated node, then "pos_test" will be true and "allocated_node" is the "initial_lesion_node"
    initial_lesion_node<-allocated_node
  } else if (neg_test){
    #If there is a positive clade outside the allocated node, then "neg_test" will be true.  In this case need to find the positive clade.
    #Do this be iteratively going from the allocated node to its ancestral node and looking for the ancestral node that contains ALL positive clades in the tree
    all_clades=unique(tree$edge[,2])
    all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))
    
    #Iteratively look back through ancestral nodes to find the one encasing other positive samples
    j=1
    repeat {
      ancestor=get_ancestor_node(allocated_node,tree,degree=j)
      ancestor_tips=getTips(tree,ancestor)
      #Look at clades that don't have any of the samples in "ancestor_tips". If this ancestor is the "initial_lesion_node", none will meet the "pos_test".
      if(!any(mut_df$pos_test[unlist(lapply(all_clade_samples, function(samples) !any(samples %in% ancestor_tips)))])) {
        break #Once this criteria is met, do not need to look back any further
      }
      else if(ancestor==tree$edge[1,1]){
        break
      }
      j=j+1
    }
    initial_lesion_node<-ancestor #The initial lesion node is therefore the most recent "ancestor" from the previous loop
  }
  return(initial_lesion_node)
}

find_MAV_lesion_node=function(node1,node2,tree,Chrom="auto") {
  if(node1==node2 & !Chrom%in%c("X","Y")) {
    stop(return(return(list(initial_lesion_node=NA,Filter="FAIL",Class="FAIL"))))
  } else {
    #get the ancestral nodes
    ancestor_1=tree$edge[tree$edge[,2]==node1,1]
    ancestor_2=tree$edge[tree$edge[,2]==node2,1]
    ancestors=c(ancestor_1,ancestor_2)
    ancestor_heights=sapply(ancestors,function(node) nodeheight(tree,node))
    
    if(ancestor_1==ancestor_2){
      Filter="PASS"
      Class="simple"
      initial_lesion_node=get_ancestor_node(node = node1,tree = tree)
    } else if(ancestor_1%in%get_all_node_children(ancestor_2,tree)|ancestor_2%in%get_all_node_children(ancestor_1,tree)){
      Filter="PASS"
      Class="removed"
      if(which.min(ancestor_heights)==1 & node2%in%get_all_node_children(node1,tree)) {
        initial_lesion_node<-node1
      } else if(which.min(ancestor_heights)==2 & node1%in%get_all_node_children(node2,tree)){
        initial_lesion_node<-node2
      } else {
        initial_lesion_node<-ancestors[which.min(ancestor_heights)]
      }
    } else {
      return(return(list(initial_lesion_node=NA,Filter="FAIL",Class="FAIL")))
    }
  }
  return(list(initial_lesion_node=initial_lesion_node,Filter=Filter,Class=Class))
}

#This function extracts a straight-forward phasing info df from the julia algorithm output
extract_phasing_info=function(list,Ref,Alt) {
  phasing_df=dplyr::bind_rows(list)
  #basects_df=Reduce(rbind,lapply(list,function(list) return(list[[2]])))
  if(is.logical(phasing_df)) {
    stop(return(NA))
  } else if(nrow(phasing_df)==0){
    stop(return(NA))
  }
  phasing_df$mut_base=sapply(strsplit(phasing_df$Mutation_allele,split = "="),function(x) x[2])
  phasing_df$snp_base=sapply(strsplit(phasing_df$SNP_allele,split = "="),function(x) x[2])
  
  SNP_sites=sort(unique(str_split(phasing_df$SNP_allele,pattern = "=",simplify = T)[,1]))
  
  phasing_by_SNP_list=lapply(SNP_sites,function(SNP_site) {
    phasing_df_snp<-phasing_df[grepl(SNP_site,phasing_df$SNP_allele),]
    
    #Summarise mut allele phasing
    if(any(phasing_df_snp$mut_base%in%Alt)) {
      alt_phasing=table(phasing_df_snp$snp_base[phasing_df_snp$mut_base%in%Alt])
      alt_phases_with_base=names(alt_phasing)[which.max(alt_phasing)]
      n_alt_reads_supporting=alt_phasing[alt_phases_with_base]
      n_alt_reads_against=(sum(alt_phasing)-n_alt_reads_supporting)
    } else {
      alt_phasing=NA
      alt_phases_with_base=NA
      n_alt_reads_supporting=0
      n_alt_reads_against=0
    }
    
    #Summarise wt allele phasing
    if(any(phasing_df_snp$mut_base==Ref)){
      ref_phasing=table(phasing_df_snp$snp_base[phasing_df_snp$mut_base==Ref])
      ref_phases_with_base=names(ref_phasing)[which.max(ref_phasing)]
      n_ref_reads_supporting=ref_phasing[ref_phases_with_base]
      n_ref_reads_against=(sum(ref_phasing)-n_ref_reads_supporting)
    } else {
      ref_phasing=NA
      ref_phases_with_base=NA
      n_ref_reads_supporting=0
      n_ref_reads_against=0
    }
    
    return(data.frame(SNP_site=SNP_site,
                      alt_phases_with_base=alt_phases_with_base,
                      n_alt_reads_supporting=n_alt_reads_supporting,
                      n_alt_reads_against=n_alt_reads_against,
                      ref_phases_with_base=ref_phases_with_base,
                      n_ref_reads_supporting=n_ref_reads_supporting,
                      n_ref_reads_against=n_ref_reads_against
    )
    )
  })
  phasing_by_SNP_df=dplyr::bind_rows(phasing_by_SNP_list)
  return(phasing_by_SNP_df)
}

#Function will look in the supplied output_dir to see if phasing output for given sample/Chrom/Pos already exists, if not will run the .jl script. Imports the data.
#Run example: get_phasing_list(samples=positive_samples1,Chrom=Chrom,Pos=Pos,project=project,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set)

get_phasing_list=function(samples,
                          Chrom,
                          Pos,
                          project,
                          tree=NULL,
                          output_dir,
                          ref_sample_set,
                          distance=1000,
                          force_rerun=F,
                          verbose=F,
                          use_tree=T,
                          phasing_dir="/lustre/scratch126/casm/team154pc/ms56/my_programs/Mike_phasing/", #directory with the phasing scripts
                          julia_path="/lustre/scratch126/casm/team154pc/ms56/my_programs/julia-1.1.0/bin/julia") {
  wd<-getwd()
  setwd(phasing_dir) #Need to be in this directory for the function
  if(is.numeric(project)) {
    phasing_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample));print(paste("Reference sample set chosen as",ref_sample_set))}
      if(!file.exists(phasing_output_file)|file.info(phasing_output_file)$size==0|force_rerun) {
        #This section is to account for the long bam headers in sample PD44579b which interfere with the script
        if(grepl("PD44579b",sample)) {
          #Import all the necessary bams with edited headers
          sapply(c(sample,unlist(strsplit(ref_sample_set,","))),function(bam_sample) {
            new_bam_path=paste0("new_bams/",bam_sample,".sample.dupmarked.bam")
            if(!file.exists(new_bam_path)) {
              print(paste("Importing bam file for",bam_sample,"and replacing header"))
              bam_path=paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",bam_sample,"/",bam_sample,".sample.dupmarked.bam")
              command=paste(julia_path,bam_path,"offending_string.txt")
              system(command)
            }
          })
          #Now run using the modified julia script to use these local files
          command=paste(julia_path,"DRIVER_phasing_specify_BAM_directory.jl",Chrom,Pos,sample,paste0(phasing_dir,"new_bams"),as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        } else {
          command=paste(julia_path,"DRIVER_phasing.jl",Chrom,Pos,sample,project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        }
      } else if(verbose) {
        print("Existing phasing files found in specified output directory")
      }
      if(file.exists(phasing_output_file)&file.info(phasing_output_file)$size!=0){
        phasing=read.table(phasing_output_file,header = T,stringsAsFactors = F,colClasses="character")
      } else {
        print("Unable to run phasing script")
        phasing=NA
      }
      return(phasing)
    })
  } else if(is.data.frame(project)) {
    phasing_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample))}
      sample_project=project$project[project$sample==sample]
      if(use_tree){
        set.seed(1)
        ref_sample_set=paste0(sample(x=tree$tip.label[tree$tip.label%in%project$sample[project$project==sample_project]],size=5),collapse=",")
      } else {
        sample_stem=stringr::str_split(sample,pattern = "_",simplify=T)[,1]
        set.seed(1)
        ref_sample_set=paste0(sample(x=project$sample[project$project==sample_project & grepl(sample_stem,project$sample)],size=5),collapse=",")
      }
      if(verbose) {print(paste("Ref sample set chosen as",ref_sample_set))}
      
      if(!file.exists(phasing_output_file)|file.info(phasing_output_file)$size==0|force_rerun) {
        command=paste(julia_path,"DRIVER_phasing.jl",Chrom,Pos,sample,sample_project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
        system(command)
      } else if(verbose) {
        print("Existing phasing files found in specified output directory")
      }
      if(file.exists(phasing_output_file)&file.info(phasing_output_file)$size!=0){
        phasing=read.table(phasing_output_file,header = T,stringsAsFactors = F,colClasses="character")
      } else {
        phasing="Unable to run phasing script"
      }
      return(phasing)
    })
  }
  setwd(wd)
  return(phasing_list)
}


get_base_counts_list=function(samples,
                              Chrom,
                              Pos,
                              project,
                              tree=NULL,
                              output_dir,
                              ref_sample_set,
                              distance=1000,
                              force_rerun=F,
                              verbose=F,
                              use_tree=T,
                              julia_path="/lustre/scratch126/casm/team154pc/ms56/my_programs/julia-1.1.0/bin/julia") {
  wd<-getwd()
  setwd("/lustre/scratch126/casm/team154/ms56/my_programs/Mike_phasing") #Need to be in this directory for the function
  if(is.numeric(project)) {
    basects_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample));print(paste("Reference sample set chosen as",ref_sample_set))}
      if(!file.exists(basects_output_file)|file.info(basects_output_file)$size==0|force_rerun) {
        if(grepl("PD44579b",sample)) {
          #Import all the necessary bams with edited headers
          sapply(c(sample,unlist(strsplit(ref_sample_set,","))),function(bam_sample) {
            new_bam_path=paste0("new_bams/",bam_sample,".sample.dupmarked.bam")
            if(!file.exists(new_bam_path)) {
              print(paste("Importing bam file for",bam_sample,"and replacing header"))
              bam_path=paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",bam_sample,"/",bam_sample,".sample.dupmarked.bam")
              command=paste(julia_path,"header_edit.jl",bam_path,"offending_string.txt")
              system(command)
            }
          })
          #Now run using the modified julia script to use these local files
          command=paste(julia_path,"DRIVER_phasing_specify_BAM_directory.jl",Chrom,Pos,sample,"/lustre/scratch126/casm/team154pc/ms56/my_programs/Mike_phasing/new_bams",as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        } else {
          command=paste(julia_path,"DRIVER_phasing.jl",Chrom,Pos,sample,project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        }
      } else if(verbose) {
        print("Existing base counts files found in specified output directory")
      }
      if(file.exists(basects_output_file)&file.info(basects_output_file)$size!=0){
        basects=read.table(basects_output_file,header = T,stringsAsFactors = F)
      } else {
        basects="Unable to run phasing script"
      }
      return(basects)
    })
  } else if(is.data.frame(project)) {
    basects_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample))}
      sample_project=project$project[project$sample==sample]
      if(use_tree){
        set.seed(1)
        ref_sample_set=paste0(sample(x=tree$tip.label[tree$tip.label%in%project$sample[project$project==sample_project]],size=5),collapse=",")
      } else {
        sample_stem=stringr::str_split(sample,pattern = "_",simplify=T)[,1]
        set.seed(1)
        ref_sample_set=paste0(sample(x=project$sample[project$project==sample_project & grepl(sample_stem,project$sample)],size=5),collapse=",")
      }
      if(verbose) {print(paste("Ref sample set chosen as",ref_sample_set))}
      if(!file.exists(basects_output_file)|file.info(basects_output_file)$size==0|force_rerun) {
        command=paste(julia_path,"DRIVER_phasing.jl",Chrom,Pos,sample,sample_project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
        system(command)
      } else if(verbose) {
        print("Existing base counts files found in specified output directory")
      }
      if(file.exists(basects_output_file)&file.info(basects_output_file)$size!=0){
        basects=read.table(basects_output_file,header = T,stringsAsFactors = F)
      } else {
        basects="Unable to run phasing script"
      }
      return(basects)
    })
  }
  setwd(wd)
  return(basects_list)
}

get_clade_base_counts=function(nodes,tree,Chrom,Pos,project,ref_sample_set,phasing_output_dir,distance=1000,force_rerun=F) {
  clade_base_counts_list=lapply(nodes,function(node) {samples=getTips(tree=tree,node=node);base_counts_list=get_base_counts_list(samples=samples,Chrom=Chrom,Pos=Pos,project=project,tree=tree,output_dir = phasing_output_dir,ref_sample_set = ref_sample_set,distance=distance,force_rerun = force_rerun);return(base_counts_list)})
  aggregated_clade_base_counts_list=lapply(clade_base_counts_list,function(list) {
    chrom_pos_df=list[[1]][,c(2,3)] #Get the co-ordinates of apparent het SNPs from the 1st in the list
    list_mod=lapply(list,function(df) {res<-left_join(chrom_pos_df,df[,-1],by=c("Chr","Pos"));res[is.na(res)]<-0;return(res[,3:7])}) #Now get just the base counts at these sites
    counts_df=Reduce(function(df1,df2) {return(df1+df2)},list_mod) #Aggregate these across a pure subclade
    return(cbind(chrom_pos_df,counts_df))
  })
  return(aggregated_clade_base_counts_list)
}

return_heterozygous_SNPs=function(base_counts_list) {
  pos_het_res=lapply(base_counts_list,function(df) {
    if(nrow(df)==0) {stop(return(NA))}
    het_SNPs=apply(df[,3:7],1,function(x) {
      counts=x; sum_counts=sum(x)
      if(sum_counts==0) {stop(return(NA))}
      het_test=sapply(counts,function(y) return(binom.test(x=y,n=sum_counts)$p.value))
      absent_test=sapply(counts,function(y) return(binom.test(x=y,n=sum_counts,p=0.01)$p.value))
      hom_test=sapply(counts,function(y) return(binom.test(x=y,n=sum_counts,p=0.99)$p.value))
      lik_hom=prod(pmax(hom_test,absent_test))
      lik_het=prod(pmax(het_test,absent_test))
      if(lik_het>lik_hom) {return(T)} else {return(F)}
    })
    return(het_SNPs)
  })
  positions=base_counts_list[[1]]$Pos
  het_positions=positions[Reduce(function(x,y) {x&y},pos_het_res)]
  return(het_positions)
}

check_matching_phasing=function(phasing_info1,phasing_info2,het_positions=NULL) {
  if(any(sapply(list(phasing_info1,phasing_info2),is.logical))) {
    result<-"Unable to confirm phasing"
  } else {
    comb_df=inner_join(phasing_info1,phasing_info2,by="SNP_site")
    comb_df$depth=comb_df$n_alt_reads_supporting.x+comb_df$n_alt_reads_supporting.y+comb_df$n_ref_reads_supporting.x+comb_df$n_ref_reads_supporting.y
    
    #Test if the called SNPs appear real, or if previously tested on the basects data, filter the included SNPs based on this
    if(is.null(het_positions)) {
      true_het1<-comb_df$alt_phases_with_base.x!=comb_df$ref_phases_with_base.x
      true_het2<-comb_df$alt_phases_with_base.y!=comb_df$ref_phases_with_base.y
      het_test=mapply(FUN=function(x,y) {vec=c(x,y);vec<-vec[!is.na(vec)];if(length(vec)==0) {return(T)} else if(all(vec)) {return(T)} else {return(F)}},x=true_het1,y=true_het2)
      
      #If no true het SNPs, stop function; else filter the comb_df
      if(!any(het_test)) {
        stop(return("No SNPs appear to be heterozygous"))
      } else {
        comb_df<-comb_df[het_test,]
      }
      
      #Although heterozygosity is likely after the above test, it is not confirmed. Test for this:
      het_confirmed=apply(comb_df[,c("ref_phases_with_base.x","ref_phases_with_base.y","alt_phases_with_base.x","alt_phases_with_base.y")],1,function(x) length(unique(x[!is.na(x)]))>1)
      if(any(het_confirmed)) {
        comb_df<-comb_df[het_confirmed,]
        het_not_confirmed<-F
      } else {
        het_not_confirmed<-T
      }
      
    } else {
      snp_positions=as.numeric(str_split(pattern=":",comb_df$SNP_site,simplify=T)[,2])
      comb_df<-comb_df[snp_positions%in%het_positions,]
      if(nrow(comb_df)==0) {
        stop(return("No SNPs appear to be heterozygous"))
      }
      het_not_confirmed<-F
    }
    
    
    #Test for either alt matching alt, or ref matching ref for any individual SNP
    matching_res=list(Matching_alt_phasing=comb_df$alt_phases_with_base.x==comb_df$alt_phases_with_base.y,
                      Matching_ref_phasing=comb_df$ref_phases_with_base.x==comb_df$ref_phases_with_base.y,
                      Non_matching_alt_ref_phasing=comb_df$alt_phases_with_base.x!=comb_df$ref_phases_with_base.y,
                      Non_matching_ref_alt_phasing=comb_df$ref_phases_with_base.x!=comb_df$alt_phases_with_base.y)
    
    matching_res=lapply(matching_res, function(vec) {
      names(vec)<-1:length(vec)
      vec_no_NAs<-vec[!is.na(vec)]
      if(length(unique(vec_no_NAs))>1) { #If there is disagreement between different SNPs, retain the highest depth ones only
        vec_no_NAs<-vec_no_NAs[as.character(which(comb_df$depth>median(comb_df$depth)))]
        vec_no_NAs<-vec_no_NAs[!is.na(vec_no_NAs)]
      }
      return(vec_no_NAs)
    })
    
    if(all(sapply(matching_res,function(x) length(x)==0))){
      result<-"Unable to confirm phasing"
    } else if(all(unlist(matching_res))) {
      if(het_not_confirmed) {
        result<-"Same phasing suggested, though SNP heterozygosity not confirmed"
      } else {
        result<-"Same phasing confirmed"
      }
    } else {
      if(het_not_confirmed) {
        result<-"Non-matching suggested, though SNP heterozygosity not confirmed"
      } else {
        result<-"Non-matching phasing confirmed"
      }
    }
  }
  return(result)
}


check_matching_phasing_non_clonal=function(phasing_info1,phasing_info2,het_positions=NULL) {
  if(any(sapply(list(phasing_info1,phasing_info2),is.logical))) {
    result<-"Unable to confirm phasing"
  } else {
    if(is.null(het_positions)){stop(return("Must supply list of confirmed heterozygous SNP positions"))}
    comb_df=inner_join(phasing_info1,phasing_info2,by="SNP_site")
    comb_df$depth=comb_df$n_alt_reads_supporting.x+comb_df$n_alt_reads_supporting.y+comb_df$n_ref_reads_supporting.x+comb_df$n_ref_reads_supporting.y
    snp_positions=as.numeric(str_split(pattern=":",comb_df$SNP_site,simplify=T)[,2])
    comb_df<-comb_df[snp_positions%in%het_positions,]
    if(nrow(comb_df)==0) {
      stop(return("No SNPs appear to be heterozygous"))
    }
    
    #Test for either alt matching alt, or ref matching ref for any individual SNP
    matching_res=list(Matching_alt_phasing=comb_df$alt_phases_with_base.x==comb_df$alt_phases_with_base.y)
    
    matching_res=lapply(matching_res, function(vec) {
      names(vec)<-1:length(vec)
      vec_no_NAs<-vec[!is.na(vec)]
      if(length(unique(vec_no_NAs))>1) { #If there is disagreement between different SNPs, retain the highest depth ones only
        vec_no_NAs<-vec_no_NAs[as.character(which(comb_df$depth>median(comb_df$depth)))]
      }
      return(vec_no_NAs)
    })
    
    if(all(sapply(matching_res,function(x) length(x)==0))){
      result<-"Unable to confirm phasing"
    } else if(all(unlist(matching_res))) {
      result<-"Same phasing confirmed"
    } else {
      result<-"Non-matching phasing confirmed"
    }
  }
  return(result)
}

assess_phasing_non_clonal=function(phasing_summaries,sample_sets,Chrom,Pos,project,tree=NULL,output_dir,ref_sample_set,distance=1000,use_tree=T){
  phase_sum1=phasing_summaries[[1]]
  phase_sum2=phasing_summaries[[2]]
  if(is.logical(phase_sum1)|is.logical(phase_sum2)) {
    stop(return("Unable to phase"))
  } else {
    clade_base_counts_list=lapply(sample_sets,function(samples) {base_counts_list=get_base_counts_list(samples=samples,Chrom=Chrom,Pos=Pos,project=project,tree=tree,output_dir = output_dir,ref_sample_set = ref_sample_set,distance=distance,use_tree=use_tree);return(base_counts_list)})
    aggregated_clade_base_counts_list=lapply(clade_base_counts_list,function(list) {
      chrom_pos_df=list[[1]][,c(2,3)] #Get the co-ordinates of apparent het SNPs from the 1st in the list
      list_mod=lapply(list,function(df) {res<-left_join(chrom_pos_df,df[,-1],by=c("Chr","Pos"));res[is.na(res)]<-0;return(res[,3:7])}) #Now get just the base counts at these sites
      counts_df=Reduce(function(df1,df2) {return(df1+df2)},list_mod) #Aggregate these across a pure subclade
      return(cbind(chrom_pos_df,counts_df))
    })
    het_positions<-return_heterozygous_SNPs(base_counts_list=aggregated_clade_base_counts_list)
    het_positions<-het_positions[het_positions!=Pos] #Exclude the position of the actual mutation
    outcome=check_matching_phasing_non_clonal(phase_sum1,phase_sum2,het_positions = het_positions)
    return(outcome)
  }
}

get_confirmed_heterozygous_SNPs=function(phasing_info1,phasing_info2){
  comb_df=full_join(phasing_info1,phasing_info2,by="SNP_site")
  comb_df$depth=comb_df$n_alt_reads_supporting.x+comb_df$n_alt_reads_supporting.y+comb_df$n_ref_reads_supporting.x+comb_df$n_ref_reads_supporting.y
  
  #Test if the called SNPs appear real
  true_het1<-comb_df$alt_phases_with_base.x!=comb_df$ref_phases_with_base.x
  true_het2<-comb_df$alt_phases_with_base.y!=comb_df$ref_phases_with_base.y
  het_test=mapply(FUN=function(x,y) {vec=c(x,y);vec<-vec[!is.na(vec)];if(length(vec)==0) {return(T)} else if(all(vec)) {return(T)} else {return(F)}},x=true_het1,y=true_het2)
  
  #If no true het SNPs, stop function; else filter the comb_df
  if(!any(het_test)) {
    stop(return("No SNPs appear to be heterozygous"))
  } else {
    comb_df<-comb_df[het_test,]
  }
  
  #Although heterozygosity is likely after the above test, it is not confirmed. Test for this:
  het_confirmed=apply(comb_df[,c("ref_phases_with_base.x","ref_phases_with_base.y","alt_phases_with_base.x","alt_phases_with_base.y")],1,function(x) length(unique(x[!is.na(x)]))>1)
  if(any(het_confirmed)) {
    comb_df<-comb_df[het_confirmed,]
    return(comb_df$SNP_site)
  } else {
    return(NA)
  }
}


check_for_both_alleles_confirming_ref=function(phasing_info) {
  if(is.logical(phasing_info)) {
    result<-"No nearby heterozygous SNPs to confirm"
  } else if(sum(phasing_info$n_ref_reads_against)>0 & sum(phasing_info$n_ref_reads_against)>(sum(phasing_info$n_ref_reads_supporting)*0.2)){
    result<-"Both alleles confirmed with reference allele"
  } else {
    result<-"May have biased allele sequencing or LOH - suggest further confirmation"
  }
  return(result)
}

create_mut_df=function(mut,tree,matrices) {
  #Create reference set of sample sets that form clades - used in assessing PVVs and MAVs
  all_clades=unique(tree$edge[,2])
  all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))
  
  mut_df<-dplyr::bind_rows(mapply(function(samples,clade) {return(data.frame(NV=sum(matrices$NV[mut,samples]),NR=sum(matrices$NR[mut,samples]),clades=clade))},samples=all_clade_samples,clade=all_clades,SIMPLIFY = FALSE))
  mut_df$pos_test<-apply(mut_df,1,function(x){(x[2]>=12 & (x[1]/x[2])>=0.25)|(x[2]>=8 & (x[1]/x[2])>=0.3)|(x[2]>=6 & (x[1]/x[2])>=0.4)})
  mut_df$neg_test<-apply(mut_df,1,function(x){x[1]==0 & (x[2])>=10})
  return(mut_df)
}

create_MAV_df=function(mut1,mut2,tree,matrices) {
  #Create reference set of sample sets that form clades - used in assessing PVVs and MAVs
  all_clades=unique(tree$edge[,2])
  all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))
  
  MAV_df<-dplyr::bind_rows(mapply(function(samples,clade) {return(data.frame(NV1=sum(matrices$NV[mut1,samples]),NV2=sum(matrices$NV[mut2,samples]),NR=(sum(matrices$NV[mut2,samples])+sum(matrices$NR[mut1,samples])),clades=clade))},samples=all_clade_samples,clade=all_clades,SIMPLIFY = FALSE))
  MAV_df$mut1_pos_test<-apply(MAV_df,1,function(x){(x[3]>=12 & (x[1]/x[3])>=0.2)|(x[3]>=7 & (x[1]/x[3])>=0.25)|(x[3]>=5 & (x[1]/x[3])>=0.5)})
  MAV_df$mut2_pos_test<-apply(MAV_df,1,function(x){(x[3]>=12 & (x[2]/x[3])>=0.2)|(x[3]>=7 & (x[2]/x[3])>=0.25)|(x[3]>=5 & (x[2]/x[3])>=0.5)})
  MAV_df$neg_test<-apply(MAV_df,1,function(x){sum(x[1:2])==0 & (x[3])>=10})
  
  return(MAV_df)
}

#Function to determine the "lesion path" and the fixed or "pure" subclades thrown off by the lesion
get_pure_subclades=function(mut1,mut2=NULL,lesion_node,tree,matrices) {
  if(is.null(mut2)) {test_type="PVV"} else {test_type="MAV"}
  print(paste("Testing",test_type))
  
  if(test_type=="PVV"){mut_df=create_mut_df(mut=mut1,tree=tree,matrices=matrices)} else {mut_df=create_MAV_df(mut1=mut1,mut2=mut2,tree=tree,matrices=matrices)}
  
  #1. get daughter nodes of lesion node
  lesion_children=get_node_children(lesion_node,tree=tree)
  if(length(lesion_children)>2) { #if initial_lesion_node is at site of polytomy, drop the negative branches of the polytomy
    print("Removing polytomy")
    if(test_type=="PVV") {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$pos_test[mut_df$clades%in%nodes]))})
    } else {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])|any(mut_df$mut2_pos_test[mut_df$clades%in%nodes]))})
    }
    lesion_children<-lesion_children[keep_children]
  }
  
  #Test these daughter nodes to see if they are "pure positive", "pure negative" or "mixed"
  if(test_type=="PVV") {types=get_node_types(lesion_children,mut_df,tree=tree)} else {types=get_MAV_node_types(lesion_children,mut_df,tree=tree)}
  if(sum(types=="mixed")>1) {
    stop(return("More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion"))
  } else if(length(unique(types))==1){
    stop(return(ifelse(test_type=="PVV","Not PVV","Not MAV")))  
  }
  
  names(lesion_children)<-types
  pure_subclades=lesion_children[names(lesion_children)!="mixed"]
  
  while(any(names(lesion_children)=="mixed")) {
    lesion_node=lesion_children["mixed"] #Get the new lesion node for this iteration (the "mixed" descendant of the previous lesion node)
    lesion_children=get_node_children(lesion_node,tree=tree)
    if(test_type=="PVV") {types=get_node_types(lesion_children,mut_df,tree=tree)} else {types=get_MAV_node_types(lesion_children,mut_df,tree=tree)}
    if(sum(types=="mixed")>1) {stop(return("More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion"))}
    names(lesion_children)=types
    if(!any(types=="mixed")) {
      pure_subclades=c(pure_subclades,lesion_children)
      break
    } else {
      lesion_node=lesion_children["mixed"]
      pure_subclades=c(pure_subclades,lesion_children[types!="mixed"])
    }
  }
  return(pure_subclades)
}

get_mixed_subclades=function(mut1,mut2=NULL,lesion_node,tree,matrices) {
  if(is.null(mut2)) {test_type="PVV"} else {test_type="MAV"}
  print(paste("Testing",test_type))
  
  if(test_type=="PVV"){mut_df=create_mut_df(mut=mut1,tree=tree,matrices=matrices)} else {mut_df=create_MAV_df(mut1=mut1,mut2=mut2,tree=tree,matrices=matrices)}
  
  #1. get daughter nodes of lesion node
  lesion_children=get_node_children(lesion_node,tree=tree)
  if(length(lesion_children)>2) { #if initial_lesion_node is at site of polytomy, drop the negative branches of the polytomy
    print("Removing polytomy")
    if(test_type=="PVV") {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$pos_test[mut_df$clades%in%nodes]))})
    } else {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])|any(mut_df$mut2_pos_test[mut_df$clades%in%nodes]))})
    }
    lesion_children<-lesion_children[keep_children]
  }
  
  #Test these daughter nodes to see if they are "pure positive", "pure negative" or "mixed"
  if(test_type=="PVV") {types=get_node_types(lesion_children,mut_df,tree=tree)} else {types=get_MAV_node_types(lesion_children,mut_df,tree=tree)}
  names(lesion_children)<-types
  if(sum(types=="mixed")>1) {
    stop(return("More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion"))
  } else if(sum(types=="mixed")==0){
    stop(return("No mixed daughters"))
  } else {
    return(lesion_children["mixed"])
  }
}

get_file_paths_and_project=function(dataset,Sample_ID,input_data_dir="/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/input_data/") {
  if(dataset=="MSC_fetal") {
    tree_file_path=paste0(input_data_dir,"/",dataset,"/Tree_",Sample_ID,".tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/Filtered_mut_set_annotated_",Sample_ID)
    project=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex=NA
  } else if(dataset=="EM") {
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,"_standard_rho01.tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/annotated_mut_set_",Sample_ID,"_standard_rho01")
    if(grepl("PD45534|KX004",Sample_ID)) {
      project=read.csv(paste0(input_data_dir,"/EM/KX004_project_reference.csv"),header=T)
      project<-project[,c("Sample","Project")]
      colnames(project)<-c("sample","project")
    } else {
      project_ref=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
      project_ref<-project_ref[,c(1,3)]
      colnames(project_ref)<-c("sample","project")
      sample=substr(Sample_ID,1,5)
      project=as.numeric(project_ref$project[project_ref$sample==sample])
    }
    
    sex=NA
  } else if(dataset=="KY") {
    tree_file_path=paste0(input_data_dir,"/",dataset,"/",Sample_ID,"_rmix_consense_tree_no_branch_lengths_1811.tree") 
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/Filtered_muts_",Sample_ID)
    project_ref=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project=as.numeric(project_ref$project[project_ref$sample==Sample_ID])
    sex=NA
  } else if(dataset=="MSC_BMT") {
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,"_m40_postMS_reduced_a_j_pval_post_mix.tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/annotated_mut_set_",Sample_ID,"_m40_postMS_reduced_a_j_pval_post_mix")
    project=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex_vec=c(Pair11="male",Pair13="male",Pair21="male",Pair28="female",Pair31="male",Pair40="male")
    sex=sex_vec[Sample_ID]
  } else if(dataset=="PR") {
    tree_file_path=paste0(input_data_dir,"/",dataset,"/",Sample_ID,"/snp_tree_with_branch_length_polytomised.tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/",Sample_ID,"/Filtered_muts_",Sample_ID)
    project=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project<-project[,c("sample","project")]
    sex=NA
  } else if(dataset=="MF"){
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,"_noMixed.tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/annotated_mut_set_",Sample_ID,"_pval")
    project=2305
    sex=NA
  } else if(dataset=="NW"){
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,".tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/filtered_muts_",Sample_ID)
    project=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex=NA
  } else if(dataset=="SN"){
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,".tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/Filtered_muts_",Sample_ID)
    project=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project<-project[,c("sample","project")]
    colnames(project)<-c("sample","project")
    sex=NA
  } else if(dataset=="MSC_chemo"){
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,"_m40_postMS_reduced_a_j_vaf_post_mix.tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/annotated_mut_set_",Sample_ID,"_m40_postMS_reduced_a_j_vaf_post_mix")
    project=1805
    sex="female"
  } else if(dataset=="DK"){
    tree_file_path=paste0(input_data_dir,"/",dataset,"/tree_",Sample_ID,"_m40_postMS_reduced_a_j_vaf_post_mix.tree")
    filtered_muts_path=paste0(input_data_dir,"/",dataset,"/annotated_mut_set_",Sample_ID,"_m40_postMS_reduced_a_j_vaf_post_mix")
    project=read.csv(paste0(input_data_dir,"/",dataset,"/Samples_project_ref_",dataset,".csv"),header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex=NA
  }
  return(list(tree_file_path=tree_file_path,filtered_muts_path=filtered_muts_path,project=project,sex=sex))
}

#Combining the MAVs into a consensus Ref and Alt covering the same positions is complicated by the fact that they may report slightly different sequences (e.g. if one is an SNV and the other an MNV)
#This function is to define a reference set
establish_ref_and_alt=function(Ref1,Ref2,Alt1,Alt2,Pos1,Pos2) {
  if(Ref1==Ref2) {
    Ref<-Ref1 
  } else if(Pos1!=Pos2){
    Ref_vec1=as.vector(str_split(Ref1,"",simplify=T))
    names(Ref_vec1)=seq_along(Ref_vec1)+Pos1-1
    
    Ref_vec2=as.vector(str_split(Ref2,"",simplify=T))
    names(Ref_vec2)=seq_along(Ref_vec2)+Pos2-1
    
    Ref_vec<-c(Ref_vec1,Ref_vec2[!names(Ref_vec2)%in%names(Ref_vec1)])
    Ref<-paste0(Ref_vec,collapse="")
    
    #Do the same for the Alt1 vector
    Alt1_vec=as.vector(str_split(Alt1,"",simplify=T))
    names(Alt1_vec)=seq_along(Alt1_vec)+Pos1-1
    Alt1_new_vec=c(Alt1_vec,Ref_vec[!names(Ref_vec)%in%names(Ref_vec1)])
    Alt1<-paste0(Alt1_new_vec,collapse="")
    
    #Do the same for the Alt1 vector
    Alt2_vec=as.vector(str_split(Alt2,"",simplify=T))
    names(Alt2_vec)=seq_along(Alt2_vec)+Pos2-1
    Alt2_new_vec=c(Ref_vec[!names(Ref_vec)%in%names(Ref_vec2)],Alt2_vec)
    Alt2<-paste0(Alt2_new_vec,collapse="")
    
  } else if(nchar(Ref2)>nchar(Ref1)){
    Ref<-Ref2 #The longer ref is the "new ref"
    
    #Set this up as a vector named by the position of each base
    Ref_vec=as.vector(str_split(Ref,"",simplify=T))
    names(Ref_vec)=seq_along(Ref_vec)+Pos2-1
    
    #Do the same for the Alt1 vector
    Alt1_vec=as.vector(str_split(Alt1,"",simplify=T))
    names(Alt1_vec)=seq_along(Alt1_vec)+Pos1-1
    
    #Make the new Alt1 by replacing the matching positions of the Ref vec
    Alt1_new_vec<-Ref_vec
    Alt1_new_vec[names(Alt1_vec)]<-Alt1_vec
    Alt1<-paste0(Alt1_new_vec,collapse="")
    
  } else if(nchar(Ref1)>nchar(Ref2)){
    Ref<-Ref1 #The longer ref is the "new ref"
    
    #Set this up as a vector named by the position of each base
    Ref_vec=as.vector(str_split(Ref,"",simplify=T))
    names(Ref_vec)=seq_along(Ref_vec)+Pos1-1
    
    #Do the same for the Alt1 vector
    Alt2_vec=as.vector(str_split(Alt2,"",simplify=T))
    names(Alt2_vec)=seq_along(Alt2_vec)+Pos2-1
    
    #Make the new Alt1 by replacing the matching positions of the Ref vec
    Alt2_new_vec<-Ref_vec
    Alt2_new_vec[names(Alt2_vec)]<-Alt2_vec
    Alt2<-paste0(Alt2_new_vec,collapse="")
  }
  return(list(Ref=Ref,Alt1=Alt1,Alt2=Alt2))
}

#This looks through the output from the "Phase_MAVs.R" script and extracts a phasing summary
extract_MAV_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$positive_subclade_res)) {
    stop(return("No result"))
  } else {
    res<-list$positive_subclade_res
  }
  
  if(class(res)=="character") {
    stop(return(res))
  } else if(class(res)=="list"){
    res_vec=unlist(res)
  }
  
  if(length(res_vec)==1) {
    stop(return(res_vec))
  } else if(length(unique(res_vec))==1) {
    return(res_vec[1])
  } else {
    res_vec_MAV<-res_vec[grepl("pure_mut1",names(res_vec))&grepl("pure_mut2",names(res_vec))]
    if(length(res_vec_MAV)==0) {
      stop(return("Unable to confirm phasing"))
    } else if(any(res_vec_MAV=="Same phasing confirmed")) {
      return("Same phasing confirmed")
    } else {
      return("Unable to confirm phasing")
    }
  }
}

#A set of slightly fiddley functions to used in the functions to extract the phasing summary from the PVV phasing info
confirm_het_SNP=function(phasing_info_df){
  if(class(phasing_info_df)=="logical") {
    return(NULL)
  } else {
    true_het<-phasing_info_df$alt_phases_with_base!=phasing_info_df$ref_phases_with_base
    result=data.frame(SNP_site=phasing_info_df$SNP_site,res=sapply(true_het,function(res) {
      if(is.na(res)) {
        return(NA)
      } else if(res) {
        return("Heterozygous")
      } else {
        return("Not heterozygous")
      }
    }))
    return(result) 
  }
}

return_het_SNPs_from_positive_clades=function(positive_subclade_phasing_info) {
  res=lapply(positive_subclade_phasing_info,confirm_het_SNP)
  res[unlist(lapply(res,is.null))]<-NULL
  comb_res=Reduce(f=function(df1,df2) {return(full_join(df1,df2,by="SNP_site"))},res)
  if(!is.null(comb_res)&!all(is.na(comb_res[,grepl("res",colnames(comb_res))]))){
    final_res=apply(comb_res[,-1,drop=F],1,function(x) {res<-x[!is.na(x)]; if(length(unique(res))==1){return(res[1])}else{return(NA)}})
    het_SNPs=comb_res$SNP_site[final_res=="Heterozygous" & !is.na(final_res)]
    if(length(het_SNPs)>0) {
      return(het_SNPs) 
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

get_alt_base=function(SNP,positive_subclade_phasing_info) {
  alts=unlist(lapply(positive_subclade_phasing_info,function(df) {
    if(class(df)=="logical") {
      return(NA)
    } else {
      return(df$alt_phases_with_base[df$SNP_site==SNP])
    }
  }))
  alts<-alts[!is.na(alts)]
  if(length(unique(alts))>1) {
    return("Conflicting results")
  } else {
    names(alts)<-NULL
    return(alts[1])
  }
}

assess_presence_of_alt_allele=function(alt_bases,negative_subclade_phasing_info) {
  het_SNP_sites=names(alt_bases)
  res=lapply(negative_subclade_phasing_info,function(df) {
    if(class(df)=="logical") {
      return("Alt allele not confirmed")
    } else if(!any(het_SNP_sites%in%df$SNP_site)) {
      return("Alt allele not confirmed")
    } else {
      res2=sapply(het_SNP_sites,function(SNP) {
        alt_base=alt_bases[SNP]
        if(SNP%in%df$SNP_site){
          if(df$ref_phases_with_base[df$SNP_site==SNP]==alt_base) {
            return("Alt allele reads present")
          } else if(df$ref_phases_with_base[df$SNP_site==SNP]!=alt_base & df$n_ref_reads_against[df$SNP_site==SNP]>0) {
            return("Alt allele reads present")
          } else {
            return("Alt allele not confirmed")
          } 
        } else {
          return(NA)
        }
      })
      res2<-res2[!is.na(res2)]
      if(any(res2=="Alt allele reads present" )) {
        return("Alt allele reads present")
      } else {
        max_reads=max(df$n_ref_reads_supporting[df$SNP_site%in%het_SNP_sites])
        return(paste("Alt allele not confirmed with maximum of",max_reads,"reads supporting the other allele"))
      } 
    }
  })
  return(res)
}

#This function extracts the positive subclade phasing info results from the results list
extract_MAV_pos_clade_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$positive_subclade_res)) {
    stop(return("No result"))
  } else {
    res<-list$positive_subclade_res
  }
  
  if(class(res)=="character") {
    stop(return(res))
  } else if(class(res)=="list"){
    res_vec=unlist(res)
  }
  
  if(length(res_vec)==1) {
    stop(return(res_vec))
  } else if(length(unique(res_vec))==1) {
    return(res_vec[1])
  } else {
    res_vec_MAV<-res_vec[grepl("pure_mut1",names(res_vec))&grepl("pure_mut2",names(res_vec))]
    if(length(res_vec_MAV)==0) {
      stop(return("Unable to confirm phasing"))
    } else if(any(res_vec_MAV=="Same phasing confirmed")) {
      return("Same phasing confirmed")
    } else {
      return("Unable to confirm phasing")
    }
  }
}

#This function extracts the positive subclade phasing info results from the results list
extract_PVV_pos_clade_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$positive_subclade_res)) {
    stop(return("No result"))
  } else {
    res_pos<-list$positive_subclade_res
  }
  
  if(class(res_pos)=="character") {
    stop(return(res_pos))
  } else if(class(res_pos)=="list"){
    res_pos_vec=unlist(res_pos)
  }
  
  if(length(res_pos_vec)==1) {
    stop(return(res_pos_vec))
  } else if(length(unique(res_pos_vec))==1) {
    return(res_pos_vec[1])
  } else {
    if(any(res_pos_vec=="Same phasing confirmed")) {
      return("Same phasing confirmed in at least one subclade")
    } else if(any(res_pos_vec=="Non-matching phasing confirmed")) {
      return("Non-matching phasing confirmed in at least one subclade")
    } else {
      return("Unable to confirm phasing")
    }
  }
}

#This function examines the read counts of the negative subclades of the PVV to see if they include reads that match
#the phasing of the mutant allele in the positive subclades. If they do this means (1) there is no LOH, (2) both alleles have been sequenced
extract_PVV_neg_clade_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$negative_subclade_res)) {
    stop(return("No result"))
  } else {
    res_neg<-list$negative_subclade_res
  }
  
  if(class(res_neg)=="character") {
    stop(return(res_neg))
  } else if(class(res_neg)=="list"){
    res_neg_vec=unlist(res_neg)
  } else if(is.na(res_neg)) {
    stop(return(NA))
  }
  
  if(length(res_neg_vec)==1) {
    res<-res_neg_vec
  } else if(length(unique(res_neg_vec))==1) {
    res<-res_neg_vec[1]
  } else if(any(res_neg_vec=="Both alleles confirmed with reference allele")){
    res<-"Both alleles confirmed with reference allele in at least one subclade"
  } else {
    res<-res_neg_vec
  } 
  
  if(any(res=="May have biased allele sequencing or LOH - suggest further confirmation")) {
    pos_clades=which(names(list$phasing_info_by_subclade)%in%c("pure_positive","pure_mut1","pure_mut2"))
    het_SNPs=return_het_SNPs_from_positive_clades(list$phasing_info_by_subclade[pos_clades])
    #print(het_SNPs)
    if(!is.null(het_SNPs)) {
      alt_bases<-sapply(het_SNPs,function(SNP) {get_alt_base(SNP,list$phasing_info_by_subclade[pos_clades])})
      if(all(alt_bases=="Conflicting results")) {
        stop(return("Positive clades have non-matching phasing"))
      }
      alt_bases<-alt_bases[!alt_bases=="Conflicting results"]
      neg_clades=which(names(list$phasing_info_by_subclade)=="pure_negative")
      res<-unlist(assess_presence_of_alt_allele(alt_bases,list$phasing_info_by_subclade[neg_clades]))
      res<-unique(res)
      if(length(res)>1) {
        if(any(res=="Alt allele reads present")) {
          res<-"Alt allele reads present in at least one subclade"
        } else {
          res<-paste(res,collapse=",")
        }
      }
    } else {
      res<-"No nearby heterozygous SNPs to confirm"
    }
  }
  return(res) 
}

#ASCAT copy number functions
get_cn=function(cn_summary_file){
  ##cat("opening ",cn_summary_file,"\n")
  if(!file.exists(cn_summary_file)){
    warning(sprintf("%s: does not exist",cn_summary_file))
    return(NULL)
  }
  cn=read.csv(cn_summary_file,header = FALSE)
  cn$start=cn$V3
  cn$end=cn$V4
  cn$chr=cn$V2
  ###V7=total copy number
  ## V8=minor allele copy number
  cn$major=cn$V7-cn$V8
  cn$minor=cn$V8
  cn[,-grep("^V",colnames(cn))]
}
get_ASCAT_minor_allele_cn=function(Chrom,Pos,sample,project){
  if(is.data.frame(project)) {
    sample_project<-project$project[project$sample==sample]<-project$project[project$sample==sample]
  } else {
    sample_project<-project
  }
  file = paste0("/nfs/cancer_ref01/nst_links/live/", sample_project, "/", sample, "/", sample, ".ascat_ngs.summary.csv")
  cn=get_cn(file)
  if(!is.null(cn)) {
    minor_allele_cn=cn$minor[cn$chr==Chrom & cn$start<Pos & cn$end>Pos]
    return(minor_allele_cn)
  } else {
    return(NA)
  }
}

get_mean_ASCAT_minor_allele_cn=function(Chrom,Pos,samples,project) {
  cn_vec=sapply(samples, function(sample) {
    #print(sample)
    cn=get_ASCAT_minor_allele_cn(Chrom = Chrom,Pos=Pos,sample=sample,project=project)
    return(cn)
  })
  #print(cn_vec)
  return(mean(cn_vec,na.rm = T))
}

##FUNCTIONS FOR THE ANALYSIS OF LESION SEGREAGATION DATA
#estimate the parameters pf
estimate_gamma_params=function(value_vec,log_rate_range=c(-2,1),shape_range=c(1,5)) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rate_vec = 10^(seq(log_rate_range[1],log_rate_range[2],by=0.1)) # rho will be bounded within 1e-6 and 0.89
  shape_vec=seq(shape_range[1],shape_range[2],0.05)
  params_grid=expand_grid(rate_vec,shape_vec)
  ll = sapply(1:nrow(params_grid), function(i) {shape=params_grid$shape_vec[i]; rate=params_grid$rate_vec[i];sum(dgamma(x=value_vec, shape=shape,rate=rate,log = T))})
  return(params_grid[which.max(ll),])
}

#Updated version of the squash_tree function that allows you to squash from the root, as well as from the tips
squash_tree=function(tree,cut_off=50,from_root=F) {
  if(from_root){
    idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
    new_edge_lengths=nodeHeights(tree)[idxs_to_squash,2]-cut_off #work-out the edge lengths that these should be such that they finish at the cut-off
    tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
    
    tree$edge.length[nodeHeights(tree)[,2]<=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
    return(tree)
  } else {
    tree$edge.length[nodeHeights(tree)[,1]>=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
    idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
    new_edge_lengths=cut_off - nodeHeights(tree)[idxs_to_squash,1] #work-out the edge lengths that these should be such that they finish at the cut-off
    tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
    return(tree)
  }
}

#First version of the sharedness stat - as per NW. However, this will give higher values of sharedness with smaller trees
calculate_sharedness_stat=function(tree) {
  prop_samples<-sapply(tree$edge[,2],function(node) {
    prop_samples<-length(getTips(tree,node))/length(tree$tip.label)
    return(prop_samples)
  })
  mean_w<-weighted.mean(x=prop_samples,w=tree$edge.length)
  return(mean_w)
}

#Second version of the sharedness stat. Minus 1 from the numerator & denominator.
calculate_sharedness_stat_2=function(tree) {
  prop_samples<-sapply(tree$edge[,2],function(node) {
    prop_samples<-(length(getTips(tree,node))-1)/(length(tree$tip.label)-1)
    return(prop_samples)
  })
  mean_w<-weighted.mean(x=prop_samples,w=tree$edge.length)
  return(mean_w)
}

#Count the number of internal nodes above a certain height i.e. for calculating the "post-developmental nodes"
count_internal_nodes=function(tree,cut_off=50,type="above"){
  nodeheights=nodeHeights(tree)
  if(type=="above") {
    internal_nodes_above_cutoff=sum(nodeheights[,2]>cut_off & !tree$edge[,2]%in%1:length(tree$tip.label))
    return(internal_nodes_above_cutoff)
  } else if(type=="below"){
    internal_nodes_below_cutoff=sum(nodeheights[,2]<cut_off & !tree$edge[,2]%in%1:length(tree$tip.label))
    return(internal_nodes_below_cutoff)
  }
}

##
get_coords_of_mutated_bases=function(mut_ref) {
  Pos<-as.numeric(str_split(mut_ref,pattern = "-",simplify=T)[,2])
  Ref<-str_split(mut_ref,pattern = "-",simplify=T)[,3]
  covered_coords<-Pos:(Pos+nchar(Ref)-1)
  return(covered_coords)
}

nodeHeights = function (tree, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  nHeight <- function(tree) {
    tree <- reorder(tree)
    edge <- tree$edge
    el <- tree$edge.length
    res <- numeric(max(tree$edge))
    for (i in seq_len(nrow(edge))) res[edge[i, 2]] <- res[edge[i, 
                                                               1]] + el[i]
    res
  }
  nh <- nHeight(tree)
  return(matrix(nh[tree$edge], ncol = 2L) + ROOT)
}

nodeheight=function (tree, node, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (node == (Ntip(tree) + 1)) 
    h <- 0
  else {
    a <- setdiff(c(getAncestors(tree, node), node), Ntip(tree) + 
                   1)
    h <- sum(tree$edge.length[sapply(a, function(x, e) which(e == 
                                                               x), e = tree$edge[, 2])])
  }
  h + ROOT
}

getAncestors=function (tree, node, type = c("all", "parent")) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  type <- type[1]
  if (type == "all") {
    aa <- vector()
    rt <- Ntip(tree) + 1
    currnode <- node
    while (currnode != rt) {
      currnode <- getAncestors(tree, currnode, "parent")
      aa <- c(aa, currnode)
    }
    return(aa)
  }
  else if (type == "parent") {
    aa <- tree$edge[which(tree$edge[, 2] == node), 1]
    return(aa)
  }
  else stop("do not recognize type")
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

#Note - only apply the "get_edge_from_tree" option if starting from an SNV only tree. Function will assume that all the existing edge length is SNVs.
correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
      nSNV=tree$edge.length[tree$edge[,2]==node]
    } else {
      nSNV = sum(details$node == node & details$Mut_type == "SNV")
    }   		
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)  
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
  }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

get_subset_tree = function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field] == value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

get_mut_burden = function(tree) {
  mut_burden = nodeHeights(tree)[tree$edge[,2] %in% 1:length(tree$tip.label),2]
  return(mut_burden)
}

get_mut_burden_stats = function(tree) {
  mut_burden = get_mut_burden(tree)
  cat(paste("Mean mutation burden is", round(mean(mut_burden),digits = 1),"\n"))
  cat(paste("Range of mutation burden is", round(range(mut_burden)[1],digits = 1),"to",round(range(mut_burden)[2],digits = 1),"\n"))
  cat(paste("Standard deviation of mutation burden is", round(sd(mut_burden),digits = 1),"\n"))
}

#Function to calculate the absolute minimum number of clones by counting the number of times a parent node is shared, but
#a daughter node is recipient only (this would give you the number of extant transplanted clones if had full phylogeny)
get_minimum_clones=function(tree,donor_ID,recip_ID){
  shared_node_test=function(tree,node,donor_ID,recip_ID) {
    node_samples=getTips(tree,node)
    n_donor=sum(grepl(donor_ID,node_samples))
    n_recip=sum(grepl(recip_ID,node_samples))
    sharing_info=ifelse(n_donor>0&n_recip>0,"shared",ifelse(n_donor>0,"donor","recipient")) 
  }
  N=dim(tree$edge)[1]
  by_node=sapply(1:N,function(i) {
    node=tree$edge[i,2]
    sharing_info=shared_node_test(tree,node,donor_ID,recip_ID)
    if(sharing_info=="shared") {
      daughter_nodes=get_node_children(node,tree = tree)
      evidence_of_clone=0
      for(i in daughter_nodes) {
        daughter_sharing=shared_node_test(tree,node=i,donor_ID,recip_ID)
        if(daughter_sharing=="recipient") {
          evidence_of_clone=sum(1,evidence_of_clone)
        }
      }
    } else {
      evidence_of_clone=0
    }
    return(evidence_of_clone)
  })
  total_clones=sum(unlist(by_node))
  return(total_clones)
}

#Setup functions for AMOVA
# function to perform amova.
amova.fn <- function(distmat, groupnames, cell_key) {
  groupnums <- length(groupnames)
  dw <- c()
  cellnums <- c()
  for (i in 1:groupnums) {
    tgroup <- groupnames[i]
    tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type==tgroup])
    cellnums <- c(cellnums, length(tcells))
    dw <- c(dw, sum(distmat[tcells, tcells]))
  }
  tdist <- distmat[colnames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames], rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames]]
  dap <- (sum(tdist) - sum(dw))/2
  
  dfAP <- length(groupnames) - 1 
  dfWP <- sum(cellnums-1)
  N <- sum(cellnums)
  SSwp <- sum(dw/(cellnums*2))
  SSap <- sum(((dw + dap)/(2*N)) - (dw/(2*cellnums)))
  
  MSwp <- SSwp/dfWP
  MSap <- SSap/dfAP
  nc <- (N - (sum(cellnums^2)/N))/dfAP
  varwp <- MSwp
  varap <- (MSap - MSwp)/nc
  obsphi <- varap/(varwp + varap)
  return(obsphi)
}

# function to randomise sample labels and repeat
randamova.fn <- function(distmat, groupnames, cell_key) {
  
  # change added 2018.01.22: when randomizing, only include the part of the distance matrix that involves the cell types being considered.
  tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames])
  #
  
  randmat <- distmat[tcells, tcells]
  colnames(randmat) <- sample(colnames(randmat))
  rownames(randmat) <- colnames(randmat)
  randphi <- amova.fn(distmat=randmat, groupnames=groupnames, cell_key=cell_key)
  return(randphi)
}

# function tying it all in together
amovapval.fn <- function(distmat, groupnames, cell_key, iterations, plottitle) {
  # calculate observed
  obsphi <- amova.fn(distmat=distmat, groupnames=groupnames, cell_key=cell_key)
  # calculate null
  randphis <- sapply(1:iterations, function(cell) randamova.fn(distmat = distmat, groupnames=groupnames, cell_key=cell_key))
  # calculate pval
  pval <- length(which(randphis>obsphi))/length(randphis) 
  
  hist(randphis, col="grey", 100, main=plottitle, xlab="Phi statistic")
  abline(v=obsphi, col="red", lwd=2)
  legend("topright", legend=paste0("Observed\n p = ", signif(pval,digits=2)), lwd=2, col="red", bty="n")
}


require("RColorBrewer")
get_idx_for_node=function(details,node){
  which(details$node==node)
}

get_edge_info=function(tree,details,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(details,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}

add_annotation=function(tree,details=NULL,matrices=NULL,annot_function,...){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree,details,matrices,tree$edge[i,2],...))
}

add_binary_proportion=function(tree,##<< enhanced phylo returned from plot_tree
                               details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                               matrices,##<< a list of matrices parallel to details with columns named by tip labels
                               node,
                               bfield,
                               b.add.line=TRUE,
                               b.add.text=FALSE,
                               ...
){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  
  bdat=details[[bfield]][info$idx]
  if(is.null(bdat) || class(bdat)!="logical"){
    stop("Error in provided bfield (does it exist and is it boolean?)")
  }
  pass=sum(bdat,na.rm=TRUE)
  fail=sum(!bdat,na.rm=TRUE)
  tot=pass+fail
  ycross=info$yb+(fail/tot)*(info$yt-info$yb)
  ##Could add in a third category NA
  #missing=sum(is.na(bdat))
  if(b.add.line){
    arrows(y0=info$yb,y1=ycross,x0=info$x,length = 0,col="black",lend=1,...)
    arrows(y0=ycross,y1=info$yt,x0=info$x,length = 0,col="red",lend=1,...)
  }
  if(b.add.text){
    text(y=ycross,x=info$x,label=pass,pos = 4,offset = 0,...)
  }
}

add_simple_labels=function(tree,##<< enhanced phylo returned from plot_tree
                           details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                           matrices,##<< a list of matrices parallel to details with columns named by tip labels
                           node,##<< Node (see details)
                           query.field,##<< Name of column in details to query against
                           query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                           label.field,##<< Name of column in details specifying the label text.
                           cex.label=1,
                           b.add.label=TRUE,
                           b.add.marker=TRUE,
                           ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>50){
    stop("too many variants to annotate")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  pch=query.allowed.df$pch[idx.match]
  
  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      y=0.5*(info$yb+info$yt)
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,...)
    }
    if(b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

add_vaf=function(tree,##<< enhanced phylo returned from plot_tree
                 details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                 matrices,##<< a list of matrices parallel to details with columns named by tip labels
                 node,
                 samples=NULL,
                 b.plot.bars=TRUE,
                 lwd.rect=1,
                 min.depth=1,
                 vc.field,
                 vc.df,
                 filter.on=NULL,
                 ...
){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  #browser()
  ##Can do additional filter e.g. missense 
  #if(!is.null(filter.on)){
  #info$idx=info$idx,3)#info$idx.in.details[which(details$VC=="missense")]
  # }
  
  if(length(info$idx)==0){
    return(NULL)
  }
  if(b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  if(is.null(samples)){
    samples=info$samples
  }
  cat(length(info$idx),"\n")
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=matrices$mtr[info$idx,samples],
                  dep=matrices$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,details[info$idx,])
  df=df[which(df$dep>=min.depth),]
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  
  df=df[order(df$vaf),]
  yd=info$yt-info$yb
  
  
  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }
  
  if(!b.plot.bars){
    r=0.8  ##r>0.5 will cause likely overlap problems
  }else{
    r=0.4
  }
  #arrows(x0=info$x-w,x1=info$x-w+df$vaf*2*w,y0=y,lend=1,length=0,col="black")
  #arrows(x0=info$x-w+df$vaf*2*w,x1=info$x+w,y0=y,lend=1,length=0,col="grey")
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c("black","white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  if( !b.plot.bars){
    
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  min.mean.vaf=0.45
  z1=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
  z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
  max.p.value=max(z1$p.value,z2$p.value)
  txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
  if(z1$p.value>0.05 & z2$p.value<0.05) {
    border.color="green"
  } else if(max.p.value<0.05){
    if(max.p.value<0.05/dim(tree$edge)[1]){
      border.color="red"
    }else{
      border.color="blue"
    }
  }else{
    border.color="darkgrey"
  }
  
  
  rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=lwd.rect)
  if(border.color!="darkgrey"){
    text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  
  #df
  
}


plotBars=function(x,y,radius,col,prop,border="black",width=1){
  #cat(prop,"\n")
  if(width<2){
    arrows(x0 = x-radius,y0=y,x1=x-radius+2*radius*prop[1],col="darkgrey",lend=2,length=0)
    arrows(x0 = x-radius+2*radius*prop[1],y0=y,x1=x-radius+2*radius,col=rgb(0.98,0.98,0.98),lend=2,length=0)
  }else{
    rect(xleft = x-radius,xright =x-radius+2*radius*prop[1],ybottom = y-width/2,ytop=y+width/2,border = NA,col="darkgrey")
    rect(xleft =  x-radius+2*radius*prop[1],xright =x+radius,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(0.98,0.98,0.98))
  }
  1
}

plotPie=function(x,y,radius,col,prop,llwd=0.5,border="black",width=NA){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  
  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(length(prop)>1){
    for(i in 2:length(cutpoint)){
      polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              border=border,col=col[i-1],lwd=llwd)
    }
  }else{
    i=2
    polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            border=border,col=col[1],lwd=llwd)
  }
  yscale
}

plot_tree_vaf=function(tree,details,matrices,samples=NULL){
  tree=plot_tree(tree)
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = FALSE)
                     }
  )
  ##Post process res to add legend
  
}


plot_tree_labels_genes=function(tree,details,query.field="GENE",label.field="GENE",genes=c("JAK2","CBL","TET2","DNMT3A"),cex.label=1){
  qdf=data.frame(value=genes,col=rainbow(length(genes)),pch=19)
  plot_tree_labels(tree,details,
                   query.allowed.df = qdf,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}



##Gets unique colour pch combos and returns in dataframe with columns "col" and "pch"
get_color_pch_df=function(n){
  pch.list=c(18,17,16,15,0:6)
  if(n>length(pch.list)*8){
    stop("Too many colours requested")
  }
  cols=rep(RColorBrewer::brewer.pal(8,"Set1"),times=length(pch.list))
  pch=rep(pch.list,each=8)
  data.frame(col=cols,pch=pch,stringsAsFactors = FALSE)[1:n,]
  
}

get_qdf=function(values){
  if(length(values)>length(unique(values))){
    stop("get_qdf: please provide values without duplication")
  }
  cbind(data.frame(value=values,stringsAsFactors = FALSE),
        get_color_pch_df(length(values)))
}

plot_tree_labels_consequence=function(tree,details,consequences,
                                      query.allowed.df=get_qdf(consequences),
                                      query.field="VC",
                                      label.field="GENE",
                                      cex.label=1){
  ##qdf=get_qdf(consequences)
  plot_tree_labels(tree,details,
                   query.allowed.df = query.allowed.df,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}

#Useful function from online for drawing background boxes to the labels (then used in the "add_simple_labels_line" function)
boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA, 
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5, 
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){
  
  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  
  ## Is y provided:
  if (missing(y)) y <- x
  
  ## Recycle coords if necessary:    
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]           
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }       
  }
  
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)            
    }        
  } else {
    adj <- c(0.5, 0.5)
  }
  
  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }       
  } else {
    offsetVec <- c(0, 0)
  }
  
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]
  
  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth    
  graphics::rect(xleft = xMid - rectWidth/2, 
                 ybottom = yMid - rectHeight/2, 
                 xright = xMid + rectWidth/2, 
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg)
  
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font, 
                 adj = c(0.5, 0.5))    
  
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }    
}

add_simple_labels_line=function(tree,##<< enhanced phylo returned from plot_tree
                                details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                                matrices,##<< a list of matrices parallel to details with columns named by tip labels
                                node,##<< Node (see details)
                                query.field,##<< Name of column in details to query against
                                query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                                label.field,##<< Name of column in details specifying the label text.
                                cex.label=2,
                                b.add.label=TRUE,
                                b.add.marker=TRUE,
                                lty=1,
                                lwd=3,
                                ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>1){
    print("Some branches have multiple variants, plotting only the first variant of each branch")
    query.value=details[[query.field]][idx]
    vlabels=paste0(details[[label.field]][idx],collapse="\n")
  } else {
    query.value=details[[query.field]][idx]
    idx.match=match(query.value,query.allowed.df$value)
    cols=query.allowed.df$col[idx.match]
    vlabels=details[[label.field]][idx]
  }
  
  ## spread out
  N=ifelse(length(idx)>0,1,0)
  ##Vertical offset so that labels sit slightly above the markers.
  if(N>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col="red",lend=1,lwd=lwd,lty=lty,...)
    if(b.add.label){
      #text(rep(info$x,N),y=info$yb+0.5*(info$yt-info$yb),labels = vlabels,pos = 2,offset = 0.25,cex=cex.label)
      boxtext(info$x-1,info$yb+runif(1,min=0.35*(info$yt-info$yb),max=0.65*(info$yt-info$yb)),col.bg="white",border.bg="black",padding = c(0.5, 5),labels = vlabels,pos=2,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}


plot_tree_labels=function(tree,details,
                          query.field="VC",
                          type="label",
                          query.allowed.df=data.frame(value=c("nonsense","frameshift"),
                                                      col=c("red","black"),pch=c(17,18)
                          ),
                          label.field="GENE",
                          cex.label=1,
                          lty=1,
                          lwd=1){
  if(type=="label") {res=add_annotation(tree,
                                        details,list(),
                                        function(tree,details,matrices,node){
                                          add_simple_labels(tree,details,matrices,node,
                                                            query.field =query.field,
                                                            query.allowed.df = query.allowed.df,
                                                            label.field = label.field,
                                                            cex.label =cex.label)})
  with(query.allowed.df,legend("topleft",legend=value,col=col,pch=pch))
  }
  if(type=="line") {res=add_annotation(tree,
                                       details,list(),
                                       function(tree,details,matrices,node){
                                         add_simple_labels_line(tree,details,matrices,node,
                                                                query.field =query.field,
                                                                query.allowed.df = query.allowed.df,
                                                                label.field = label.field,
                                                                cex.label=cex.label,
                                                                lty=lty,
                                                                lwd=lwd)})}
}


plot_tree_vaf=function(tree,details,matrices,samples=NULL,b.plot.bars =TRUE,filter.on=NULL){
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = b.plot.bars,filter.on=filter.on)})
}

#add_vaf_bar
#add_vof_pie
#add_label

#add_var_col function which colours each mutation line according to a numeric variable (scaled between 0 and 1)
library(dichromat)
add_var_col=function(tree, ##<< enhanced phylo returned from plot_tree
                     details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                     matrices,##<< a list of matrices parallel to details with columns named by tip labels
                     node,
                     var_field,
                     pval_based=FALSE,
                     b.add.line=TRUE,
                     colours = c("black","green","red"),
                     scale_muts_to_branch=TRUE,
                     ...){
  
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  muts_on_edge=length(info$idx.in.details)
  edge_length=tree$edge.length[tree$edge[,2]==node]
  
  if(muts_on_edge > 0 & edge_length>0) {
    if(var_field == "vaf") {
      NV_vec = rowSums(matrices$mtr[info$idx.in.details,info$samples,drop =FALSE])
      NR_vec = rowSums(matrices$dep[info$idx.in.details,info$samples,drop = FALSE])
      chroms = unlist(lapply(strsplit(names(NV_vec),split = "-"),function(x) return(x[1])))
      if(pval_based) {
        bdat=sapply(1:length(info$idx.in.details), function(i) binom.test(NV_vec[i],NR_vec[i],p = ifelse(chroms[i] %in% c("X","Y"),0.95,0.5), alternative = "two.sided")$p.value)
      } else {
        bdat = NV_vec/NR_vec
        bdat[chroms %in% c("X","Y")] <- bdat[chroms %in% c("X","Y")]/2
      }
    } else {
      bdat=details[[var_field]][info$idx]
      if(is.null(bdat) || class(bdat)!="numeric"){
        stop("Error in provided bfield (does it exist and is it numeric?)")
      }
    }
    bdat = sort(bdat, decreasing = TRUE)
    if(scale_muts_to_branch) {
      mut_unit_of_edge=edge_length/muts_on_edge
    } else {
      mut_unit_of_edge=1
    }
    ##Could add in a third category NA
    #missing=sum(is.na(bdat))
    if(b.add.line){
      y0_next = info$yt
      for(i in 1:muts_on_edge) {
        arrows(y0=y0_next,y1=(y0_next - mut_unit_of_edge),x0=info$x,length = 0,col=col.scale[ceiling(100*bdat[i])],lend=1,...)
        y0_next = y0_next - mut_unit_of_edge
      }
    }
  }
}



highlight_nodes=function(tree,details,matrices,node,nodes,col="red",...) {
  info=get_edge_info(tree,details,node=node)
  if(node %in% nodes){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=col,lend=1,...)
  }
}

#The "plot_node_number" function
plot_node_number = function(tree,details,matrices,node,cex=0.4) {
  info=get_edge_info(tree,details,node)
  text(info$x,info$yb,node,cex = cex,col="black",font=2)
}

#Plot the tip point colour according to whether is donor or recipient
plot_d_or_r_tip_point = function(sample,tree,details,donor_ID,recip_ID,cols=c("dark green","red")) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  tip_col=ifelse(grepl(donor_ID,sample),cols[1],cols[2])
  points(x=info$x,y=info$yb,type="p",pch=20,bg=tip_col,col=tip_col)
}

plot_category_tip_point = function(sample_ID,tree,details=NULL,cat_df,cat_name="cat",cols=RColorBrewer::brewer.pal(8,"Set1"),col="black",...) {
  all_categories<-cat_df%>%pull(cat_name)%>%unique()
  
  if(!all(all_categories%in%names(cols))) {
    cols=cols[1:length(cat_df%>%pull(cat_name)%>%unique())]
    names(cols)<-sort(all_categories)
  }
  
  node=which(tree$tip.label==sample_ID)
  info=get_edge_info(tree,details,node)
  tip_col=cols[cat_df%>%filter(sample==sample_ID)%>%pull(cat_name)]
  points(x=info$x,y=info$yb,type="p",pch=21,bg=tip_col,col=col,...)
}

plot_postGT_tree=function(tree,details,matrices,node,highlight="post",sharing_cols=c("black","gray92"),cat_df){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  info=get_edge_info(tree,details,node=node)
  if(highlight=="pre"){
    n=sum(cat_df%>%dplyr::filter(Sample%in%info$samples)%>%pull(Time_point)==0)
  } else if(highlight=="post"){
    n=sum(cat_df%>%dplyr::filter(Sample%in%info$samples)%>%pull(Time_point)>0)
  }
  sharing_info=ifelse(n>0,"highlight","lowlight")
  names(sharing_cols)=c("highlight","lowlight")
  #if(length(tree$edge.length[tree$edge[,2]==node])>0){
  arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1)
  #}
}

plot_sharing_info=function(tree,details,matrices,node,donor_ID,recip_ID,sharing_cols=c("black","dark green","red"),...){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  info=get_edge_info(tree,details,node=node)
  n_donor=sum(grepl(donor_ID,info$samples))
  n_recip=sum(grepl(recip_ID,info$samples))
  sharing_info=ifelse(n_donor>0&n_recip>0,"shared",ifelse(n_donor>0,"donor","recipient"))
  names(sharing_cols)=c("shared","donor","recipient")
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1,lwd=ifelse(sharing_info=="shared",1.5,1),...)
  }
}

plot_sharing_multiple=function(tree,details,matrices,node,sharing_cols=c("black","dark green","red","orange","brown","green"),...){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  categories=unique(tree$tip.label)
  if(length(categories)>5) {stop("Too many tip label categories")}
  info=get_edge_info(tree,details,node=node)
  if(length(unique(info$samples))>1) {
    sharing_info <- "shared"
  } else {
    sharing_info <- unique(info$samples)
  }
  sharing_cols=sharing_cols[1:(1+length(categories))]  	
  names(sharing_cols)=c("shared",categories)
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1,...)
  }
}

plot_mut_vaf_by_branch=function(tree,
                                details,
                                matrices,
                                node,
                                mut,
                                colours=c("black","green","red"),
                                cex=0.4,
                                #show_pval=FALSE,
                                ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  #Get the allocated node number
  allocated_node=details$node[details$mut_ref==mut]
  expected_samples=get_edge_info(tree,details,allocated_node)$samples
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  variant_reads=sum(matrices$NV[mut,samples])
  total_reads=sum(matrices$NR[mut,samples])
  vaf=variant_reads/total_reads
  
  #Plot the vaf on the tree using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=col.scale[ceiling(100*vaf)],lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(variant_reads>0|node %in% which(tree$tip.label %in% expected_samples))
    text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
  
  #Print the mutation name and log10(pvalue)
  text(1,1,pos=4,mut)
  #if(show_pval&"pval"%in%colnames(details)){text(1,50,pos=4,paste0("Log10 p-value for allocated node:",round(log10(details$pval[details$mut_ref==mut]))))} 
}

plot_MAV_mut=function(tree,
                      details,
                      matrices,
                      node,
                      lesion_node=NA,
                      mut1,
                      mut2=NULL,
                      colours=c("darkgray","#B3262A","#008BCE"),
                      cex=0.4,
                      #show_pval=FALSE,
                      ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  mut_colfunc = colorRampPalette(colours[-1])
  mut_colscale = mut_colfunc(11)
  
  #Get the tips within the lesion node
  if(!is.na(lesion_node)) {
    expected_samples=getTips(tree,lesion_node)
  }
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  mut1_reads=sum(matrices$NV[mut1,samples])
  if(is.null(mut2)){mut2_reads=0}else{mut2_reads=sum(matrices$NV[mut2,samples])}
  variant_reads=mut1_reads+mut2_reads
  total_reads=variant_reads+(sum(matrices$NR[mut1,samples]) - mut1_reads)
  wt_reads=total_reads-variant_reads
  
  if((mut1_reads+mut2_reads)!=0) {
    if(is.null(mut2)) {
      base_col=colours[2]
    } else {
      base_col=mut_colscale[1+round(10*mut1_reads/(mut1_reads+mut2_reads),digits=0)] #How red or blue should the "mut" element of the colour be
    }
    final_colfunc=colorRampPalette(c("light gray",base_col))
    final_colscale=final_colfunc(101)
    branch_col=final_colscale[1+round(100*(mut1_reads+mut2_reads)/total_reads)] #Now how "concentrated" should the mut colour be
  } else {
    branch_col="light gray"
  }  
  
  #Plot the branches using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=branch_col,lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(!is.na(lesion_node)){
    if(node %in% which(tree$tip.label %in% expected_samples)) {
      text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
    }
  }
  
  #Print the mutation name
  if(node==1) { #Do for node==1 so that only prints the name once
    if(is.null(mut2)) {
      text(1,1.025*max(tree$coords$a1),pos=4,mut1)
    } else {
      text(1,1.025*max(tree$coords$a1),pos=4,paste0(mut1,"/",strsplit(x=mut2,split="-")[[1]][4]))
    }
  }
}

confirm_PVV_phylogeny=function(tree,
                               details,
                               matrices,
                               node,
                               mut,
                               PVV_mut=NULL,
                               lesion_node,
                               colours=c("black","green","red"),
                               cex=0.4,
                               #show_pval=FALSE,
                               ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  #Get the samples within the lesion node
  lesion_node_samples=getTips(tree,lesion_node)
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  variant_reads=sum(matrices$NV[mut,samples])
  total_reads=sum(matrices$NR[mut,samples])
  vaf=variant_reads/total_reads
  
  #Plot the vaf on the tree using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=col.scale[ceiling(100*vaf)],lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(variant_reads>0|node %in% which(tree$tip.label %in% lesion_node_samples))
    text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
  
  #Print the mutation name and log10(pvalue)
  if(node==1) {text(1,1,pos=4,paste("Confirmation of phylogeny for",PVV_mut,":", mut))}
  #if(show_pval&"pval"%in%colnames(details)){text(1,50,pos=4,paste0("Log10 p-value for allocated node:",round(log10(details$pval[details$mut_ref==mut]))))} 
}


###FUNCTIONS TO GO WITH THE APE "PLOT.PHYLO" FUNCTION

highlight_groups=function(tree,group1,group2,cols=c("#17698E","#17A258")) {
  any_group=c(group1,group2) #make combined list of samples assigned to at least one group
  
  #Iterate through all the nodes & work if (1) all the tips are in group1 (2) all are in group2 (3) mixture of both
  col_vec=sapply(tree$edge[,2],function(node) {
    node_daughters=getTips(tree,node)
    if(any(node_daughters%in%any_group)) {
      node_daughters<-node_daughters[node_daughters%in%any_group] #only include those tips that are assigned to a specific group
      if(all(node_daughters%in%group1)) {
        col<-cols[1]
      } else if(all(node_daughters%in%group2)) {
        col<-cols[2]
      } else {
        col<-"black" #If has tips in both groups, colour black
      }
    } else {
      col<-"lightgrey"
    }
  })
  return(col_vec)
}

highlight_samples=function(tree,samples) {
  tips=which(tree$tip.label%in%samples)
  edge_cols=sapply(tree$edge[,2],function(node) ifelse(node%in%tips,"red","black"))
  return(edge_cols)
}


##FUNCTION FOR THE SIMULATED TREES
drivers_per_sample=function(tree){
  if(is.null(tree$events)){stop(print("Need tree with events matrix recording driver information."))}
  driver_nodes=tree$events$node[tree$events$value==1 & tree$events$driverid>0]
  n_drivers<-sapply(tree$tip.label[-1],function(tip) { #Don't include the 'ancestral tip'
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    return(sum(driver_nodes%in%tip_nodes))
  })
  driver_ids<-sapply(tree$tip.label[-1],function(tip) { #Don't include the 'ancestral tip'
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    driver_ids<-tree$events%>%dplyr::filter(driverid!=0 & node %in% tip_nodes)%>%pull(driverid)
    return(paste0(driver_ids,collapse=","))
  })
  return(data.frame(Sample=tree$tip.label[-1],n_drivers=n_drivers,driver_ids=ifelse(nchar(driver_ids)==0,NA,driver_ids)))
}

##FUNCTION FOR THE DATA TREES
drivers_per_sample_data=function(tree,details){
  require(dplyr)
  if(is.null(details$is.driver)){stop(print("Need variable 'is.driver' in the details matrix"))}
  driver_nodes=details%>%dplyr::filter(is.driver==1)%>%pull(node)
  n_drivers<-sapply(tree$tip.label,function(tip) {
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    return(sum(driver_nodes%in%tip_nodes))
  })
  driver_ids<-sapply(tree$tip.label,function(tip) {
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    driver_ids<-details%>%dplyr::filter(node %in% tip_nodes & is.driver==1)%>%pull(variant_ID)
    return(paste0(driver_ids,collapse=","))
  })
  out_df<-data.frame(Sample=tree$tip.label,n_drivers=n_drivers,driver_ids=ifelse(nchar(driver_ids)==0,NA,driver_ids))
  out_df<-out_df%>%dplyr::filter(Sample!="Ancestral")
  return(out_df)
}

set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    #cat("\nsetting parent - child ",parent,child,"\n")
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    #cat("before:",length(which(tree$cedge>0)),"\n")
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    #cat("after:",length(which(tree$cedge>0)),"\n")
    #cat(parent,child,":",tree$cedge)
    tree=set_cedge(child,tree)
  }
  tree
}

##Not very efficient -- recursively calculates height as average of children's height.
get_height=function(tree,node){
  #if(is.null(tree$height)){
  #  tree$height=rep(NA,1:length(tree$edge.length))
  #}
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

##horizontal elbow
elbow=function(x0,x1,y0,y1,...){ 
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}
##vertical elbow
elbowv=function(x0,x1,y0,y1,...){ 
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)
  
}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}

plot_tree=function(tree,direction="down",cex.label=5,offset=0,plot_axis=T,title=NULL,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,vspace.reserve=0,cex.axis=1,tck=-0.02,scale=NULL){
  par(mar=c(1, 1, 1, 3) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    ymin=0-ymax*0.05-vspace.reserve*ymax
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(ymin,ymax),xlab="",ylab="")
  }else{
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,srt=270,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=c("darkgrey", "blueviolet","deeppink")[Y_loss], cex=cex.terminal.dots,pch=15)
    }
  }
  tree$direction=direction
  tree$top=top
  #scale =10
  if(is.null(scale)) {
    scales=c(0,0.5,2,5,10,100,200,500,1000,2000,5000)
    scale=scales[max(which(ymax/4>scales))]
  }
  #scale=scales[max(which(ymax/2>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  if(plot_axis){
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2, cex.axis = cex.axis, lwd = 1, lwd.ticks = 1, col = "black",tck=tck) 
  }
  if(!is.null(title)) {
    text(x = 0,y=ymax,pos = 4,labels = title)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("Mutation number",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,border = "black",lwd=0.25,col="darkred")
    
  }
  tree$ymax=ymax
  tree$vspace.reserve=vspace.reserve
  tree
}

## Add heatmap with additional information under tree
add_heatmap=function(tree,heatmap,heatvals=NULL,border="white",cex.label=2){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-0.05*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
    }
    top=bot
  }
  tree
}


plot_tree_old=function(tree,direction="down",cex.label=1,offset=0,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,plot_axis=TRUE){
  par(mar=c(1, 1, 1, 3) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
  }else{
    
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,srt=90,pos = 1,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
    }
  }
  tree$direction=direction
  tree$top=top
  scales=c(0,1,10,50,100,200,500,1000,2000,5000)
  #scale=scales[max(which(ymax/2>scales))]
  scale=scales[max(which(ymax/5>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  if(plot_axis) {
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")
    
  }
  tree
}


get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}
get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}

get_sample_mutations = function(sample,tree,details,vcf_file=FALSE) {
  terminal_node <- which(tree$tip.label == sample)
  all_nodes <- get_ancestral_nodes(terminal_node,tree$edge)
  mutations <- details[details$node %in% all_nodes,]
  if(vcf_file==TRUE) {
    vcf_file = create_vcf_files(mat=mutations)
    return(vcf_file)
  } else {
    return(mutations)
  }
}

validate_colony_muts=function(colony,tree,details,NV,NR,validatable_depth_cutoff=8,pval_cutoff=0.05) {
  tree$tip.label=gsub("_hum","",tree$tip.label)
  colony_muts = get_sample_mutations(colony,tree,details = details)$mut_ref
  n_colony = length(colony_muts)
  NR_colony=NR[colony_muts,colony]
  NR_colony[NR_colony==0] <- 1
  NV_colony=NV[colony_muts,colony]
  n_colony_validatable = sum(NR_colony >= 8)
  p_colony_validatable = n_colony_validatable/n_colony
  pvals = sapply(1:length(NV_colony), function(i) {binom.test(NV_colony[i],NR_colony[i],p=ifelse(grepl("X",colony_muts[i])|grepl("Y",colony_muts[i]),0.95,0.5),alternative="less")$p.value})
  pvals_v = pvals[NR_colony>=8]
  n_colony_validated = sum(pvals_v > 0.02)
  p_colony_validated = n_colony_validated/n_colony_validatable
  
  private_muts=details_targ_full[details$node==which(tree$tip.label==colony),"mut_ref"]
  n_private = length(private_muts)
  if(n_private>0){
    NR_private = NR[private_muts,colony]
    NR_private[NR_private==0] <-1
    NV_private = NV[private_muts,colony]
    n_private_validatable = sum(NR_private >= validatable_depth_cutoff)
    p_private_validatable = n_private_validatable/n_private
    pvals = sapply(1:length(NV_private), function(i) {binom.test(NV_private[i],NR_private[i],p=ifelse(grepl("X",private_muts[i])|grepl("Y",private_muts[i]),0.95,0.5),alternative="less")$p.value})
    pvals_v = pvals[NR_private>=validatable_depth_cutoff]
    n_private_validated = sum(pvals_v > pval_cutoff)
    p_private_validated = n_private_validated/n_private_validatable
    n_likely_subclonal=sum(pvals<pval_cutoff & NR_private>=validatable_depth_cutoff & NV_private>=3) #Good number of reads to support variant being real, but not consistent with true somatic probability distribution
  } else{
    n_private_validatable=n_private_validated=n_likely_subclonal=0
    p_private_validatable=p_private_validated=NA
  }
  colony_out=data.frame(sample=colony,
                        mean_cov=mean(NR[,colony]),
                        on_existing_phylogeny=colony%in%tree$tip.label,
                        n_colony=n_colony,
                        n_colony_validatable=n_colony_validatable,
                        p_colony_validatable=p_colony_validatable,
                        n_colony_validated=n_colony_validated,
                        p_colony_validated = p_colony_validated,
                        n_private = n_private,
                        n_private_validatable=n_private_validatable,
                        p_private_validatable=p_private_validatable,
                        n_private_validated=n_private_validated,
                        p_private_validated=p_private_validated,
                        n_likely_subclonal=n_likely_subclonal)
  return(colony_out)
}

validate_mutation = function(mutation,tree,details,NV,NR,pval_cutoff=0.05,vaf_cutoff=0.3,depth_cutoff=8,counts_only=FALSE) {
  require(ape)
  if(mutation %in% rownames(NV)){
    mut_node=details$node[details$mut_ref==mutation]
    if(mut_node>length(tree$tip.label)){
      samples=gsub(pattern="_hum",replacement="",extract.clade(tree,mut_node)$tip.label)
      if(any(samples %in% colnames(NV))) {
        NV_mut=sum(NV[mutation,samples[samples %in% colnames(NV)]])
        NR_mut=sum(NR[mutation,samples[samples %in% colnames(NV)]])
        if(counts_only) {stop(return(c(NV_mut,NR_mut)))}
      }else{
        validation_result <- "No targeted samples expected to have mutation"
        stop(return(validation_result))
      }
    }else{
      samples=gsub(pattern="_hum",replacement="",tree$tip.label[mut_node])
      if(any(samples %in% colnames(NV))){
        NV_mut=NV[mutation,samples]
        NR_mut=NR[mutation,samples]
        if(counts_only) {stop(return(c(NV_mut,NR_mut)))}
      }else{
        validation_result<-"No targeted samples expected to have mutation"
        stop(return(validation_result))
      }
    }
    if(NR_mut>=depth_cutoff) {
      pval=binom.test(NV_mut,NR_mut,p=ifelse(grepl("X",mutation)|grepl("Y",mutation),0.95,0.5),alternative="less")$p.value
      vaf=sum(NV_mut)/sum(NR_mut) #Use alternative vaf cutoff for shared mutations.  Shared mutations have higher depth, and even a small bias in capture/ sequencing rates can result in failing the binomial test.
      if(pval>pval_cutoff|(vaf>=vaf_cutoff&length(samples>1)&sum(NV_mut)>=10)){validation_result <- "PASS"} else {validation_result <- "FAIL"}
    }else{
      validation_result <- "Inadequate depth"
    }
  }else{
    validation_result<-"Not in bait set"
  }
  return(validation_result)
}

view_validation_vaf_plot=function(sample, tree, details,NV,NR,depth_cutoff=0,include_muts="all") {
  tree$tip.label=gsub("_hum","",tree$tip.label)
  if(include_muts=="all") {
    mutations=get_sample_mutations(sample,tree,details)$mut_ref
  }else if(include_muts=="private"){
    mutations=details[details$node==which(tree$tip.label==sample),"mut_ref"]
  }
  mutations<-mutations[mutations%in%rownames(NV)] #Include only mutations in the bait set
  mutations<-mutations[NR[mutations,sample]>=depth_cutoff]
  if(length(mutations)<2){stop(return("<2 mutations meeting criteria - indequate to create plot\n"))}
  n_fail=sum(details$validation_results[details$mut_ref %in% mutations] == "FAIL")
  n_pass=sum(details$validation_results[details$mut_ref %in% mutations] == "PASS")
  sample_mean_depth = mean(NR[mutations,sample])
  dens <- density((NV/NR)[mutations,sample])
  plot(dens, main = sample,xlim=c(-0.05,1.05))
  abline(v = dens$x[which.max(dens$y)])
  text(0.7,
       max(dens$y) - 0.2,
       paste("Peak VAF dens=",
             round(dens$x[which.max(dens$y)], digits = 2),
             "\nDepth cutoff for validation=",
             depth_cutoff,
             "\nMean coverage of included mutations=",
             round(sample_mean_depth,digits =2),
             "\nMutations included=",
             length(mutations),
             "\nFail variants:",
             n_fail,
             "\nPass variants:",
             n_pass), col = "red", cex = 0.7)
}

#These are functions required to calculate the true daughter nodes from a polytomous tree that is structured as a dichotomous tree
get_level_daughters = function(node,tree){
  ancestral_node_height=nodeHeights(tree)[tree$edge[,2]==node,2]
  all_daughters=get_all_node_children(node,tree)
  all_daughters_ordered=tree$edge[tree$edge[,2]%in%all_daughters,2]
  all_daughters_heights=nodeHeights(tree)[tree$edge[,2]%in%all_daughters,2]
  return(all_daughters_ordered[all_daughters_heights==ancestral_node_height])
}
get_direct_daughters=function(ancestral_node,tree){
  level_ancestors=c(ancestral_node,get_level_daughters(ancestral_node,tree_targ))
  all_daughters=tree$edge[tree$edge[,1]%in%level_ancestors,2]
  true_daughters=all_daughters[!all_daughters%in%level_ancestors]
  return(true_daughters)
}

#Print the "lineage loss stats" (how much the daughter node cell fractions add up to the parent node cell fraction)
print_lineage_loss_stats = function(ancestral_node,sample,tree,details,matrices=list()) {
  ancestral_node_frac = get_node_cell_frac(ancestral_node,sample,tree=tree,details=details,matrices=matrices)
  daughters=get_direct_daughters(ancestral_node,tree=tree_targ)
  daughter_node_fracs = rep(0,length(daughters))
  for(i in 1:length(daughters)) {
    cell_frac=get_node_cell_frac(daughters[i],sample=sample,details = details,tree=tree,matrices = matrices)
    daughter_node_fracs[i] <- cell_frac
  }
  names(daughter_node_fracs) <- daughters
  print("Parent node cell fraction:")
  print(round(ancestral_node_frac,digits=3))
  print("Individual daughter node cell fractions:")
  print(round(daughter_node_fracs,digits=3))
  total_daughter_frac=sum(daughter_node_fracs)
  print("Summed daughter node cell fraction is:")
  print(round(total_daughter_frac,digits=3))
  print("Daughters/ Parents ratio is:")
  print(round(as.numeric(total_daughter_frac/ancestral_node_frac),digits=3))
}

#Function to calculate the cell fraction in a bulk sample for the mutations leading to a specific node
get_node_cell_frac = function(node,sample,tree,details,matrices){
  node_muts=details$mut_ref[details$node==node]
  if(length(node_muts)==0){
    return(NA)
  } else {
    auto_muts=node_muts[!grepl("X",node_muts)&!grepl("Y",node_muts)]
    xy_muts=node_muts[grepl("X",node_muts)|grepl("Y",node_muts)]
    auto_cell_NV = sum(matrices$mtr[auto_muts,sample])
    auto_cell_NR = sum(matrices$dep[auto_muts,sample])
    xy_cell_NV = sum(matrices$mtr[xy_muts,sample])
    xy_cell_NR = sum(matrices$dep[xy_muts,sample])
    if(auto_cell_NR+xy_cell_NR == 0) {
      return(NA)
    } else {
      cell_frac_CI = binom.test(auto_cell_NV+xy_cell_NV,auto_cell_NR+2*xy_cell_NR)$conf.int[1:2] * 2
      overall_cell_frac=sum(c(auto_cell_NV,xy_cell_NV))/sum(c((auto_cell_NR/2),xy_cell_NR))
      attr(overall_cell_frac,"conf.int") <- cell_frac_CI
      return(overall_cell_frac)
    }
  }
}

#PLOT NODE LABELS FUNCTIONS (to add to Nick's functions)
#The "plot_node_cell_frac" function
plot_node_cell_frac = function(node,sample,tree,details,matrices,cex=0.6) {
  node_cell_frac=get_node_cell_frac(node,sample,tree,details,matrices)
  info=get_edge_info(tree,details,node)
  if(!is.na(node_cell_frac) & node_cell_frac>0.002) {
    text(info$x,info$yb,round(node_cell_frac,digits=3),cex = cex,col="black",font=2)
  }
}
#The "plot_node_number" function
plot_node_number = function(tree,details,matrices,node,cex=0.4) {
  info=get_edge_info(tree,details,node)
  text(info$x,info$yb,node,cex = cex,col="black",font=2)
}
#The add_annotation_targeted is for functions that require a SAMPLE name to list the bulk sample that is being overlayed on tree (and who's read counts are in the NV/ NR matrices)
add_annotation_targeted=function(sample,tree,details,matrices,annot_function,plot_sample_name=TRUE,...){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree$edge[i,2],sample,tree,details,matrices,...))
  if(plot_sample_name) {text(20,2,sample)}
}

#Colour the mutations based on a categorial variable - in the first case a validation result
add_categorical_col=function(tree, ##<< enhanced phylo returned from plot_tree
                             details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                             matrices,##<< a list of matrices parallel to details with columns named by tip labels
                             node,
                             var_field,
                             annot=list(),
                             b.add.line=TRUE,
                             ...){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  if(length(info$idx.in.details) > 0) {
    bdat=details[[var_field]][info$idx]
    bdat_filt=bdat[bdat %in% names(annot)]
    if(is.null(bdat)){
      stop("Error in provided bfield (does it exist and is it numeric?)")
    }
    bdat_filt = sort(bdat_filt,decreasing=TRUE)
    if(length(bdat_filt)>0){
      edge_length=length(bdat_filt)
      ##Could add in a third category NA
      #missing=sum(is.na(bdat))
      if(b.add.line){
        y0_next = info$yt
        for(i in 1:edge_length) {
          arrows(y0=y0_next,y1=(y0_next - 1),x0=info$x,length = 0,col=annot[[bdat_filt[i]]],lend=1,...)
          y0_next = y0_next - 1
        }
      }
    }
  }
}

#Plot just the tip labels (or a point) for an individual sample
plot_sample_tip_label = function(tree,details,sample) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  text(sample,x=info$x,y=info$yb,cex = 0.7)
}

plot_sample_tip_point = function(tree,details,sample) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  points(x=info$x,y=info$yb,type="p",pch=21,bg="black",col="black")
}

plot_d_or_r_tip_point = function(sample,tree,details,donor_ID,recip_ID) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  tip_col=ifelse(grepl(donor_ID,sample),"dark green","red")
  points(x=info$x,y=info$yb,type="p",pch=20,bg=tip_col,col=tip_col)
}

#Nick's "get_ancestral_nodes" function. Needed for the validation functions above
get_ancestral_nodes= function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}

#Function to extract the node counts from mutations on a given branch, for a given sample
#Divides counts into auto & XY mutations, and returns as a df
get_node_read_counts=function(node,sample,tree,details,matrices,exclude_mut_indexes=NULL) {
  info=get_edge_info(tree.multi,details,node)
  
  #Exclude mutations from a branch if they match those given
  if(!is.null(exclude_mut_indexes)) {
    info$idx.in.details<-info$idx.in.details[!info$idx.in.details%in%exclude_mut_indexes]
  }
  
  #Separate out the indexes of the autosomal and XY chromosomes
  info$idx.in.details.AUTO=info$idx.in.details[!grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  info$idx.in.details.XY=info$idx.in.details[grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  
  #Store the read counts of each in the df
  NV_auto=sum(matrices$NV[info$idx.in.details.AUTO,sample])
  NR_auto=sum(matrices$NR[info$idx.in.details.AUTO,sample])
  NV_xy=sum(matrices$NV[info$idx.in.details.XY,sample])
  NR_xy=sum(matrices$NR[info$idx.in.details.XY,sample])
  
  #Return as a vector
  #return(c(NV_auto,NR_auto,NV_xy,NR_xy))
  
  #Or as a df
  df=data.frame(node=node,NV_auto=NV_auto,NR_auto=NR_auto,NV_xy=NV_xy,NR_xy=NR_xy)
  return(df)
}

#Function to bootstrap the counts (of format above), to return cell fraction bootstrapped estimates
bootstrap_counts = function(node_counts,boot_straps=1000) {
  NR_auto=node_counts$NR_auto;NV_auto=node_counts$NV_auto
  NR_xy=node_counts$NR_xy;NV_xy=node_counts$NV_xy
  if((NR_auto+NR_xy)==0) {
    cell_fracs<-rep(0,boot_straps)
  } else if(NR_xy==0){
    cell_fracs=(rbinom(boot_straps,NR_auto,prob = NV_auto/NR_auto))/(NR_auto/2)
  } else if(NR_auto==0) {
    cell_fracs=(rbinom(boot_straps,NR_xy,prob=NV_xy/NR_xy))/NR_xy
  } else {
    cell_fracs=(rbinom(boot_straps,NR_auto,prob = NV_auto/NR_auto) +
                  rbinom(boot_straps,NR_xy,prob=NV_xy/NR_xy))/sum(c((NR_auto/2),NR_xy))
  }
}

#Function to check the distribution for a given sample
check_branch_distribution=function(sample,node,tree,details,matrices,return_counts=FALSE) {
  #Get node info
  info=get_edge_info(tree,details,node)
  
  #Separate out the indexes of the autosomal and XY chromosomes
  info$idx.in.details.AUTO=info$idx.in.details[!grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  info$idx.in.details.XY=info$idx.in.details[grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  
  if(length(info$idx.in.details.AUTO)==0) {
    x_auto<-NULL
    n_auto<-NULL
  } else {
    x_auto=matrices$NV[info$idx.in.details.AUTO,sample]
    n_auto=matrices$NR[info$idx.in.details.AUTO,sample]
  }
  if(length(info$idx.in.details.XY)==0) {
    x_xy<-NULL
    n_xy<-NULL
  } else {
    x_xy=matrices$NV[info$idx.in.details.XY,sample]
    n_xy=matrices$NR[info$idx.in.details.XY,sample]
  }
  x=c(x_auto,x_xy)
  n=c(n_auto,2*n_xy)
  muts=c(details$mut_ref[info$idx.in.details.AUTO],details$mut_ref[info$idx.in.details.XY])
  idx=c(info$idx.in.details.AUTO,info$idx.in.details.XY)
  df=data.frame(mut_ref=muts,idx.in.details=idx,NV=x,NR=n)
  n[n==0]<-1
  if(length(x)>1) {
    pval=prop.test(x,n)$p.value
  } else{
    pval<-NA
  }
  names(pval)=sample
  if(return_counts){
    return(df)
  } else {
    return(pval)
  }
}

#Function to look through the muts in a branch and define counts that come from the earliest (or not earliest) binomial distribution
find_early_muts_from_branch=function(sample,node,tree,details,matrices,return_late_muts=FALSE,cols=c("red","orange","blue","purple")) {
  print(node)
  df_node=check_branch_distribution(sample,node,tree,details,matrices,return_counts = TRUE)
  
  #Perform binomial mixture model to work out which mutations are the early ones on each branch
  hist(df_node$NV/df_node$NR,col='gray',freq=F,xlab="VAF",main=paste("Node VAFs for node",node))
  lines(density(df_node$NV/df_node$NR),lwd=2,lty='dashed')
  
  if(nrow(df_node)==2) {
    max_p = which.max(df_node$NV/df_node$NR)
    early_idx=c(1:length(df_node$NV)%in%max_p)
  } else {
    res = binom_mix(df_node$NV,df_node$NR,nrange=2:min(3,length(df_node$NV)))
    for (i in 1:res$n){
      meancov = round(mean(df_node$NR))
      lines(x=(0:meancov)/meancov,
            y=meancov*res$prop[i]*dbinom(0:meancov,meancov,prob=res$p[i]),
            type="l",col=cols[i],lwd=2)
    }
    max_p=which.max(res$p)
    early_idx=res$Which_cluster==max_p
  }
  
  if(return_late_muts) {
    max_p_idx=df_node$idx.in.details[!early_idx]
  } else {
    max_p_idx=df_node$idx.in.details[early_idx]
  }
  return(max_p_idx)
}

node_lineage_loss = function(node,sample,tree,details,matrices,boot_straps,CI=0.95,display_vafs=FALSE,return_ancestral_cell_frac=FALSE) { #tree must be a multifurcating tree
  sample=gsub("_comb","",sample)
  ancestral_node_counts=get_node_read_counts(node,sample,tree,details,matrices)
  ancestral_cell_fracs = bootstrap_counts(ancestral_node_counts,boot_straps = boot_straps)
  if(display_vafs) {hist(ancestral_cell_fracs,main="Ancestral node cell fraction distribution");abline(v=median(ancestral_cell_fracs))}
  daughters=tree$edge[tree$edge[,1]==node,2]   #Find the daughter nodes - this is why must be multifurcating tree
  #Get counts from each individual node by apply lapply; then Reduce to a df
  daughter_node_counts=Reduce(rbind,lapply(daughters,get_node_read_counts,sample=sample,tree=tree,details=details,matrices=matrices))
  #Now bootstrap across the rows of the df to get df of bootstrapped cell fractions
  daughter_cell_fracs=Reduce(rbind,lapply(1:nrow(daughter_node_counts), function(i) {bootstrap_counts(daughter_node_counts[i,])}))
  
  if(display_vafs) {apply(daughter_cell_fracs,1,hist,main="Daughter cell fraction distributions")}
  
  #default value to output
  lineages_captured=colSums(daughter_cell_fracs)/ancestral_cell_fracs
  lineages_captured[lineages_captured>1]<-1
  if(display_vafs) {hist(lineages_captured,main = "Proportion of daughter lineages captured in phylogeny")}
  
  #Alternative stat for "cell fraction lost at node" - i.e. absolute values, not relative
  # cell_frac_lost=ancestral_cell_fracs-colSums(daughter_cell_fracs)
  # cell_frac_lost[cell_frac_lost<0]<-0
  # hist(cell_frac_lost)
  
  if(return_ancestral_cell_frac) {
    ancestral_cell_fracs[ancestral_cell_fracs>1] <- 1
    df=data.frame(median=median(ancestral_cell_fracs,na.rm=TRUE),lower_CI=quantile(ancestral_cell_fracs,(1-CI)/2,na.rm=TRUE),upper_CI=quantile(ancestral_cell_fracs,1-(1-CI)/2,na.rm=TRUE))
  } else {
    df=data.frame(median=median(lineages_captured,na.rm=TRUE),lower_CI=quantile(lineages_captured,(1-CI)/2,na.rm=TRUE),upper_CI=quantile(lineages_captured,1-(1-CI)/2,na.rm=TRUE))
  }
  return(df)
}

plotDonut=function(x,y,median=NA,radius,col,prop,llwd=0.5,border=NA,width=NA,plotPie=FALSE){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  
  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(!any(is.na(cutpoint))&!any(is.nan(cutpoint))) {
    if(length(prop)>1){
      for(i in 2:length(cutpoint)){
        if(cutpoint[i-1]==cutpoint[i]) {next}
        polygon(x+c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
                y+yscale*c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
                border=border,col=col[i-1],lwd=llwd)
      }
    }else{
      i=2
      polygon(x+c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
              y+yscale*c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
              border=border,col=col[1],lwd=llwd)
    }
    yscale
    #Cut out the "donut" from the centre
    if(!plotPie) {
      polygon(x+c((radius/3)*cos(seq(from=0,to=2*pi,by=(2*pi/N)))),
              y+yscale*c((radius/3)*sin(seq(0,2*pi,by=(2*pi/N)))),
              border=border,col="white",lwd=llwd) 
    }
    #Now to draw the line for the median
    if(median>0 & !is.na(median)) {
      lines(x+c((0.66*radius)*sin(seq(from=0,to=2*median*pi,by=(2*pi/N)))),
            y+yscale*c((0.66*radius)*cos(seq(0,2*median*pi,by=(2*pi/N)))),
            col="black",lwd=3,lend=1)
    }
  }
}

###THE MAIN PLOTTING FUNCTION FOR THE TREES##
#Not entirely self-sufficient.  Needs various appropriately labelled objects in the environment

generate_targ_seq_plots=function(samples,
                                 tree,
                                 details_targ,
                                 matrices,
                                 post_prob_type=c("raw","clean"),
                                 info_type=c("post.prob","cell_frac","log_cell_frac"),
                                 prob_threshold_to_include=0.5, #Probability threshold from the post.prob matrix for plotting
                                 plot_cell_frac=TRUE,
                                 plot_donut=TRUE,
                                 donut_info="cell_frac", #other option is "lineages_lost"
                                 CI=0.8,  #Confidence intervals on the pie chart, default = 80% CI
                                 radius=3.5,  #Radius of the pie charts on the plot
                                 scale_muts_to_branch=FALSE) {
  require(plotrix)
  if(post_prob_type=="raw") {
    post.prob.mat <-post.prob[details_targ$mut_ref,]
  } else if(post_prob_type=="clean"){
    post.prob.mat <-clean.post.prob[details_targ$mut_ref,]
  }
  if(info_type=="post.prob") {
    post.prob.mat[post.prob.mat<0.05] <- 0
    details_targ_full=cbind(details_targ,post.prob.mat)
  } else {
    cell_frac=calculate_cell_frac(matrices$NV,matrices$NR)
    cell_frac_present<-cbind(cell_frac,cell_frac[,gsub("_comb","",colnames(post.prob.mat)[grep("_comb",colnames(post.prob.mat))])]) #To double up the "_comb" results, so that matched the post.prob.mat
    cell_frac_present=cell_frac_present[details_targ$mut_ref,]
    colnames(cell_frac_present)<- colnames(post.prob)
    cell_frac_present[post.prob.mat<prob_threshold_to_include]<-0
    if(info_type=="cell_frac") {
      details_targ_full=cbind(details_targ,cell_frac_present)
    } else if(info_type=="log_cell_frac") {
      log_cell_frac_present=cell_frac_present
      log_cell_frac_present[cell_frac_present != 0] <- log(log_cell_frac_present[cell_frac_present != 0])#change all the non 0 vafs to the log of their vaf
      log_cell_frac_present_scaled=log_cell_frac_present
      scale_range=c(0.01,1)
      log_cell_frac_present_scaled[cell_frac_present!=0] = plotrix::rescale(log_cell_frac_present[cell_frac_present!=0],newrange = scale_range) #scale these figures between 0 and 1
      details_targ_full=cbind(details_targ,log_cell_frac_present_scaled)
    }
  }
  lims=par("usr")
  sapply(samples, function(sample) {
    tree=plot_tree(tree, cex.label = 0,lwd=0.5,plot_axis = TRUE,default_edge_color="lightgrey")
    lims=par("usr")
    if(grepl("PD",sample)) {
      text(0,lims[4]-1,paste0(lcm_smry$Tissue[lcm_smry$Sample_ID==sample]," (",sample,"): Mean depth is ",round(mean(NR[,gsub("_comb","",sample)]),digits = 2)),cex=1,pos=4)
    } else {
      text(0,lims[4]-0.5,paste0(sample,": Mean depth is ",round(mean(NR[,gsub("_comb","",sample)]),digits = 2)),cex=1,pos=4)
    }
    if(grepl("PD",sample)) {
      text(0,lims[4]-3,paste0("Library concentration was ",round(lcm_smry$Conc[lcm_smry$Sample_ID==sample],digits=0)),pos=4)
    }
    
    
    add_annotation(tree=tree,
                   details=details_targ_full,
                   matrices,
                   annot_function=function(tree,details,matrices,node) {
                     add_var_col(tree,
                                 details,
                                 matrices,
                                 node,
                                 var_field = sample,
                                 pval_based=FALSE,
                                 lwd = 5,
                                 colours=colour.scale,
                                 scale_muts_to_branch=scale_muts_to_branch)
                   }
    )
    if(plot_cell_frac) {
      add_annotation_targeted(sample,
                              tree=tree,
                              details=details_targ_full,
                              matrices=matrices,
                              annot_function=function(node,sample,tree,details,matrices,cex=0.6) {
                                node_cell_frac=get_node_cell_frac(node,gsub("_comb","",sample),tree,details,matrices)
                                info=get_edge_info(tree,details,node)
                                if(!is.na(node_cell_frac) & any(post.prob.mat[info$idx.in.details,sample]>prob_threshold_to_include)) {
                                  text(info$x,info$yb,round(node_cell_frac,digits=3),cex = cex,col="black",font=2)
                                }
                              })
    }
    if(plot_donut) {
      #Detect which nodes to plot donuts for - do it for branches with a mean clean.post.prob of >0.5
      nodes_to_check=unique(tree$edge[,2])[!unique(tree$edge[,2])%in%1:length(tree$tip.label)]
      nodes_to_include=sapply(nodes_to_check,function(node) {
        if(mean(post.prob.mat[details_targ$node==node,sample])>0.5){return(node)}else{return(NA)}
      })
      nodes_to_include<-nodes_to_include[!is.na(nodes_to_include)]
      
      if(donut_info=="cell_frac") {
        #Iterate through these nodes and plot the donuts
        for(node in nodes_to_include) {
          print(node)
          data=node_lineage_loss(node=node,sample=sample,tree = tree,details = details_targ,matrices=matrices,boot_straps = 10000,CI=CI,return_ancestral_cell_frac = TRUE)
          #Make the pie chart for plotting on the node
          df2<-data%>%dplyr::select(-median)%>%gather(key="category",value="count")
          df2$ymax = df2$count
          df2$ymin = c(0, head(df2$ymax, n=-1))
          df2=rbind(df2,data.frame(category="lineages_lost",count=(1-df2$ymax[df2$category=="upper_CI"]),ymax=1,ymin=(df2$ymax[df2$category=="upper_CI"])))
          df2$prop=df2$ymax-df2$ymin
          
          #Get the node co-ordinates
          info=get_edge_info(node=node,tree=tree,details=details_targ)
          #Plot the pie chart
          #plotDonut(info$x,mean(c(info$yb,info$yt)),median=data$median,radius=radius,col=c( "#8D8DCB" ,"#C6C6E5", "#FFFFFF"),prop=df2$prop,border="black")
          
          median_only_prop=c(data$median,1-data$median)
          
          plotDonut(info$x,mean(c(info$yb,info$yt)),radius=radius,col=c( "#08306B" ,"#FFFFFF"),prop=median_only_prop,border="black",plotPie = TRUE)
        }
      } else if(donut_info=="lineages_lost") {
        
        #Iterate through these nodes and plot the donuts
        for(node in nodes_to_include) {
          print(node)
          data=node_lineage_loss(node=node,sample=sample,tree = tree,details = details_targ,matrices=matrices,boot_straps = 10000,CI=CI,return_ancestral_cell_frac = FALSE)
          #Make the pie chart for plotting on the node
          df2<-data%>%select(-median)%>%gather(key="category",value="count")
          df2$ymax = df2$count
          df2$ymin = c(0, head(df2$ymax, n=-1))
          df2=rbind(df2,data.frame(category="lineages_lost",count=(1-df2$ymax[df2$category=="upper_CI"]),ymax=1,ymin=(df2$ymax[df2$category=="upper_CI"])))
          df2$prop=df2$ymax-df2$ymin
          
          #Get the node co-ordinates
          info=get_edge_info(node=node,tree=tree,details=details_targ)
          #Plot the pie chart
          plotDonut(info$x,info$yb,median=data$median,radius=radius,col=c( "#8D8DCB" ,"#C6C6E5", "#FFFFFF"),prop=df2$prop,border="black",plotPie = TRUE)
        }
      }
    }
  }
  )
}

#The "squash tree" function to cut the tree at any given node height.  Tree structure is maintained, but edge lengths are shortened.
squash_tree=function(tree,cut_off=50) {
  tree$edge.length[nodeHeights(tree)[,1]>=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
  idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
  new_edge_lengths=cut_off - nodeHeights(tree)[idxs_to_squash,1] #work-out the edge lengths that these should be such that they finish at the cut-off
  tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
  return(tree)
}

#Function to calculate vaf, but avoiding dividing by zero by skipping sites with depth of 0.
calculate_vaf=function(NV,NR){
  NR[NR==0]<-1
  vaf=NV/NR
  return(vaf)
}

#Function for getting trees from different stages of development - set the edges to the appropriate length and
#this function prunes the branches that should now be removed
prune_tree_of_zero_tips = function(tree) {
  current_tree=NULL
  for(i in 1:10) {
    print(i)
    if(is.null(current_tree)){current_tree<-tree} else {current_tree<-tree_pruned}
    current_tree$tip.label<-1:length(current_tree$tip.label)
    tips_to_remove=NULL
    for(j in 1:length(current_tree$tip.label)) {
      private_branch_length=current_tree$edge.length[current_tree$edge[,2]==j]
      if(private_branch_length>0) {
        next
      } else {
        ancestor=current_tree$edge[current_tree$edge[,2]==j,1]
        if(ancestor==current_tree$edge[1,1]) {
          ancestral_branch_length<-1
        } else {
          ancestral_branch_length=current_tree$edge.length[current_tree$edge[,2]==ancestor]
        }
        all_daughters=current_tree$edge[current_tree$edge[,1]==ancestor,2]
        daughter_branch_lengths<-current_tree$edge.length[current_tree$edge[,2]%in%all_daughters]
        if(ancestral_branch_length==0|all(daughter_branch_lengths==0)) {
          tips_to_remove=c(tips_to_remove,j)
        }
      }
    }
    if(is.null(tips_to_remove)) {stop(return(current_tree))}
    tree_pruned=drop.tip(current_tree,trim.internal=FALSE,current_tree$tip.label[tips_to_remove])
  }
}

#The "clean up" function to make it fit better with phylogeny
clean_up_post=function(post.prob,details,tree) { #post.prob is the named vector of posterior probs for a given sample (i.e. vector)
  
  #Define a "node map" to allow tracking to true parent nodes in polytomous tree that is structured as dichotomous
  node_map=list()
  for(node in unique(details$node)) {node_map[[node]] = get_direct_daughters(node,tree_targ)}
  
  #Set up the output prob vector to be the same as the input, but can be edited
  clean.post=post.prob
  
  for(mut in details$mut_ref) {
    #Get the node for the current mutation
    allocated_node=details$node[details$mut_ref==mut]
    info=get_edge_info(tree,details,allocated_node)
    #The the daughter node numbers
    daughters=get_direct_daughters(allocated_node,tree=tree)
    #Get the ancestral node number
    ancestor=which(unlist(lapply(node_map,function(x) allocated_node%in%x)))
    
    #Other mutations on the current branch
    other_muts_on_same_branch=details$mut_ref[details$node==allocated_node & details$mut_ref!=mut]
    mean_current_branch=mean(post.prob[other_muts_on_same_branch])
    mean_ancestral_branch_prob=mean(post.prob[details$mut_ref[details$node==ancestor]])
    daughter_probs=post.prob[details$mut_ref[details$node%in%daughters]]
    
    #Remove artefacts, i.e. mutations that don't make sense from rest of phylogeny
    if((length(ancestor)==0|mean_ancestral_branch_prob<0.1) & all(daughter_probs<0.1)) {clean.post[mut] <- 0}
    
    #Boost post.prob of mutations that have confident downstream probabilities
    if((length(ancestor)==0|mean_ancestral_branch_prob>0.5) & (mean_current_branch>0.5|length(other_muts_on_same_branch)==0) & any(daughter_probs >0.8)) {clean.post[mut] <-(post.prob[mut]+1)/2}
  }
  return(clean.post)
}


##Tim's Binomial mixture model functions

## Expectation step
estep = function(x,size,p.vector,prop.vector,ncomp){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

## Maximisation step
mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)  
}

## EM algorithm
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC)
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}

#This is the function to use for males (as is the case for both foetuses)
calculate_cell_frac=function(NV,NR) {
  #Remove 0's from the depth, to avoid dividing by 0
  NR[NR==0]<-1
  #Get vectors to select out the autosomal and XY chromosomal mutations
  XY_muts=grepl("X",rownames(NR))|grepl("Y",rownames(NR))
  #For autosomal mutations, cell frac = NV/ (NR/2)
  cell_frac=NV/(NR/2)
  #For XY mutations, cell frac = NV/NR -> replace these accordingly
  cell_frac[XY_muts,]<-(NV[XY_muts,])/(NR[XY_muts,])
  #Cell frac cannot be greater than 1, therefore if comes out as > 1 (which can happen when dividing the NR by 2), coerce to 1
  cell_frac[cell_frac>1]<-1
  return(cell_frac)
}

#Adjust the add_heatmap function to allow large heatmap under the tree
add_mut_heatmap=function(tree,heatmap,heatvals=NULL,border="white",heatmap_bar_height=0.05,cex.label=2,label.cols=NA){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-heatmap_bar_height*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      if(!is.na(label.cols[1])) {
        text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label,col=label.cols[i])
      } else {
        text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
      }
      
    }
    top=bot
  }
  tree
}


plot_96profile=function(df,colnames=c("sampleID","chr","pos","ref","alt"),genomeFile,savefile=NULL) {
  require(Rsamtools)
  require(GenomicRanges)
  require(IRanges)
  mutations = data.frame(sampleID=df[[colnames[1]]],
                         chr=df[[colnames[2]]],
                         pos=df[[colnames[3]]],
                         ref=df[[colnames[4]]],
                         trinuc_ref= "-",
                         mut=df[[colnames[5]]])
  mutations_GRange<-GRanges(mutations$chr, IRanges::IRanges(mutations$pos-1, mutations$pos+1))
  mutations$trinuc_ref = as.vector(Rsamtools::scanFa(file=genomeFile,mutations_GRange ))
  
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  if(is.null(savefile)) {
    colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
    y = freqs_full; maxy = max(y)
    h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations")
    for (j in 1:length(sub_vec)) {
      xpos = h[c((j-1)*16+1,j*16)]
      rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
      text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
    }
  } else {
    pdf(paste0(savefile,".pdf"),width = 15,height=6)
    colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
    y = freqs_full; maxy = max(y)
    h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations")
    for (j in 1:length(sub_vec)) {
      xpos = h[c((j-1)*16+1,j*16)]
      rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
      text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
    }
    dev.off()
  }
}

find_latest_acquisition_node=function(tree,pos_samples){
  #Get list of ancestral nodes for all samples
  ancestral_nodes_list=lapply(pos_samples,function(Sample) {
    get_ancestral_nodes(node = which(tree$tip.label==Sample),edge=tree$edge,exclude_root = F)
  })
  #Find nodes that are ancestral to all the samples
  common_nodes=Reduce(intersect,ancestral_nodes_list)
  #Which of these is the most recent (i.e. has the maximum node height)
  nodeheights<-sapply(common_nodes,function(node) nodeheight(tree = tree,node = node))
  MRCA_node<-common_nodes[which.max(nodeheights)]
  return(MRCA_node)
}


compare_nodes=function(x, y) 
{
  tree1 <- deparse(substitute(x))
  tree2 <- deparse(substitute(y))
  n1 <- Ntip(x)
  n2 <- Ntip(y)
  
  key1 <- makeNodeLabel(x, "md5sum")$node.label
  key2 <- makeNodeLabel(y, "md5sum")$node.label
  mk12 <- match(key1, key2)
  mk21 <- match(key2, key1)
  if (any(tmp <- is.na(mk12))) {
    nk <- sum(tmp)
  }
  if (any(tmp <- is.na(mk21))) {
    nk <- sum(tmp)
  }
  nodes1 <- which(!is.na(mk12))
  nodes2 <- mk12[!is.na(mk12)]
  
  NODES <- data.frame(nodes1 + n1,nodes2 + n2)
  names(NODES) <- c(tree1, tree2)
  return(NODES)
}

##Function to give the tree comparison stats using Quartet
get_tree_comp_stats=function(tree1,tree2) {
  if(!is.null(tree1)) {
    #Harmonise the trees by only keeping samples in both trees
    tree1<-drop.tip(tree1,tree1$tip.label[!tree1$tip.label%in%tree2$tip.label])
    tree2<-drop.tip(tree2,tree2$tip.label[!tree2$tip.label%in%tree1$tip.label])
    
    #Calculate and print the similarity metrics for the iqtree tree
    splits_comp<-Quartet::SplitStatus(tree1,tree2)
    RF_diff<-round(Quartet::SimilarityMetrics(splits_comp)[,"SymmetricDifference"],digits=5)
    
    #Quartet is unable to run the quartet status for very large trees (includes KX004)
    if(length(tree1$tip.label)>477) {
      Quartet_diff<-NA
    } else {
      quartets_comp<-Quartet::QuartetStatus(tree1,tree2)
      Quartet_diff<-round(Quartet::SimilarityMetrics(quartets_comp)[,"SymmetricDifference"],digits = 5)
    }
  } else {
    RF_diff<-NA;Quartet_diff<-NA
  }
  return(data.frame(RF_diff=RF_diff,Quartet_diff=Quartet_diff))
}

add_mut_heatmap_PVV= function(tree,heatmap,heatvals=NULL,border="white",heatmap_bar_height=0.05,cex.label=2,label.cols=NA){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-heatmap_bar_height*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      if(!is.na(label.cols[1])) {
        text(labels[i],x=-0.02,y=0.5*(top+bot),pos = 2,cex = cex.label,col=label.cols[i])
      } else {
        text(labels[i],x=-0.02,y=0.5*(top+bot),pos = 2,cex = cex.label)
      }
      
    }
    top=bot
  }
  tree
}

#Define function required later in script for comparing phylogenies
comparePhylo_and_plot=function(tree1,tree2,names,lwd=1){
  plot_comp_tree=function(tree,comp,title,col="red",lwd=lwd,tree_pos){
    tree_name=deparse(substitute(tree))
    shared_clades=comp[,tree_pos]
    edge_width=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),lwd,2*lwd))
    edge_col=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),"black",col))
    plot(tree,show.tip.label=F,direction="downwards",edge.color=edge_col,edge.width=edge_width,main=title)
  }
  comp<-compare_nodes(tree1,tree2)
  #par(mfrow=c(1,2))
  plot_comp_tree(tree1,comp=comp,title=names[1],tree_pos = 1,lwd = lwd)
  plot_comp_tree(tree2,comp=comp,title=names[2],tree_pos = 2,lwd = lwd)
}

get_RF_dist=function(tree1,tree2){
  comp<-comparePhylo(tree1,tree2)
  RF_dist=1-(length(comp$NODES[,1])/length(unique(tree$edge[,1])))
  return(RF_dist)
}

create_vaf_mat=function(NV,NR) {NR[NR==0]<-1;return(NV/NR)}

#Function to import the different trees, either on the farm, or via the mounted directories
import_tree_file=function(SampleID,data_set,type) {
  require(ape)
  if(!type%in%c("mpboot","iqtree","scite")) {stop(print("Type must be either mpboot, iqtree, or scite"))}
  root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
  
  if(type=="mpboot") {
    files=list.files(paste0(root_dir,"/lesion_segregation/input_data/",data_set,"/"),pattern = ".tree",full.names = T)
    relevant_tree_path=grep(paste0("_",SampleID),x=files,value = T)
  } else {
    if(data_set=="MSC_BMT") {
      HSCT_dir=paste0(root_dir,"/Zur_HSCT/final_versions")
      relevant_pattern=ifelse(type=="iqtree",".fa.allocated_muts.tree","scite_tree_allocated_muts_")
      relevant_tree_files=list.files(path=HSCT_dir,pattern=relevant_pattern,full.names = T)
      relevant_tree_path=grep(SampleID,relevant_tree_files,value=T)
      
    } else if(data_set=="MF") {
      relevant_tree_path=ifelse(type=="iqtree",paste0(root_dir,"/Marga_benchmarking/",SampleID,"/iqtree_tree_allocated_muts_",SampleID,".tree"),paste0(root_dir,"/Marga_benchmarking/",SampleID,"/scite_tree_allocated_muts_",SampleID,".tree"))
      tree=ape::read.tree(relevant_tree_path)
    } else if(data_set=="EM") {
      shortID=stringr::str_split(SampleID,pattern = "_",simplify = T)[,1]
      relevant_tree_path=ifelse(type=="iqtree",paste0(root_dir,"/Emily_benchmarking/",shortID,"/vaf_filtered/iqtree_tree_allocated_muts_",SampleID,".tree"),paste0(root_dir,"/Emily_benchmarking/",shortID,"/vaf_filtered/scite_tree_allocated_muts_",SampleID,".tree"))
    }
  }
  
  #Return the relevant tree if it exists
  if(file.exists(relevant_tree_path)) {
    print(paste("The",type,"tree for",SampleID,"is being read in."))
    tree=read.tree(relevant_tree_path)
    return(tree)
  } else {
    print(paste("The",type,"tree for",SampleID,"is not found at the expected location."))
    return(NULL)
  }
}

import_muts_file=function(SampleID,data_set) {
  root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/Mounts/lustre2","/lustre/scratch126/casm/team154pc/ms56")
  
  data_set_muts_files=list.files(path=paste0(root_dir,"/lesion_segregation/input_data/",data_set),pattern="annotated",full.names = T)
  filepath=grep(paste0("_",SampleID),data_set_muts_files,value=T)
  load(filepath)
  return(filtered_muts)
}

get_direct_daughters=function(ancestral_node,tree){
  level_ancestors=c(ancestral_node,get_level_daughters(ancestral_node,tree))
  all_daughters=tree$edge[tree$edge[,1]%in%level_ancestors,2]
  true_daughters=all_daughters[!all_daughters%in%level_ancestors]
  return(true_daughters)
}

#Function to create the heatmap used below the trees to show the VAFs
create_mutation_heatmap=function(mut_refs1,mut_refs2=NULL,matrices,tree,cols=c("#FF7F00","#377EB8")) {
  vaf_mat=create_vaf_mat(matrices$NV,matrices$NR)
  combined_muts=c(mut_refs1,mut_refs2)
  hm=matrix(NA,nrow=length(combined_muts),ncol=length(tree$tip.label),dimnames = list(combined_muts,tree$tip.label))
  
  muttype1_cols=colorRampPalette(colors=c("white",cols[1]))(101)
  muttype2_cols=colorRampPalette(colors=c("white",cols[2]))(101)
  for(i in 1:nrow(hm)) {
    if(i %in% 1:length(mut_refs1)) {palette<-muttype1_cols} else {palette <-muttype2_cols}
    for(j in 1:ncol(hm)) {
      vaf<-vaf_mat[combined_muts[i],tree$tip.label[j]]
      hm[i,j]<-palette[1+round(100*vaf)]
    }
  }
  return(hm)
}

#Altered version of the MutationalPatterns signature to view just a single context
plot_192_profile_single_context=function (mut_matrix, select_context="C>T",col = "red", text_col="white",ymax = 0.3, condensed = FALSE) {
  freq <- full_context <- substitution <- context <- strand <- full_context_strand <- NULL
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context_strand") %>% 
    tidyr::separate(full_context_strand, into = c("full_context", 
                                                  "strand"), sep = "-") %>% dplyr::mutate(substitution = stringr::str_replace(full_context, 
                                                                                                                              "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
                                                                                                                                                                                         "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context, -strand), 
                        names_to = "sample", values_to = "freq") %>%
    dplyr::mutate(sample = factor(sample,levels = unique(sample)))%>%
    mutate(strand=stringr::str_to_sentence(strand))
  width <- 0.6
  spacing <- 0.5

  plot <- ggplot(data = tb%>%filter(substitution==select_context), aes(x = context, y = freq, 
                                  fill = substitution, width = width, alpha = strand)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_alpha_discrete(range = c(0.1, 1)) +
    scale_fill_manual(values = col) +
    facet_grid( ~ substitution) +
    labs(y="Relative contribution",x="Sequence context",alpha="Strand") +
    coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0,ymax, 0.1))+
    guides(fill = "none") + theme_bw() + 
      theme(axis.title.y = element_text(size = 12, vjust = 1), 
            axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(size = 5, angle = 90,vjust = 0.5),
            strip.text.x = element_text(size = 10,colour = text_col),
            panel.grid.major.x = element_blank(),
            strip.background = element_rect(fill=col),
            legend.title = element_text(size=9),
            panel.spacing.x = unit(spacing, "lines"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank())

  return(plot)
}

plot_lesion_segregation_MOD=function (vcf, per_chrom = FALSE, sample_name = NA, min_muts_mean = 10, 
          chromosomes = NA, subsample = NA) {
  vcf_list <- vcf
  max_pos <- start_mb <- notused <- n <- y <- NULL
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }
  GenomeInfoDb::genome(vcf_list) <- NA
  GenomeInfoDb::seqlevelsStyle(vcf_list) <- "NCBI"

  gr <- MutationalPatterns:::.get_strandedness_gr(vcf_list)
  tb <- MutationalPatterns:::.get_strandedness_tb(gr)
  tb$sample <- "MySample"
  sample_names <- "MySample"
  nr_samples <- 1
  max_nr_mutations <- length(vcf_list)  

  tb <- dplyr::mutate(tb, sample = factor(sample, levels = unique(sample)))
  chromosomes <- GenomeInfoDb::seqlevelsInUse(vcf_list)
  if (!length(chromosomes)) {
    chromosomes <- GenomeInfoDb::seqlevels(vcf_list)
  }
  tb$seqnames <- factor(tb$seqnames, levels = chromosomes)

  chroms <- levels(tb$seqnames)
  nr_chroms <- length(chroms)
  chrom_lengths <- GenomeInfoDb::seqlengths(vcf_list)[chroms]
  chr_starts <- tibble::tibble(start_mb = rep(1, nr_chroms), 
                               seqnames = chroms)
  chr_ends <- tibble::tibble(start_mb = chrom_lengths/1e+06, 
                             seqnames = chroms)
  chr_limits <- rbind(chr_starts, chr_ends)
  tb_limits <- tibble::tibble(sample = factor(rep(sample_names, 
                                                  each = 2 * nr_chroms), levels = levels(tb$sample)), start_mb = rep(chr_limits$start_mb, 
                                                                                                                     nr_samples), seqnames = factor(rep(chr_limits$seqnames, 
                                                                                                                                                        nr_samples), levels = chroms), y = 0.5)
  tb_mean <- tb %>% dplyr::group_by(seqnames, sample) %>% dplyr::summarize(y = mean(y), 
                                                                           n = dplyr::n(), .groups = "drop_last") %>% dplyr::ungroup() %>% 
    dplyr::right_join(tb_limits[, c("sample", "seqnames", 
                                    "start_mb")], by = c("seqnames", "sample")) %>% dplyr::filter(n >= 
                                                                                                    min_muts_mean) %>% dplyr::mutate(sample = factor(sample, 
                                                                                                                                                     levels = levels(tb$sample)))
  point_size <- 200/max_nr_mutations
  if (per_chrom == TRUE) {
    point_size <- point_size * 5
  }
  if (point_size > 2) {
    point_size <- 2
  }
  else if (point_size < 0.02) {
    point_size <- 0.02
  }
  fig <- plot_lesion_segregation_gg_MOD(tb, tb_mean, tb_limits, 
                                        point_size, sample_name)
  return(fig)
}

plot_lesion_segregation_gg_MOD=function (tb, tb_mean, tb_limits, point_size, sample_name) 
{
  y <- start_mb <- NULL
  if (is.na(sample_name)) {
    my_labs <- labs(y = "Strand", x = "Coordinate (mb)")
  }
  else {
    my_labs <- labs(y = "Strand", x = "Coordinate (mb)", 
                    title = sample_name)
  }
  fig <- ggplot(mapping = aes(y = y, x = start_mb, color = y)) + 
    geom_jitter(data = tb,width = 0.05, height = 0.1, size = point_size) +
    facet_grid( ~ seqnames, scales = "free_x", space = "free_x") +
    geom_blank(data = tb_limits) + 
    scale_y_continuous(breaks = c(0, 1), labels = c("-","+"), limits = c(-1, 2)) +
    my_labs +
    theme_bw() + 
    theme(text = element_text(size = 8), axis.text.x = element_text(size=6,angle=90), 
          axis.text.y = element_text(size = 8), 
          legend.position = "none") + 
    scale_color_gradientn(colors = c("#008BCE", "grey90", 
                                     "#E5018C"), limits = c(0, 1))
  return(fig)
}

#Functions for handling the HDP objects
create_exposures_df=function(HDP_multi,trinuc_mut_mat,key_table,minimum_branch_muts=50,sep="_") {
  sample_remove=rownames(trinuc_mut_mat)[rowSums(trinuc_mut_mat)<minimum_branch_muts]
  trinuc_mut_mat=trinuc_mut_mat[!rownames(trinuc_mut_mat)%in%sample_remove,]
  key_table=key_table[!key_table$Sample%in%sample_remove,]
  freq=nrow(trinuc_mut_mat)
  
  dp_distn <- comp_dp_distn(HDP_multi)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(trinuc_mut_mat),,drop=FALSE])
  colnames(exposures)=rownames(trinuc_mut_mat)
  rownames(exposures)<-paste0("N",rownames(exposures))
  sigs=rownames(exposures)
  sig_profiles=HDP_multi@comp_categ_distn$mean
  
  exposures_df<-as.data.frame(t(exposures),stringsAsFactors=F)%>%
    tibble::rownames_to_column("branch")%>%
    tidyr::separate(col="branch",into=c("node","Sample_ID"),sep=sep)%>%
    mutate(node=as.numeric(node))
  return(exposures_df)
}

mut_mat_HDP_comp=function(HDP_multi,ymax=0.2,plot=T){
  require(MutationalPatterns)
  sig_profiles=t(mut_example_multi@comp_categ_distn$mean)
  colnames(sig_profiles)<-paste0("N",0:(ncol(sig_profiles)-1))
  bases=c("A","C","G","T")
  subs=c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")
  rownames(sig_profiles)<-paste0(rep(rep(bases,each=4),times=6),rep(subs,each=16),rep(bases,times=24))
  if(plot){
    plot_96_profile(sig_profiles,ymax=ymax,condensed = T) 
  }
  return(sig_profiles)
}
